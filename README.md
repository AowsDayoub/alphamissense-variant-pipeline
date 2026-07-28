# Cardiovascular VUS Reclassification Pipeline

A ten-stage Python pipeline that narrows 2.3 million ClinVar variants of uncertain
significance down to 306 cardiovascular candidates for reclassification, by combining
population allele frequencies from gnomAD with AlphaMissense pathogenicity predictions
across 71.7 million scored variants.

Built in April 2026 on a consumer laptop.

## The question

Variants of uncertain significance are the dead end of clinical genetics. A patient's
cardiac gene panel returns a variant that has never been functionally characterised,
carries no population frequency data for their ancestry group, and cannot be acted on.

Middle Eastern populations are severely underrepresented in the reference databases
that drive this classification, which means variants circulating in these populations
are disproportionately likely to remain uncertain.

This pipeline asks a narrower, tractable version of that problem: among cardiovascular
VUS in ClinVar, which are absent or ultra-rare in Middle Eastern populations, and of
those, which does a structural pathogenicity predictor flag as likely pathogenic?
Those are the variants worth taking to a wet lab.

## Results

| Stage | Input | Output |
|---|---|---|
| ClinVar extraction | 8,907,730 rows | 2,299,605 VUS (GRCh38) |
| Gene inventory | 2,299,605 VUS | 21,157 unique genes |
| Cardiovascular panel filter | 2,299,605 VUS | 85,340 cardiovascular VUS |
| gnomAD retrieval | 93 genes | 551,501 variants |
| Column trimming | 37 columns | 13 columns |
| ClinVar × gnomAD merge | 85,340 + 551,501 | 85,340 annotated |
| ACMG stratification | 85,340 | 336 cleared, 85,004 retained |
| AlphaMissense scan | 71,697,556 predictions | 70,087 targets scored |
| **Final candidate set** | | **306 variants** |

### Stratification detail

Of 85,340 cardiovascular VUS, 2,028 had a recorded Middle Eastern allele count and
83,312 did not.

ACMG frequency criteria cleared 336 variants: one exceeding 5% global frequency (BA1)
and 335 exceeding the disease-specific threshold for their gene (BS1). The remaining
85,004 advanced as targets, split into 1,692 borderline carriers and 83,312 variants
absent from Middle Eastern data.

### AlphaMissense yield

Of 85,004 targets, 70,087 received a pathogenicity score. The 14,917 without a score
are predominantly frameshift and other non-missense variants outside the model's scope.

| | Absent from ME data | Borderline carriers |
|---|---|---|
| Scored | 68,612 | 1,475 |
| Likely pathogenic | 22,715 | **306** |
| Ambiguous | 8,337 | 192 |
| Likely benign | 37,560 | 977 |
| Mean pathogenicity | 0.4206 | 0.3222 |

Variants absent from Middle Eastern data score higher on average than those circulating
at low frequency (0.4206 vs 0.3222) — the direction predicted if the stratification is
capturing genuine selection against damaging alleles rather than sorting noise.

The 306 borderline carriers classified as likely pathogenic form the final candidate
set: variants present but ultra-rare in Middle Eastern populations, where predicted
structural damage and epidemiological rarity agree.

## Approach

**Stage 1 — ClinVar extraction and narrowing** (scripts 01–03)

`variant_summary.txt` is read in 500,000-row chunks and filtered to GRCh38 assembly
rows classified as "Uncertain significance". A gene inventory is taken across all
retained VUS, then the set is narrowed to a 93-gene cardiovascular panel spanning four
disease categories: channelopathies, cardiomyopathies and structural disease,
aortopathies and connective tissue disorders, and inherited lipid disorders.

**Stage 2 — gnomAD population frequencies** (scripts 04–07)

Allele frequency data is pulled from the gnomAD v4 GraphQL API gene by gene, with a
6.5-second delay between requests to stay inside the rate limit. The query retrieves
variant identifiers, HGVS notation, consequence, global allele counts and frequencies,
homozygote counts, and a full per-population breakdown.

A pilot script validates the API contract on three genes before the production run.

The resulting master dataset is trimmed from 37 columns to 13 — core annotation,
global frequencies, and Middle Eastern frequencies — with missing Middle Eastern counts
normalised to zero rather than dropped, so that global frequency data remains available
for every variant downstream.

**Stage 3 — Integration** (script 08)

ClinVar and gnomAD use different variant identifier conventions. ClinVar coordinates are
reassembled into gnomAD's `chrom-pos-ref-alt` format and the two datasets are joined on
that key with a left join, preserving every cardiovascular VUS.

**Stage 4 — ACMG epidemiological stratification** (script 09)

Each variant is assigned a maximum credible population allele frequency based on the
prevalence of the disease its gene causes: 0.1% for channelopathy and cardiomyopathy
genes, 0.01% for aortopathy genes, 0.5% for lipid disorder genes.

Three ACMG criteria are then applied:

- **BA1** — global allele frequency above 5%. Stand-alone benign; discarded.
- **BS1** — Middle Eastern frequency above the gene's disease-specific threshold. Too
  common in this population to cause the disease; discarded.
- **PM2** — absent from Middle Eastern data. Retained as a candidate.

Variants present but below their gene's threshold fall between BS1 and PM2 and are
retained as borderline carriers. Both retained groups pass to the final stage.

**Stage 5 — AlphaMissense interrogation** (script 10)

The full AlphaMissense hg38 release — 71,697,556 predicted pathogenicity scores — is
streamed in 5,000,000-row chunks with only six columns parsed. Each chunk is matched
against the target set by reconstructed variant identifier.

Variants scored across multiple transcripts are ranked by pathogenicity and pivoted into
wide format, so every transcript-level prediction is preserved rather than collapsed to
a single value. Scored predictions are then merged back onto the ACMG classifications.

## Performance

Runtimes on the hardware below:

| Script | Runtime |
|---|---|
| 01 — ClinVar extraction | 2.30 min |
| 02 — Gene inventory | 9.32 s |
| 03 — Cardiovascular filter | 34.91 s |
| 05 — gnomAD ME frequencies | ~3 min (API-bound) |
| 06 — gnomAD global master | ~10 min (API-bound) |
| 07 — Column trimming | 9.18 s |
| 08 — ClinVar × gnomAD merge | 6.31 s |
| 09 — ACMG stratification | 4.44 s |
| 10 — AlphaMissense scan | 2.96 min |

The gnomAD stages are bound not by compute but by the API rate limit: 93 sequential
requests at 6.5 seconds apart.

| Hardware | |
|---|---|
| Machine | HP Victus laptop |
| CPU | Intel Core i5, 13th generation |
| Memory | 32 GB |
| Storage | 2 TB NVMe SSD |
| Compute | CPU only |

Chunk sizes reflect the constraint directly. The ClinVar stages run at 500,000 rows;
by the AlphaMissense stage this had risen tenfold to 5,000,000, once memory behaviour
under load was understood.

## Repository structure

```
alphamissense-variant-pipeline/
├── README.md
├── LICENSE
├── .gitignore
├── pipeline/
│   ├── 01_clinvar_vus_extractor.py           Chunked GRCh38 + VUS extraction
│   ├── 02_clinvar_unique_gene_extractor.py   Gene inventory across all VUS
│   ├── 03_clinvar_vus_cardiology_filter.py   93-gene cardiovascular panel filter
│   ├── 04_gnomad_pilot_batch.py              API validation, 3 genes
│   ├── 05_gnomad_middle_eastern_variants.py  ME frequencies, full panel
│   ├── 06_gnomad_global_master.py            Full population breakdown, full panel
│   ├── 07_gnomad_column_trimmer.py           Column reduction and normalisation
│   ├── 08_clinvar_gnomad_alignment.py        Identifier reconciliation and join
│   ├── 09_acmg_epidemiology_filter.py        BA1 / BS1 / PM2 stratification
│   └── 10_alphamissense_interrogation.py     71.7M-row scan and transcript pivot
└── assets/
    ├── 00_cardiovascular_gene_panel.txt      The 93-gene panel
    ├── sample_variant_summary.txt            Format sample, ClinVar input
    ├── sample_alphamissense_hg38.txt         Format sample, AlphaMissense input
    └── 10_final_candidate_variants.txt       Final candidate set, 306 variants
```

Intermediate outputs are written as tab-separated text files, each stage consuming the
previous stage's output. These are not included here — several exceed 100 MB — but are
fully reproducible by running the pipeline in order.

Script 02 produces a gene inventory used to survey the VUS landscape before selecting
the cardiovascular panel. Its output does not feed later stages.

## Requirements

```
python >= 3.9
pandas
requests
```

```bash
pip install pandas requests
```

## Running the pipeline

Each script defines its input and output paths at the top of the file, written as
`...[YOUR_PATH]...\alphamissense-variant-pipeline\Assets\`. Replace this with your own
directory before running.

Download the two source datasets into that directory first (see **Data sources** below):

- `00_source_variant_summary.txt` — ClinVar
- `00_source_alphamissense_hg38.tsv` — AlphaMissense

If the gnomAD API is not directly reachable from your network, scripts 04–06 read an
optional proxy from the environment:

```bash
export GNOMAD_PROXY_URL="http://user:password@host:port"
```

Then run in order:

```bash
python pipeline/01_clinvar_vus_extractor.py
python pipeline/02_clinvar_unique_gene_extractor.py
python pipeline/03_clinvar_vus_cardiology_filter.py
python pipeline/04_gnomad_pilot_batch.py
python pipeline/05_gnomad_middle_eastern_variants.py
python pipeline/06_gnomad_global_master.py
python pipeline/07_gnomad_column_trimmer.py
python pipeline/08_clinvar_gnomad_alignment.py
python pipeline/09_acmg_epidemiology_filter.py
python pipeline/10_alphamissense_interrogation.py
```

Script 04 is a pilot and can be skipped; it is retained to document the API validation
step.

## Data sources

The source datasets are not redistributed here. Sample files showing the expected input
format are provided in `assets/`.

**ClinVar** — `variant_summary.txt`, National Center for Biotechnology Information.
https://www.ncbi.nlm.nih.gov/clinvar/

**gnomAD v4** — Genome Aggregation Database, Broad Institute, accessed via the public
GraphQL API. https://gnomad.broadinstitute.org/

**AlphaMissense** — `AlphaMissense_hg38.tsv`, Google DeepMind. Consult the official
repository for current licensing terms before use.
https://github.com/google-deepmind/alphamissense

If you use the AlphaMissense data, cite:

> Cheng, J., Novati, G., Pan, J., Bycroft, C., Žemgulytė, A., Applebaum, T.,
> Pritzel, A., Wong, L.H., Zielinski, M., Sargeant, T., Schneider, R.G.,
> Senior, A.W., Jumper, J., Hassabis, D., Kohli, P., & Avsec, Ž. (2023).
> Accurate proteome-wide missense variant effect prediction with AlphaMissense.
> *Science*, 381(6664), eadg7492. https://doi.org/10.1126/science.adg7492

## Disclaimer

AlphaMissense predictions are computational predictions with varying confidence and have
not been validated for clinical use. The ACMG criteria applied here are limited to
population frequency evidence and do not constitute a complete variant classification.
Nothing produced by this pipeline is suitable for diagnosis or clinical decision-making.
Outputs are research hypotheses requiring functional validation.

## Licence

Code released under the MIT Licence. This covers the pipeline only, not the source
datasets, which remain subject to their own terms.

## Author

**Aows Dayoub** — MD Candidate, Faculty of Medicine, Damascus University

- Google Scholar: https://scholar.google.com/citations?user=ur9NdVwAAAAJ
- ORCID: https://orcid.org/0009-0009-7896-9130
- Email: AowsDayoub02@gmail.com
