import requests
import time
import pandas as pd
import os

# Optional proxy, read from the environment. Set GNOMAD_PROXY_URL if the
# gnomAD API is not directly reachable from your network.
PROXY_URL = os.environ.get("GNOMAD_PROXY_URL")
proxies = {"http": PROXY_URL, "https": PROXY_URL} if PROXY_URL else None

# Define where you want to save the final downloaded data!
OUTPUT_PATH = r"...[YOUR_PATH]...\alphamissense-variant-pipeline\Assets\05_gnomad_middle_eastern_frequencies.txt"

URL = "https://gnomad.broadinstitute.org/api"

# The Full Cardiovascular Trawl Net (93 Genes)
CARDIO_GENES = [
    "SCN5A", "SCN10A", "SCN1B", "SCN2B", "SCN3B", "SCN4B", "GPD1L", "SNTA1", "RANGRF",
    "KCNQ1", "KCNH2", "KCNJ2", "KCNJ5", "KCNJ8", "KCNE1", "KCNE2", "KCNE3", "KCND3", "ABCC9",
    "CACNA1C", "CACNB2", "CACNA2D1", "RYR2", "CASQ2", "CALM1", "CALM2", "CALM3", "TRDN", "TECRL", "AKAP9", "ANK2",
    "HCN4", "MYH7", "MYBPC3", "TNNT2", "TNNI3", "TNNC1", "TPM1", "ACTC1", "MYL2", "MYL3", "MYPN", "NEXN", "ALPK3",
    "PKP2", "DSP", "DSG2", "DSC2", "JUP", "TMEM43", "TTN", "DES", "FLNC", "BAG3", "VCL", "CSRP3", "TCAP", "LDB3",
    "LMNA", "RBM20", "EMD", "PRKAG2", "LAMP2", "GLA", "TTR", "GAA", "PLN",
    "FBN1", "FBN2", "COL3A1", "MFAP5", "LOX", "TGFBR1", "TGFBR2", "TGFB2", "TGFB3", "SMAD3", "SMAD4", "SMAD6",
    "ACTA2", "MYH11", "MYLK", "PRKG1", "LDLR", "APOB", "PCSK9", "LDLRAP1", "APOE", "LPL", "APOC2", "GPIHBP1", "ABCG5", "ABCG8"
]

GRAPHQL_QUERY = """
query VariantsInGene {
  gene(gene_symbol: "TARGET_GENE", reference_genome: GRCh38) {
    variants(dataset: gnomad_r4) {
      variant_id
      exome {
        populations {
          id
          ac
          an
        }
      }
    }
  }
}
"""

print(f"🌐 Initiating Full Extraction for {len(CARDIO_GENES)} genes...")

all_variants_data = []

# Loop through all 93 genes
for i, gene in enumerate(CARDIO_GENES, 1):
    print(f"[{i}/{len(CARDIO_GENES)}] Fetching Middle Eastern data for {gene}...")
    
    exact_query = GRAPHQL_QUERY.replace("TARGET_GENE", gene)
    payload = {"query": exact_query}
    
    try:
        response = requests.post(URL, json=payload, proxies=proxies)
        
        if response.status_code == 200:
            data = response.json()
            try:
                variants = data['data']['gene']['variants']
                for variant in variants:
                    var_id = variant['variant_id']
                    mid_af = 0.0 
                    
                    if variant.get('exome') and variant['exome'].get('populations'):
                        for pop in variant['exome']['populations']:
                            if pop['id'] == 'mid':
                                ac = pop.get('ac', 0)
                                an = pop.get('an', 0)
                                if an > 0:
                                    mid_af = ac / an
                                break
                    
                    all_variants_data.append({
                        "Gene": gene,
                        "Variant_ID": var_id,
                        "Mid_Allele_Freq": mid_af
                    })
            except (TypeError, KeyError):
                print(f"  ⚠️ No variants found for {gene}.")
        else:
            print(f"  ❌ API Error for {gene}: {response.status_code}")
            
    except Exception as e:
        print(f"  ❌ Connection error on {gene}: {e}")
    
    # Polite pause to avoid IP ban (2 seconds * 93 genes = ~3 minutes total)
    time.sleep(2)

if all_variants_data:
    # Convert to DataFrame
    df_mid = pd.DataFrame(all_variants_data)
    
    # SAVE THE FILE TO YOUR HARD DRIVE
    df_mid.to_csv(OUTPUT_PATH, sep='\t', index=False)

    print("\n✅ PHASE 2 COMPLETE!")
    print(f"Successfully downloaded and saved {len(df_mid):,} Middle Eastern frequencies.")
    print(f"File saved permanently to: {OUTPUT_PATH}")
else:
    print("\n⚠️ Extraction failed. No data to save.")