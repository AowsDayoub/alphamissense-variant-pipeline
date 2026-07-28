import requests
import time
import pandas as pd
import os

# Optional proxy, read from the environment. Set GNOMAD_PROXY_URL if the
# gnomAD API is not directly reachable from your network.
PROXY_URL = os.environ.get("GNOMAD_PROXY_URL")
proxies = {"http": PROXY_URL, "https": PROXY_URL} if PROXY_URL else None

# The destination for your massive dataset
OUTPUT_PATH = r"...[YOUR_PATH]...\alphamissense-variant-pipeline\Assets\06_gnomad_global_master.txt"
URL = "https://gnomad.broadinstitute.org/api"

# The Full 93-Gene Cardiovascular Master Panel
CARDIO_GENES = [
    # Electric Heart
    "SCN5A", "SCN10A", "SCN1B", "SCN2B", "SCN3B", "SCN4B", "GPD1L", "SNTA1", "RANGRF",
    "KCNQ1", "KCNH2", "KCNJ2", "KCNJ5", "KCNJ8", "KCNE1", "KCNE2", "KCNE3", "KCND3", "ABCC9",
    "CACNA1C", "CACNB2", "CACNA2D1", "RYR2", "CASQ2", "CALM1", "CALM2", "CALM3", "TRDN", "TECRL", "AKAP9", "ANK2",
    "HCN4", 
    # Structural Heart
    "MYH7", "MYBPC3", "TNNT2", "TNNI3", "TNNC1", "TPM1", "ACTC1", "MYL2", "MYL3", "MYPN", "NEXN", "ALPK3",
    "PKP2", "DSP", "DSG2", "DSC2", "JUP", "TMEM43", 
    "TTN", "DES", "FLNC", "BAG3", "VCL", "CSRP3", "TCAP", "LDB3",
    "LMNA", "RBM20", "EMD", 
    "PRKAG2", "LAMP2", "GLA", "TTR", "GAA", 
    "PLN",
    # Vascular Heart
    "FBN1", "FBN2", "COL3A1", "MFAP5", "LOX", 
    "TGFBR1", "TGFBR2", "TGFB2", "TGFB3", "SMAD3", "SMAD4", "SMAD6",
    "ACTA2", "MYH11", "MYLK", "PRKG1", 
    # Metabolic Heart
    "LDLR", "APOB", "PCSK9", "LDLRAP1", "APOE", "LPL", "APOC2", "GPIHBP1", "ABCG5", "ABCG8"
]

# The GREEDY GraphQL Query
GRAPHQL_QUERY = """
query VariantsInGene {
  gene(gene_symbol: "TARGET_GENE", reference_genome: GRCh38) {
    variants(dataset: gnomad_r4) {
      variant_id
      rsids
      hgvsc
      hgvsp
      consequence
      exome {
        ac
        an
        af
        ac_hom
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

print(f"🌐 Initiating Full OMNIVORE Extraction for {len(CARDIO_GENES)} genes...")
start_time = time.time()

all_variants_data = []

for i, gene in enumerate(CARDIO_GENES, 1):
    print(f"[{i}/{len(CARDIO_GENES)}] Downloading full variant database for {gene}...")
    
    exact_query = GRAPHQL_QUERY.replace("TARGET_GENE", gene)
    payload = {"query": exact_query}
    
    try:
        response = requests.post(URL, json=payload, proxies=proxies)
        
        if response.status_code == 200:
            data = response.json()
            try:
                variants = data['data']['gene']['variants']
                for variant in variants:
                    
                    # 1. Base Clinical Data (Raw strings straight from Broad)
                    v_dict = {
                        "Gene": gene,
                        "Variant_ID": variant.get('variant_id'),
                        "RSID": ", ".join(variant.get('rsids', [])), 
                        "HGVSc": variant.get('hgvsc'),
                        "HGVSp": variant.get('hgvsp'),
                        "Consequence": variant.get('consequence')
                    }
                    
                    exome = variant.get('exome')
                    if exome:
                        # 2. Global Totals (Raw counts)
                        v_dict['Global_AC'] = exome.get('ac', 0)
                        v_dict['Global_AN'] = exome.get('an', 0)
                        v_dict['Global_AF'] = exome.get('af', 0.0)
                        v_dict['Global_Homozygotes'] = exome.get('ac_hom', 0)
                        
                        # 3. Dynamic Population Unpacking
                        if exome.get('populations'):
                            for pop in exome['populations']:
                                pop_id = pop['id']
                                ac = pop.get('ac', 0)
                                an = pop.get('an', 0)
                                
                                # Save raw numbers directly
                                v_dict[f'{pop_id}_AC'] = ac
                                v_dict[f'{pop_id}_AN'] = an
                                
                                # Local calculation for Frequency
                                if an > 0:
                                    v_dict[f'{pop_id}_AF'] = ac / an
                                else:
                                    v_dict[f'{pop_id}_AF'] = 0.0
                                    
                    all_variants_data.append(v_dict)
                    
            except (TypeError, KeyError):
                print(f"  ⚠️ No variants found for {gene}.")
        else:
            print(f"  ❌ API Error for {gene}: {response.status_code}")
            print(f"  Server Response: {response.text}")
            
    except Exception as e:
        print(f"  ❌ Connection error on {gene}: {e}")
    
    # Strictly obeying the 10-requests-per-minute speed limit
    time.sleep(6.5)

end_time = time.time()
total_minutes = (end_time - start_time) / 60

if all_variants_data:
    df_master = pd.DataFrame(all_variants_data)
    df_master.to_csv(OUTPUT_PATH, sep='\t', index=False)

    print("\n✅ PHASE 2 COMPLETE: MASTER DATASET ACQUIRED!")
    print(f"Successfully downloaded {len(df_master):,} completely detailed variants.")
    print(f"Total Columns extracted per variant: {len(df_master.columns)}")
    print(f"File saved permanently to: {OUTPUT_PATH}")
    print(f"⏱️ Total extraction time: {total_minutes:.2f} minutes.")
else:
    print("\n⚠️ Extraction failed. No data to save.")