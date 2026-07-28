import requests
import time
import pandas as pd
import os

# Optional proxy, read from the environment. Set GNOMAD_PROXY_URL if the
# gnomAD API is not directly reachable from your network.
PROXY_URL = os.environ.get("GNOMAD_PROXY_URL")
proxies = {"http": PROXY_URL, "https": PROXY_URL} if PROXY_URL else None


# The official gnomAD API endpoint
URL = "https://gnomad.broadinstitute.org/api"

# Pilot Batch: 3 High-Yield Cardiology Genes
PILOT_GENES = ["SCN5A", "MYH7", "KCNQ1"]

# The CORRECTED GraphQL Query requesting 'ac' and 'an'
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

print("🌐 Establishing secure connection to the Broad Institute...")

all_variants_data = []

# Loop through our pilot genes
for gene in PILOT_GENES:
    print(f"Fetching Middle Eastern data for {gene}...")
    
    # Swap our placeholder for the actual gene name
    exact_query = GRAPHQL_QUERY.replace("TARGET_GENE", gene)
    payload = {"query": exact_query}
    
    try:
        response = requests.post(URL, json=payload, proxies=proxies)
        
        if response.status_code == 200:
            data = response.json()
            
            # Dig into the nested JSON response
            try:
                variants = data['data']['gene']['variants']
                
                for variant in variants:
                    var_id = variant['variant_id']
                    mid_af = 0.0 # Default frequency is 0
                    
                    # Check if exome data exists for this variant
                    if variant.get('exome') and variant['exome'].get('populations'):
                        # Search specifically for the 'mid' (Middle Eastern) ID
                        for pop in variant['exome']['populations']:
                            if pop['id'] == 'mid':
                                ac = pop.get('ac', 0)
                                an = pop.get('an', 0)
                                
                                # Calculate the Allele Frequency (prevent division by zero)
                                if an > 0:
                                    mid_af = ac / an
                                break
                    
                    # Save the variant and its Middle Eastern frequency
                    all_variants_data.append({
                        "Gene": gene,
                        "Variant_ID": var_id,
                        "Mid_Allele_Freq": mid_af
                    })
                    
            except (TypeError, KeyError):
                print(f"⚠️ No data returned for {gene}. Check gene name or response format.")
        else:
            print(f"❌ API Error for {gene}: {response.status_code}")
            print(f"Server Response: {response.text}") 
            
    except Exception as e:
        print(f"❌ An error occurred: {e}")
    
    # The Polite Pause
    time.sleep(2)

# If we successfully gathered data, build the dataframe
if all_variants_data:
    df_mid = pd.DataFrame(all_variants_data)
    df_mid_sorted = df_mid.sort_values(by="Mid_Allele_Freq", ascending=False)

    print("\n✅ API Pilot Complete!")
    print(f"Successfully downloaded {len(df_mid):,} total variants for the 3 pilot genes.")
    print("\nTop 10 variants with the highest Middle Eastern presence (The 'Alibis'):")
    print(df_mid_sorted.head(10).to_string(index=False))
else:
    print("\n⚠️ No data was extracted. Please review the errors above.")