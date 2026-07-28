import pandas as pd
import time

# Our new starting point and the destination for the gene list
INPUT_PATH = r"...[YOUR_PATH]...\alphamissense-variant-pipeline\Assets\01_clinvar_vus_extracted.txt"
OUTPUT_PATH = r"...[YOUR_PATH]...\Research\AI Genetics Project\alphamissense-variant-pipeline\Assets\02_unique_vus_genes.txt"

print("🔍 Scanning for unique genes...")
start_time = time.time()

# We only load the 'GeneSymbol' column. This makes the process lightning fast.
df = pd.read_csv(INPUT_PATH, sep='\t', usecols=['GeneSymbol'], low_memory=False)

# Extract unique genes, drop any blank rows, and sort them alphabetically
unique_genes = df['GeneSymbol'].dropna().unique()
unique_genes.sort()

# Write them to a simple text file, one gene per line
with open(OUTPUT_PATH, 'w') as file:
    for gene in unique_genes:
        file.write(f"{gene}\n")

end_time = time.time()

print(f"✅ Done in {(end_time - start_time):.2f} seconds!")
print(f"Found {len(unique_genes):,} unique genes with VUSs.")
print(f"Saved to: {OUTPUT_PATH}")