import pandas as pd
import time

# ==========================================
# FILE PATHS
# ==========================================
CLINVAR_PATH = r"...[YOUR_PATH]...\alphamissense-variant-pipeline\Assets\03_clinvar_vus_cardiology.txt"
GNOMAD_PATH = r"...[YOUR_PATH]...\alphamissense-variant-pipeline\Assets\07_gnomad_trimmed.txt"
OUTPUT_PATH = r"...[YOUR_PATH]...\alphamissense-variant-pipeline\Assets\08_clinvar_gnomad_merged.txt"

print("⚔️ ROUND 1: Loading databases into memory...")
start_time = time.time()

# Load the datasets (low_memory=False handles mixed data types safely)
df_clinvar = pd.read_csv(CLINVAR_PATH, sep='\t', low_memory=False)
df_gnomad = pd.read_csv(GNOMAD_PATH, sep='\t', low_memory=False)

print(f"   ClinVar Suspects loaded: {len(df_clinvar):,}")
print(f"   gnomAD reference variants loaded: {len(df_gnomad):,}")

# ==========================================
# STEP 1: FORGE THE PRIMARY KEY IN CLINVAR
# ==========================================
print("\n🔗 Forging Primary Keys for matching...")

# Ensure all pieces are treated as strings so we can stitch them together
df_clinvar['Chromosome'] = df_clinvar['Chromosome'].astype(str)
df_clinvar['PositionVCF'] = df_clinvar['PositionVCF'].astype(str)
df_clinvar['ReferenceAlleleVCF'] = df_clinvar['ReferenceAlleleVCF'].astype(str)
df_clinvar['AlternateAlleleVCF'] = df_clinvar['AlternateAlleleVCF'].astype(str)

# Construct the exact gnomAD format: "Chrom-Pos-Ref-Alt"
df_clinvar['Variant_ID'] = (
    df_clinvar['Chromosome'] + "-" + 
    df_clinvar['PositionVCF'] + "-" + 
    df_clinvar['ReferenceAlleleVCF'] + "-" + 
    df_clinvar['AlternateAlleleVCF']
)

# ==========================================
# STEP 2: THE MERGE (MORTAL KOMBAT)
# ==========================================
print("💥 Executing the Merge...")

# We do a 'left' merge to keep ALL ClinVar suspects, attaching gnomAD data where available
df_final = pd.merge(
    df_clinvar, 
    df_gnomad, 
    on='Variant_ID',
    how='left'
)

# If a ClinVar variant wasn't in our Middle Eastern list, its frequency is zero.
df_final['mid_AF'] = df_final['mid_AF'].fillna(0.0)
df_final['mid_AC'] = df_final['mid_AC'].fillna(0)

# ==========================================
# STEP 3: SORT AND SAVE
# ==========================================
print("🧹 Cleaning up and saving the Master Database...")

# Sort by highest Middle Eastern frequency first
df_final = df_final.sort_values(by='mid_AF', ascending=False)

# Save the ultimate file
df_final.to_csv(OUTPUT_PATH, sep='\t', index=False)

end_time = time.time()

# Quick calculation to see how many suspects actually had an alibi
alibis_found = len(df_final[df_final['mid_AF'] > 0])

print("\n✅ MORTAL KOMBAT COMPLETE!")
print(f"Time taken: {(end_time - start_time):.2f} seconds.")
print("-" * 50)
print(f"Total Cardiovascular VUSs Assessed: {len(df_final):,}")
print(f"Suspects with Middle Eastern Alibis (AF > 0): {alibis_found:,}")
print(f"Suspects with NO Middle Eastern footprint (AF = 0): {(len(df_final) - alibis_found):,}")
print(f"💾 The Ultimate Master File is ready: {OUTPUT_PATH}")