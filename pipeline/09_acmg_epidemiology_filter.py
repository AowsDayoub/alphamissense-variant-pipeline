import pandas as pd
import time

# ==========================================
# FILE PATHS
# ==========================================
INPUT_PATH = r"...[YOUR_PATH]...\alphamissense-variant-pipeline\Assets\08_clinvar_gnomad_merged.txt"

# We are going to split the data into two files: The "Trash" and the "Hit List"
ALIBI_OUTPUT = r"...[YOUR_PATH]...\alphamissense-variant-pipeline\Assets\09a_acmg_discarded_benign.txt"
HITLIST_OUTPUT = r"...[YOUR_PATH]...\alphamissense-variant-pipeline\Assets\09b_acmg_target_hitlist.txt"

print("🧬 Initiating ACMG Epidemiological Stratification...")
start_time = time.time()

# Load the Master Database
df = pd.read_csv(INPUT_PATH, sep='\t', low_memory=False)

# Ensure our frequency columns are treated as math-ready decimals
df['Global_AF'] = df['Global_AF'].fillna(0.0).astype(float)
df['mid_AF'] = df['mid_AF'].fillna(0.0).astype(float)

# ==========================================
# STEP 1: DEFINE THE DISEASE CATEGORIES
# ==========================================
ELECTRIC_GENES = ["SCN5A", "SCN10A", "SCN1B", "SCN2B", "SCN3B", "SCN4B", "GPD1L", "SNTA1", "RANGRF", "KCNQ1", "KCNH2", "KCNJ2", "KCNJ5", "KCNJ8", "KCNE1", "KCNE2", "KCNE3", "KCND3", "ABCC9", "CACNA1C", "CACNB2", "CACNA2D1", "RYR2", "CASQ2", "CALM1", "CALM2", "CALM3", "TRDN", "TECRL", "AKAP9", "ANK2", "HCN4"]
STRUCTURAL_GENES = ["MYH7", "MYBPC3", "TNNT2", "TNNI3", "TNNC1", "TPM1", "ACTC1", "MYL2", "MYL3", "MYPN", "NEXN", "ALPK3", "PKP2", "DSP", "DSG2", "DSC2", "JUP", "TMEM43", "TTN", "DES", "FLNC", "BAG3", "VCL", "CSRP3", "TCAP", "LDB3", "LMNA", "RBM20", "EMD", "PRKAG2", "LAMP2", "GLA", "TTR", "GAA", "PLN"]
VASCULAR_GENES = ["FBN1", "FBN2", "COL3A1", "MFAP5", "LOX", "TGFBR1", "TGFBR2", "TGFB2", "TGFB3", "SMAD3", "SMAD4", "SMAD6", "ACTA2", "MYH11", "MYLK", "PRKG1"]
METABOLIC_GENES = ["LDLR", "APOB", "PCSK9", "LDLRAP1", "APOE", "LPL", "APOC2", "GPIHBP1", "ABCG5", "ABCG8"]

# Map each gene to its specific ACMG BS1 Threshold (AF_max)
threshold_map = {}
for g in ELECTRIC_GENES + STRUCTURAL_GENES: threshold_map[g] = 0.001   # 0.1%
for g in VASCULAR_GENES: threshold_map[g] = 0.0001                     # 0.01%
for g in METABOLIC_GENES: threshold_map[g] = 0.005                     # 0.5%

# Attach the strict threshold to every variant row based on its Gene
df['Disease_AF_max'] = df['GeneSymbol'].map(threshold_map)

# ==========================================
# STEP 2: THE ACMG BUCKETING LOGIC
# ==========================================
print("⚖️ Applying Rules BA1, BS1, and PM2...")

def categorize_variant(row):
    global_af = row['Global_AF']
    mid_af = row['mid_AF']
    af_max = row['Disease_AF_max']
    
    # RULE BA1: Ironclad Global Alibi (>5% anywhere in the world)
    if global_af > 0.05:
        return 'Bucket 1: BA1 (Global Benign)'
        
    # RULE BS1: Strong Middle Eastern Alibi (> Disease AF_max)
    if mid_af > af_max:
        return 'Bucket 1: BS1 (Levantine Benign)'
        
    # THE SUSPECT ZONE
    if mid_af > 0:
        return 'Bucket 2: Borderline Carrier'
    
    # RULE PM2: Completely absent from the Middle East
    return 'Bucket 3: PM2 (Ghost Variant)'

# Apply the logic engine to every row
df['ACMG_Classification'] = df.apply(categorize_variant, axis=1)

# ==========================================
# STEP 3: SPLIT AND SAVE THE TARGETS
# ==========================================
print("📂 Splitting data into Alibis and Targets...")

# Separate the innocent bystanders from the highly suspicious suspects
df_alibis = df[df['ACMG_Classification'].str.startswith('Bucket 1')]
df_targets = df[df['ACMG_Classification'].str.startswith('Bucket 2') | df['ACMG_Classification'].str.startswith('Bucket 3')]

# Save them locally
df_alibis.to_csv(ALIBI_OUTPUT, sep='\t', index=False)
df_targets.to_csv(HITLIST_OUTPUT, sep='\t', index=False)

end_time = time.time()

# Tally up the final results for the console
ba1_count = len(df[df['ACMG_Classification'] == 'Bucket 1: BA1 (Global Benign)'])
bs1_count = len(df[df['ACMG_Classification'] == 'Bucket 1: BS1 (Levantine Benign)'])
borderline_count = len(df[df['ACMG_Classification'] == 'Bucket 2: Borderline Carrier'])
ghost_count = len(df[df['ACMG_Classification'] == 'Bucket 3: PM2 (Ghost Variant)'])

print("\n✅ STRATIFICATION COMPLETE!")
print(f"Time taken: {(end_time - start_time):.2f} seconds.")
print("-" * 50)
print("🗑️ THE ALIBIS (Discarded from study):")
print(f"  - BA1 (Globally >5%): {ba1_count:,}")
print(f"  - BS1 (Too common for disease): {bs1_count:,}")
print(f"    Total Cleared: {(ba1_count + bs1_count):,}")
print("\n🎯 THE HIT LIST (Heading to AlphaMissense):")
print(f"  - Borderline Carriers (Ultra-rare in ME): {borderline_count:,}")
print(f"  - PM2 Ghosts (Absent from ME): {ghost_count:,}")
print(f"    Total Targets: {(borderline_count + ghost_count):,}")
print("-" * 50)
print(f"💾 Target List saved to: {HITLIST_OUTPUT}")