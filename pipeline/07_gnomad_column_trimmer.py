import pandas as pd
import time

# Define your input (the 140MB Master File) and your new Output file
INPUT_PATH = r"...[YOUR_PATH]...\alphamissense-variant-pipeline\Assets\06_gnomad_global_master.txt"
OUTPUT_PATH = r"...[YOUR_PATH]...\alphamissense-variant-pipeline\Assets\07_gnomad_trimmed.txt"

print("🔪 Initiating Data Trimming Protocol...")
start_time = time.time()

# Load the massive 140MB file into memory
# low_memory=False prevents pandas from throwing warnings about mixed data types
df = pd.read_csv(INPUT_PATH, sep='\t', low_memory=False)

original_row_count = len(df)
original_col_count = len(df.columns)

# ==========================================
# STEP 1: VERTICAL FILTRATION (Columns)
# ==========================================
# We define the core columns we absolutely want to keep
core_columns = ['Gene', 'Variant_ID', 'RSID', 'HGVSc', 'HGVSp', 'Consequence']

# We use a "List Comprehension" to dynamically find the columns we want to keep
columns_to_keep = [
    col for col in df.columns 
    if col in core_columns 
    or col.startswith('Global_') 
    or col.startswith('mid_')
]

# Overwrite the dataframe with ONLY the columns we kept
df = df[columns_to_keep].copy()

# ==========================================
# STEP 2: NORMALISE MIDDLE EASTERN COUNTS
# ==========================================
# A variant with no Middle Eastern entry is recorded as zero rather than missing,
# so that global frequencies stay available for every variant downstream.
if 'mid_AC' in df.columns:
    df['mid_AC'] = df['mid_AC'].fillna(0)
    df['mid_AF'] = df['mid_AF'].fillna(0.0)
else:
    print("⚠️ CRITICAL ERROR: Could not find the 'mid_AC' column in the dataset!")

# ==========================================
# STEP 3: SAVE THE LEAN DATASET
# ==========================================
# Sort the final list by Middle Eastern Frequency (Highest to Lowest)
df = df.sort_values(by='mid_AF', ascending=False)

df.to_csv(OUTPUT_PATH, sep='\t', index=False)

end_time = time.time()

print("\n✅ TRIMMING COMPLETE!")
print(f"Time taken: {(end_time - start_time):.2f} seconds.")
print("-" * 40)
print(f"📉 Columns Reduced: From {original_col_count} down to {len(df.columns)}.")
print(f"📊 Rows retained: {len(df):,} of {original_row_count:,} (Middle Eastern counts normalised).")
print(f"💾 File saved to: {OUTPUT_PATH}")