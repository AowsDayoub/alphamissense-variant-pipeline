import pandas as pd
import time

# ==========================================
# FILE PATHS
# ==========================================
HITLIST_PATH = r"...[YOUR_PATH]...\alphamissense-variant-pipeline\Assets\09b_acmg_target_hitlist.txt"
ALPHAMISSENSE_PATH = r"...[YOUR_PATH]...\alphamissense-variant-pipeline\Assets\00_source_alphamissense_hg38.tsv"
FINAL_SOLVED_PATH = r"...[YOUR_PATH]...\alphamissense-variant-pipeline\Assets\10_final_candidate_variants.txt"

CHUNK_SIZE = 5_000_000 

print("🧠 Initiating AlphaMissense Biophysical Interrogation...")
start_time = time.time()

# 1. Load our epidemiological targets
df_targets = pd.read_csv(HITLIST_PATH, sep='\t', low_memory=False)
target_ids = set(df_targets['Variant_ID'].unique())
print(f"   Loaded {len(target_ids):,} unique targets to hunt for.")

print(f"\n📥 Scanning DeepMind's 71.7 Million Predictions (Chunk Size: {CHUNK_SIZE:,})...")

# ==========================================
# STEP 1: CHUNKING & FILTERING
# ==========================================
filtered_chunks = []
total_processed = 0

columns_to_read = ['#CHROM', 'POS', 'REF', 'ALT', 'am_pathogenicity', 'am_class']

for chunk in pd.read_csv(ALPHAMISSENSE_PATH, sep='\t', skiprows=3, usecols=columns_to_read, chunksize=CHUNK_SIZE, low_memory=False):
    
    chunk['#CHROM'] = chunk['#CHROM'].astype(str).str.replace('chr', '', regex=False)
    chunk['POS'] = chunk['POS'].astype(str).str.replace(r'\.0$', '', regex=True)
    chunk['REF'] = chunk['REF'].astype(str)
    chunk['ALT'] = chunk['ALT'].astype(str)

    chunk['Variant_ID'] = chunk['#CHROM'] + "-" + chunk['POS'] + "-" + chunk['REF'] + "-" + chunk['ALT']

    matched_chunk = chunk[chunk['Variant_ID'].isin(target_ids)]

    if not matched_chunk.empty:
        filtered_chunks.append(matched_chunk[['Variant_ID', 'am_pathogenicity', 'am_class']])

    total_processed += len(chunk)
    print(f"   ... Scanned {total_processed:,} rows ...")

# ==========================================
# STEP 2: DYNAMIC WIDE-COLUMN PIVOT
# ==========================================
print("\n🧮 Processing Transcripts and Forging Wide Columns...")

if filtered_chunks:
    df_am_matches = pd.concat(filtered_chunks, ignore_index=True)
    df_am_matches['am_pathogenicity'] = pd.to_numeric(df_am_matches['am_pathogenicity'], errors='coerce')

    df_am_matches = df_am_matches.sort_values(by=['Variant_ID', 'am_pathogenicity'], ascending=[True, False])
    df_am_matches['transcript_rank'] = df_am_matches.groupby('Variant_ID').cumcount() + 1

    df_wide = df_am_matches.pivot(index='Variant_ID', columns='transcript_rank', values=['am_pathogenicity', 'am_class'])
    df_wide.columns = [f"{col[0]}_{col[1]}" if col[1] > 1 else col[0] for col in df_wide.columns]
    df_wide = df_wide.reset_index()
else:
    # Fallback in the impossible event of zero matches
    df_wide = pd.DataFrame(columns=['Variant_ID', 'am_pathogenicity', 'am_class'])

# ==========================================
# STEP 3: THE FINAL MERGE
# ==========================================
print("💥 Executing the Final Merge...")

df_final = pd.merge(
    df_targets, 
    df_wide, 
    on='Variant_ID', 
    how='left'
)

# ==========================================
# STEP 4: DESCRIPTIVE STATISTICS (The "Results" Section)
# ==========================================
print("\n📊 CALCULATING DESCRIPTIVE STATISTICS...")

# Isolate variants that actually got a score (Missense) vs those that didn't (Non-missense/Frameshifts)
df_scored = df_final.dropna(subset=['am_pathogenicity'])
non_missense_count = len(df_final) - len(df_scored)

# Split into the two clinical cohorts
df_ghosts = df_scored[df_scored['ACMG_Classification'] == 'Bucket 3: PM2 (Ghost Variant)']
df_borderline = df_scored[df_scored['ACMG_Classification'] == 'Bucket 2: Borderline Carrier']

# Helper function to safely get counts
def get_class_count(df, class_name):
    return len(df[df['am_class'] == class_name])

# ==========================================
# STEP 5: ISOLATE THE SOLVED CASES
# ==========================================
df_solved = df_final[
    (df_final['ACMG_Classification'] == 'Bucket 2: Borderline Carrier') & 
    (df_final['am_class'] == 'likely_pathogenic')
]

df_solved = df_solved.sort_values(by='am_pathogenicity', ascending=False)
df_solved.to_csv(FINAL_SOLVED_PATH, sep='\t', index=False)

end_time = time.time()
minutes = (end_time - start_time) / 60

# ==========================================
# FINAL CONSOLE OUTPUT
# ==========================================
print("\n✅ PROJECT PIPELINE COMPLETE!")
print(f"Time taken: {minutes:.2f} minutes.")
print("=" * 60)
print(" 🧬 OVERALL BIOINFORMATIC YIELD")
print("=" * 60)
print(f"Total Targets Interrogated: {len(df_targets):,}")
print(f"Variants without an AI Score (Likely Frameshifts/Non-missense): {non_missense_count:,}")
print(f"Total Missense Variants Scored by AI: {len(df_scored):,}\n")

print("-" * 60)
print(" 👻 BUCKET 3: THE GHOSTS (PM2 - Absent in Middle East)")
print("-" * 60)
print(f"Total Scored: {len(df_ghosts):,}")
print(f"  - Likely Pathogenic: {get_class_count(df_ghosts, 'likely_pathogenic'):,}")
print(f"  - Ambiguous:         {get_class_count(df_ghosts, 'ambiguous'):,}")
print(f"  - Likely Benign:     {get_class_count(df_ghosts, 'likely_benign'):,}")
if not df_ghosts.empty:
    print(f"  - Mean Pathogenicity Score: {df_ghosts['am_pathogenicity'].mean():.4f}")

print("\n" + "-" * 60)
print(" 🎯 BUCKET 2: BORDERLINE CARRIERS (Circulating in Middle East)")
print("-" * 60)
print(f"Total Scored: {len(df_borderline):,}")
print(f"  - Likely Pathogenic: {get_class_count(df_borderline, 'likely_pathogenic'):,}  <-- CANDIDATE SET")
print(f"  - Ambiguous:         {get_class_count(df_borderline, 'ambiguous'):,}")
print(f"  - Likely Benign:     {get_class_count(df_borderline, 'likely_benign'):,}")
if not df_borderline.empty:
    print(f"  - Mean Pathogenicity Score: {df_borderline['am_pathogenicity'].mean():.4f}")

print("=" * 60)
print(f"💾 Final candidate dataset written to: {FINAL_SOLVED_PATH}")