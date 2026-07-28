import pandas as pd
import time

# Define your exact local paths
INPUT_PATH = r"...[YOUR_PATH]...\alphamissense-variant-pipeline\Assets\00_source_variant_summary.txt"
OUTPUT_PATH = r"...[YOUR_PATH]...\alphamissense-variant-pipeline\Assets\01_clinvar_vus_extracted.txt"

# 32GB RAM can comfortably handle 500,000 rows in memory at once
CHUNK_SIZE = 500000

print("🚀 Starting the ClinVar extraction pipeline...")
start_time = time.time()

total_processed = 0
total_kept = 0
first_chunk = True

# Read the TSV file in chunks
# sep='\t' tells pandas it's a Tab-Separated file
# low_memory=False stops pandas from guessing column datatypes and slowing down
for chunk in pd.read_csv(INPUT_PATH, sep='\t', chunksize=CHUNK_SIZE, low_memory=False):
    
    # 1. Apply the Double Filter (GRCh38 AND Pure VUS)
    filtered_chunk = chunk[
        (chunk['Assembly'] == 'GRCh38') & 
        (chunk['ClinicalSignificance'] == 'Uncertain significance')
    ]
    
    # 2. Save to Output
    # If it's the first chunk, write the header. Otherwise, append to the bottom.
    if first_chunk:
        filtered_chunk.to_csv(OUTPUT_PATH, sep='\t', index=False, mode='w')
        first_chunk = False
    else:
        filtered_chunk.to_csv(OUTPUT_PATH, sep='\t', index=False, mode='a', header=False)
        
    # Update counters and print live progress
    total_processed += len(chunk)
    total_kept += len(filtered_chunk)
    print(f"Processed: {total_processed:,} rows | Extracted: {total_kept:,} pure VUSs...")

end_time = time.time()
minutes = (end_time - start_time) / 60

print("\n✅ Extraction Complete!")
print(f"Time taken: {minutes:.2f} minutes")
print(f"Final Count: Reduced the raw file into {total_kept:,} high-value targets.")