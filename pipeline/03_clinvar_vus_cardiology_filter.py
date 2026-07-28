import pandas as pd
import time

# Define your input (the 2.3M row file) and the new output file
INPUT_PATH = r"...[YOUR_PATH]...\alphamissense-variant-pipeline\Assets\01_clinvar_vus_extracted.txt"
OUTPUT_PATH = r"...[YOUR_PATH]...\Research\AI Genetics Project\alphamissense-variant-pipeline\Assets\03_clinvar_vus_cardiology.txt"

# The Master Cardiovascular "Trawl Net" (93 Genes)
CARDIO_GENES = {
    # 1. The Electric Heart
    "SCN5A", "SCN10A", "SCN1B", "SCN2B", "SCN3B", "SCN4B", "GPD1L", "SNTA1", "RANGRF",
    "KCNQ1", "KCNH2", "KCNJ2", "KCNJ5", "KCNJ8", "KCNE1", "KCNE2", "KCNE3", "KCND3", "ABCC9",
    "CACNA1C", "CACNB2", "CACNA2D1", "RYR2", "CASQ2", "CALM1", "CALM2", "CALM3", "TRDN", "TECRL", "AKAP9", "ANK2",
    "HCN4",
    
    # 2. The Structural Heart
    "MYH7", "MYBPC3", "TNNT2", "TNNI3", "TNNC1", "TPM1", "ACTC1", "MYL2", "MYL3", "MYPN", "NEXN", "ALPK3",
    "PKP2", "DSP", "DSG2", "DSC2", "JUP", "TMEM43",
    "TTN", "DES", "FLNC", "BAG3", "VCL", "CSRP3", "TCAP", "LDB3",
    "LMNA", "RBM20", "EMD",
    "PRKAG2", "LAMP2", "GLA", "TTR", "GAA",
    "PLN",
    
    # 3. The Vascular Heart
    "FBN1", "FBN2", "COL3A1", "MFAP5", "LOX",
    "TGFBR1", "TGFBR2", "TGFB2", "TGFB3", "SMAD3", "SMAD4", "SMAD6",
    "ACTA2", "MYH11", "MYLK", "PRKG1",
    
    # 4. The Metabolic Heart
    "LDLR", "APOB", "PCSK9", "LDLRAP1", "APOE", "LPL", "APOC2", "GPIHBP1", "ABCG5", "ABCG8"
}

CHUNK_SIZE = 500000

print("🎣 Dropping the Cardiovascular Trawl Net...")
start_time = time.time()

total_scanned = 0
total_cardio_vus = 0
first_chunk = True

# Read through the Extracted VUS file
for chunk in pd.read_csv(INPUT_PATH, sep='\t', chunksize=CHUNK_SIZE, low_memory=False):
    
    # Filter: Keep row ONLY if the GeneSymbol is exactly in our CARDIO_GENES set
    filtered_chunk = chunk[chunk['GeneSymbol'].isin(CARDIO_GENES)]
    
    # Save to Output
    if first_chunk:
        filtered_chunk.to_csv(OUTPUT_PATH, sep='\t', index=False, mode='w')
        first_chunk = False
    else:
        filtered_chunk.to_csv(OUTPUT_PATH, sep='\t', index=False, mode='a', header=False)
        
    total_scanned += len(chunk)
    total_cardio_vus += len(filtered_chunk)
    print(f"Scanned: {total_scanned:,} rows | Captured Cardiology VUSs: {total_cardio_vus:,}...")

end_time = time.time()

print("\n✅ Filtration Complete!")
print(f"Time taken: {(end_time - start_time):.2f} seconds")
print(f"Final Count: You have isolated {total_cardio_vus:,} pure cardiovascular targets.")