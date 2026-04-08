import pandas as pd

# Load ClinVar data
df = pd.read_csv('../raw_data/clinvar_gch1.tsv', sep='\t')
print(f"Total rows: {df.shape[0]}")
print(f"Total columns: {df.shape[1]}")
print("\nColumn names:")
for col in df.columns.tolist():
    print(f"  {col}")

# Load FASTA sequence
with open('../raw_data/gch1_uniprot.fasta') as f:
    lines = f.readlines()

sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
print(f"\nGCH1 sequence length: {len(sequence)} amino acids")
print(f"First 30 aa: {sequence[:30]}")

# Clinical significance breakdown
print("\nClinical significance breakdown:")
print(df['Germline classification'].value_counts())