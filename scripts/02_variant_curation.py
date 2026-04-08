import pandas as pd
import re

# ── Load raw data ──────────────────────────────────────────
df = pd.read_csv('raw_data/clinvar_gch1.tsv', sep='\t')
print(f"Starting with {len(df)} total variants")

# ── Step 1: Keep only single nucleotide variants ───────────
df_snv = df[df['Variant type'] == 'single nucleotide variant']
print(f"After keeping SNVs only: {len(df_snv)}")

# ── Step 2: Keep only rows with a protein change (p.) ──────
df_mis = df_snv[df_snv['Name'].str.contains(r'p\.', na=False)]
print(f"After keeping variants with protein change: {len(df_mis)}")

# ── Step 3: Assign labels ──────────────────────────────────
path_terms = ['Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic']
ben_terms  = ['Benign', 'Likely benign', 'Benign/Likely benign']

df_path = df_mis[df_mis['Germline classification'].isin(path_terms)].copy()
df_ben  = df_mis[df_mis['Germline classification'].isin(ben_terms)].copy()
print("\nSample of benign variant names:")
print(df_ben['Name'].head(20).to_string())
df_vus  = df_mis[df_mis['Germline classification'] == 'Uncertain significance'].copy()

df_path['label'] = 1
df_ben['label']  = 0

print(f"\nPathogenic variants: {len(df_path)}")
print(f"Benign variants: {len(df_ben)}")
print(f"VUS (set aside): {len(df_vus)}")

# ── Step 4: Combine pathogenic and benign ──────────────────
df_clean = pd.concat([df_path, df_ben], ignore_index=True)
print(f"\nCombined dataset: {len(df_clean)} variants")

# ── Step 5: Parse HGVS strings into R198Q format ──────────
AA3 = {
    'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C',
    'Gln':'Q','Glu':'E','Gly':'G','His':'H','Ile':'I',
    'Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P',
    'Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V'
}

def parse_hgvs(hgvs_str):
    """Convert p.Arg198Gln style to R198Q, extract position and residues."""
    if not isinstance(hgvs_str, str):
        return None, None, None, None
    m = re.search(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvs_str)
    if not m:
        return None, None, None, None
    wt3, pos, mut3 = m.group(1), m.group(2), m.group(3)
    wt  = AA3.get(wt3)
    mut = AA3.get(mut3)
    if not (wt and mut):
        return None, None, None, None
    return f'{wt}{pos}{mut}', int(pos), wt, mut

# Apply parser to every row
parsed = df_clean['Name'].apply(lambda x: pd.Series(parse_hgvs(x),
         index=['mutation', 'pos', 'wt_aa', 'mut_aa']))
df_clean = pd.concat([df_clean, parsed], axis=1)

# ── Step 6: Drop rows that couldn't be parsed ──────────────
before = len(df_clean)
df_clean = df_clean.dropna(subset=['mutation'])
print(f"Dropped {before - len(df_clean)} unparseable rows")
print(f"Final clean dataset: {len(df_clean)} variants")

# ── Step 7: Check for duplicates ──────────────────────────
dupes = df_clean[df_clean.duplicated(subset=['mutation'], keep=False)]
if len(dupes) > 0:
    print(f"\nWarning: {len(dupes)} duplicate mutations found:")
    print(dupes[['mutation', 'Germline classification', 'label']])
    # Keep the first occurrence
    df_clean = df_clean.drop_duplicates(subset=['mutation'], keep='first')
    print(f"After removing duplicates: {len(df_clean)} variants")
else:
    print("\nNo duplicates found")

# ── Step 8: Keep only useful columns ──────────────────────
cols = ['mutation', 'pos', 'wt_aa', 'mut_aa', 'label',
        'Germline classification', 'Name', 'Condition(s)']
df_final = df_clean[cols].copy()

# ── Step 9: Print a sample to verify it looks right ────────
print("\nSample of cleaned data:")
print(df_final.head(10).to_string())

print("\nLabel counts in final dataset:")
print(df_final['label'].value_counts())

# ── Step 10: Save outputs ──────────────────────────────────
df_final.to_csv('processed_data/variants_clean.csv', index=False)
df_vus.to_csv('processed_data/variants_vus.csv', index=False)
print("\nSaved:")
df_final.to_csv('processed_data/variants_clean.csv', index=False)
df_vus.to_csv('processed_data/variants_vus.csv', index=False)

#finding if there are genuine missense benign variants hiding in the data
print("\nProtein change column for benign variants:")
print(df_ben[['Name', 'Protein change', 'Germline classification']].head(20).to_string())