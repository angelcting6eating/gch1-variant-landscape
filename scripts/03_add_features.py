import pandas as pd

# ── Load cleaned variants ──────────────────────────────────
df = pd.read_csv('processed_data/variants_clean.csv')
print(f"Loaded {len(df)} variants")
print(df[['mutation', 'wt_aa', 'mut_aa', 'label']].head())

# ── BLOSUM62 matrix ────────────────────────────────────────
from Bio.Align import substitution_matrices
blosum62 = substitution_matrices.load("BLOSUM62")

def get_blosum62(wt, mut):
    try:
        return blosum62[wt][mut]
    except:
        return None

# ── Grantham distance matrix ───────────────────────────────
GRANTHAM = {
    ('A','R'):112,('A','N'):111,('A','D'):126,('A','C'):195,
    ('A','Q'):91, ('A','E'):107,('A','G'):60, ('A','H'):86,
    ('A','I'):94, ('A','L'):96, ('A','K'):106,('A','M'):84,
    ('A','F'):113,('A','P'):27, ('A','S'):99, ('A','T'):58,
    ('A','W'):148,('A','Y'):112,('A','V'):64,
    ('R','N'):86, ('R','D'):96, ('R','C'):180,('R','Q'):43,
    ('R','E'):54, ('R','G'):125,('R','H'):29, ('R','I'):97,
    ('R','L'):102,('R','K'):26, ('R','M'):91, ('R','F'):97,
    ('R','P'):103,('R','S'):110,('R','T'):71, ('R','W'):101,
    ('R','Y'):77, ('R','V'):96,
    ('N','D'):23, ('N','C'):139,('N','Q'):46, ('N','E'):42,
    ('N','G'):80, ('N','H'):68, ('N','I'):149,('N','L'):153,
    ('N','K'):94, ('N','M'):142,('N','F'):158,('N','P'):91,
    ('N','S'):46, ('N','T'):65, ('N','W'):174,('N','Y'):143,
    ('N','V'):133,
    ('D','C'):154,('D','Q'):61, ('D','E'):45, ('D','G'):94,
    ('D','H'):81, ('D','I'):168,('D','L'):172,('D','K'):101,
    ('D','M'):160,('D','F'):177,('D','P'):108,('D','S'):65,
    ('D','T'):85, ('D','W'):181,('D','Y'):160,('D','V'):152,
    ('C','Q'):154,('C','E'):170,('C','G'):159,('C','H'):174,
    ('C','I'):198,('C','L'):198,('C','K'):202,('C','M'):196,
    ('C','F'):205,('C','P'):169,('C','S'):112,('C','T'):149,
    ('C','W'):215,('C','Y'):194,('C','V'):192,
    ('Q','E'):29, ('Q','G'):87, ('Q','H'):24, ('Q','I'):109,
    ('Q','L'):113,('Q','K'):53, ('Q','M'):101,('Q','F'):116,
    ('Q','P'):76, ('Q','S'):68, ('Q','T'):42, ('Q','W'):130,
    ('Q','Y'):99, ('Q','V'):96,
    ('E','G'):98, ('E','H'):40, ('E','I'):134,('E','L'):138,
    ('E','K'):56, ('E','M'):126,('E','F'):140,('E','P'):93,
    ('E','S'):80, ('E','T'):65, ('E','W'):152,('E','Y'):122,
    ('E','V'):121,
    ('G','H'):98, ('G','I'):135,('G','L'):138,('G','K'):127,
    ('G','M'):127,('G','F'):153,('G','P'):42, ('G','S'):56,
    ('G','T'):59, ('G','W'):184,('G','Y'):147,('G','V'):109,
    ('H','I'):94, ('H','L'):99, ('H','K'):32, ('H','M'):87,
    ('H','F'):100,('H','P'):77, ('H','S'):89, ('H','T'):47,
    ('H','W'):115,('H','Y'):83, ('H','V'):84,
    ('I','L'):5,  ('I','K'):102,('I','M'):10, ('I','F'):21,
    ('I','P'):95, ('I','S'):142,('I','T'):89, ('I','W'):61,
    ('I','Y'):33, ('I','V'):29,
    ('L','K'):107,('L','M'):15, ('L','F'):22, ('L','P'):98,
    ('L','S'):145,('L','T'):92, ('L','W'):61, ('L','Y'):36,
    ('L','V'):32,
    ('K','M'):95, ('K','F'):102,('K','P'):103,('K','S'):121,
    ('K','T'):78, ('K','W'):110,('K','Y'):85, ('K','V'):97,
    ('M','F'):28, ('M','P'):87, ('M','S'):135,('M','T'):81,
    ('M','W'):67, ('M','Y'):36, ('M','V'):21,
    ('F','P'):114,('F','S'):155,('F','T'):103,('F','W'):40,
    ('F','Y'):22, ('F','V'):50,
    ('P','S'):74, ('P','T'):38, ('P','W'):147,('P','Y'):110,
    ('P','V'):68,
    ('S','T'):58, ('S','W'):177,('S','Y'):144,('S','V'):124,
    ('T','W'):128,('T','Y'):102,('T','V'):69,
    ('W','Y'):37, ('W','V'):88,
    ('Y','V'):55
}

def get_grantham(wt, mut):
    if wt == mut:
        return 0
    key = (wt, mut) if (wt, mut) in GRANTHAM else (mut, wt)
    return GRANTHAM.get(key, None)

# ── Biochemical property dictionaries ─────────────────────
CHARGE = {
    'D':-1,'E':-1,
    'K':1,'R':1,'H':0.1,
    'A':0,'C':0,'F':0,'G':0,'I':0,'L':0,
    'M':0,'N':0,'P':0,'Q':0,'S':0,'T':0,
    'V':0,'W':0,'Y':0
}

HYDRO = {
    'I':4.5,'V':4.2,'L':3.8,'F':2.8,'C':2.5,
    'M':1.9,'A':1.8,'G':-0.4,'T':-0.7,'S':-0.8,
    'W':-0.9,'Y':-1.3,'P':-1.6,'H':-3.2,'E':-3.5,
    'Q':-3.5,'D':-3.5,'N':-3.5,'K':-3.9,'R':-4.5
}

def delta_charge(wt, mut):
    return CHARGE.get(mut, 0) - CHARGE.get(wt, 0)

def delta_hydro(wt, mut):
    return HYDRO.get(mut, 0) - HYDRO.get(wt, 0)

# ── Apply all features ─────────────────────────────────────
df['blosum62']     = df.apply(lambda r: get_blosum62(r['wt_aa'], r['mut_aa']), axis=1)
df['grantham']     = df.apply(lambda r: get_grantham(r['wt_aa'], r['mut_aa']), axis=1)
df['delta_charge'] = df.apply(lambda r: delta_charge(r['wt_aa'], r['mut_aa']), axis=1)
df['delta_hydro']  = df.apply(lambda r: delta_hydro(r['wt_aa'], r['mut_aa']), axis=1)

# ── Also add R198Q separately as case study ────────────────
r198q = {
    'mutation': 'R198Q',
    'pos': 198,
    'wt_aa': 'R',
    'mut_aa': 'Q',
    'label': -1,  # -1 means case study, not in training set
    'Germline classification': 'Uncertain significance',
    'Name': 'NM_000161.3(GCH1):c.593G>A (p.Arg198Gln)',
    'Condition(s)': 'GTP cyclohydrolase I deficiency',
    'blosum62':     get_blosum62('R', 'Q'),
    'grantham':     get_grantham('R', 'Q'),
    'delta_charge': delta_charge('R', 'Q'),
    'delta_hydro':  delta_hydro('R', 'Q')
}

df_r198q = pd.DataFrame([r198q])
df_all = pd.concat([df, df_r198q], ignore_index=True)

# ── Print results ──────────────────────────────────────────
print("\nFeature summary for pathogenic variants:")
print(df[['mutation','blosum62','grantham','delta_charge','delta_hydro']].to_string())

print("\nR198Q feature values:")
print(df_r198q[['mutation','blosum62','grantham','delta_charge','delta_hydro']].to_string())

print("\nAny missing values?")
print(df[['blosum62','grantham','delta_charge','delta_hydro']].isnull().sum())

# ── Save ───────────────────────────────────────────────────
df_all.to_csv('processed_data/features_table.csv', index=False)
print("\nSaved to processed_data/features_table.csv")