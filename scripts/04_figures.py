import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ── Load data ──────────────────────────────────────────────
df = pd.read_csv('processed_data/features_table.csv')

# Separate pathogenic and R198Q
path = df[df['label'] == 1].copy()
r198q = df[df['mutation'] == 'R198Q'].iloc[0]

print(f"Pathogenic variants: {len(path)}")
print(f"R198Q: {r198q['mutation']}")

# ── Define features and labels ─────────────────────────────
features = ['blosum62', 'grantham', 'delta_charge', 'delta_hydro']
titles = [
    'BLOSUM62 Score\n(lower = less tolerated)',
    'Grantham Distance\n(higher = bigger physicochemical change)',
    'Δ Charge\n(mutant − wild-type)',
    'Δ Hydrophobicity\n(Kyte-Doolittle, mutant − wild-type)'
]

# ── Create 2x2 figure ──────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('Physicochemical Landscape of Pathogenic GCH1 Missense Variants\nwith R198Q in context',
             fontsize=14, fontweight='bold', y=1.01)

for ax, feat, title in zip(axes.flat, features, titles):

    values = path[feat].dropna()
    r_val  = r198q[feat]

    # ── Histogram of pathogenic distribution ──
    ax.hist(values, bins=15, color='#f06470', alpha=0.7,
            edgecolor='white', label='Pathogenic variants')

    # ── Mark R198Q as vertical line ───────────
    ax.axvline(r_val, color='#1a1a2e', linewidth=2.5,
               linestyle='--', label=f'R198Q ({r_val})')

    # ── Percentile annotation ──────────────────
    pct = (values < r_val).mean() * 100
    ax.text(0.97, 0.95, f'R198Q = {r_val:.1f}\n{pct:.0f}th percentile',
            transform=ax.transAxes, ha='right', va='top',
            fontsize=9, color='#1a1a2e',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='#1a1a2e', alpha=0.8))

    ax.set_title(title, fontsize=10, fontweight='bold')
    ax.set_xlabel(feat, fontsize=9)
    ax.set_ylabel('Number of variants', fontsize=9)
    ax.legend(fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('figures/r198q_landscape.png', dpi=150, bbox_inches='tight')
print("Saved figures/r198q_landscape.png")

# ── Figure 2: R198Q feature profile bar chart ─────────────
fig2, ax2 = plt.subplots(figsize=(8, 5))

# Normalise each feature to 0-1 scale for comparison
def normalise(series):
    mn, mx = series.min(), series.max()
    if mx == mn:
        return series * 0
    return (series - mn) / (mx - mn)

norm_path = pd.DataFrame({
    feat: normalise(path[feat]) for feat in features
})

r198q_norm = {
    feat: (r198q[feat] - path[feat].min()) /
          (path[feat].max() - path[feat].min())
    if path[feat].max() != path[feat].min() else 0
    for feat in features
}

x = np.arange(len(features))
width = 0.35

bars1 = ax2.bar(x - width/2,
                [norm_path[f].mean() for f in features],
                width, label='Pathogenic mean',
                color='#f06470', alpha=0.8, edgecolor='white')

bars2 = ax2.bar(x + width/2,
                [r198q_norm[f] for f in features],
                width, label='R198Q',
                color='#1a1a2e', alpha=0.8, edgecolor='white')

ax2.set_xticks(x)
ax2.set_xticklabels(['BLOSUM62', 'Grantham', 'Δ Charge', 'Δ Hydro'],
                    fontsize=10)
ax2.set_ylabel('Normalised value (0-1 scale)', fontsize=10)
ax2.set_title('R198Q vs Mean Pathogenic Variant Profile\n(features normalised to 0-1)',
              fontsize=11, fontweight='bold')
ax2.legend(fontsize=10)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('figures/r198q_profile.png', dpi=150, bbox_inches='tight')
print("Saved figures/r198q_profile.png")

# ── Print R198Q percentiles ────────────────────────────────
print("\nR198Q percentiles within pathogenic distribution:")
for feat in features:
    val = r198q[feat]
    pct = (path[feat] < val).mean() * 100
    print(f"  {feat}: {val:.2f} → {pct:.0f}th percentile")

print("\nPathogenic variant summary statistics:")
print(path[features].describe().round(2))