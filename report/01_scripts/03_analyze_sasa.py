
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

# ===============================================
# 2MYO - Relative Solvent Accessibility Analysis
# ===============================================

REPORT = Path(__file__).resolve().parent.parent

INPUT = REPORT / "01_utilities" / "06_sasa_2myo_per_residue.csv"
OUTPUT = REPORT / "01_utilities" / "07_sasa_2myo_rsa.csv"
FIGURE = REPORT / "01_utilities" / "08_rsa_distribution.png"

FIGURE.parent.mkdir(parents=True, exist_ok=True)

# Theoretical maximum ASA values (Tien et al., 2013)
MAX_ASA = {
    "ALA": 129, "ARG": 274, "ASN": 195,
    "ASP": 193, "CYS": 167, "GLN": 225,
    "GLU": 223, "GLY": 104, "HIS": 224,
    "ILE": 197, "LEU": 201, "LYS": 236,
    "MET": 224, "PHE": 240, "PRO": 159,
    "SER": 155, "THR": 172, "TRP": 285,
    "TYR": 263, "VAL": 174,
}

# Ankyrin-like structural regions
REGIONS = [
    ("N-terminal", 1, 27, "#8094d8"),
    ("Repeat 1", 28, 60, "#df8585"),
    ("Repeat 2", 61, 93, "#88c98b"),
    ("C-terminal", 94, 118, "#bc91ce"),
]

def get_region(resid):
    for name, start, end, _ in REGIONS:
        if start <= resid <= end:
            return name
    raise ValueError(f"Unexpected residue number: {resid}")

# Load original SASA data
df = pd.read_csv(INPUT)
df = df.sort_values("Residue").reset_index(drop=True)

assert len(df) == 118, "Expected 118 residues."
assert df["Residue"].tolist() == list(range(1, 119)), \
    "Residue numbering is incomplete or duplicated."

# Calculate RSA
df["MaxASA_A2"] = df["ResName"].map(MAX_ASA)

if df["MaxASA_A2"].isna().any():
    raise ValueError("Unknown amino acid in CSV.")

df["RSA"] = df["SASA_A2"] / df["MaxASA_A2"]

df["Region"] = df["Residue"].apply(get_region)

def classify(rsa):
    if rsa < 0.05:
        return "Buried"
    if rsa < 0.25:
        return "Partially exposed"
    return "Exposed"

df["Exposure"] = df["RSA"].apply(classify)

# Save annotated CSV
df.to_csv(OUTPUT, index=False, float_format="%.4f")

# Print results
print(f"Residues analyzed: {len(df)}")
print(f"Total SASA: {df['SASA_A2'].sum():.2f} A^2")

print("\nExposure categories:")
print(df["Exposure"].value_counts())

print("\nMean RSA by region:")
print(df.groupby("Region", sort=False)["RSA"].mean().round(3))

if (df["RSA"] > 1).any():
    print("\nWarning: Some RSA values exceed 1.")
    print("Check differences in SASA calculation conventions.")

# ===============================================
# Plot RSA along the protein sequence
# ===============================================

fig, ax = plt.subplots(figsize=(12, 4.5))

# Shade ankyrin-like regions
for name, start, end, color in REGIONS:
    ax.axvspan(
        start - 0.5,
        end + 0.5,
        color=color,
        alpha=0.20,
        label=name,
    )

# Plot RSA
ax.plot(
    df["Residue"],
    df["RSA"],
    color="#294c73",
    linewidth=1.3,
    marker="o",
    markersize=2.5,
)

# Reference threshold for solvent exposure
ax.axhline(
    0.25,
    color="gray",
    linestyle="--",
    linewidth=1,
    label="RSA = 0.25",
)

ax.set_xlabel("Residue position")
ax.set_ylabel("Relative Solvent Accessibility (RSA)")
ax.set_title("Myotrophin (2MYO) - Solvent Accessibility")

ax.set_xlim(1, 118)
ax.set_ylim(0, max(1.0, df["RSA"].max() * 1.05))

ax.legend(fontsize=8, ncol=5, loc="upper center",
          bbox_to_anchor=(0.5, -0.17))
ax.grid(axis="y", alpha=0.2)

fig.tight_layout()
fig.savefig(FIGURE, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"\nCSV saved: {OUTPUT}")
print(f"Figure saved: {FIGURE}")