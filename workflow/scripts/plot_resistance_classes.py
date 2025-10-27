# ../scripts/plot_resistance_classes.py
import os
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # évite les erreurs de backend sur serveur
import matplotlib.pyplot as plt

input_file = snakemake.input["abricate"]
output_file = snakemake.output["plot"]

os.makedirs(os.path.dirname(output_file), exist_ok=True)

def save_no_data(msg):
    """Créer un graphique 'No data' en cas d'erreur ou fichier vide"""
    plt.figure(figsize=(6, 4))
    plt.text(0.5, 0.5, msg, ha="center", va="center", fontsize=12)
    plt.axis("off")
    plt.savefig(output_file, dpi=200, bbox_inches="tight")
    plt.close()

# 1️⃣ Vérifie si le fichier existe et non vide
if (not os.path.exists(input_file)) or os.path.getsize(input_file) == 0:
    save_no_data("No data: empty or missing file")
    sys.exit(0)

print(f"Reading file: {input_file}")

# 2️⃣ Lecture du fichier avec gestion des erreurs
try:
    df = pd.read_csv(input_file, sep="\t")
    df.columns = df.columns.str.strip().str.upper()
except Exception as e:
    print(f"Error reading file: {e}")
    save_no_data(f"No data (read error): {e}")
    sys.exit(0)

# 3️⃣ Vérifie la présence des colonnes nécessaires
if not {"GENE", "RESISTANCE"}.issubset(df.columns):
    print("⚠️ Missing required columns after normalization.")
    print("Found columns:", df.columns.tolist())
    save_no_data("No data: missing 'GENE' or 'RESISTANCE' column")
    sys.exit(0)

# 4️⃣ Gère le cas de dataframe vide
if df.empty:
    save_no_data("No data: empty dataframe")
    sys.exit(0)

# 5️⃣ Prépare les données pour le comptage
df["RESISTANCE"] = df["RESISTANCE"].astype(str).str.split(";")
df = df.explode("RESISTANCE")
df["RESISTANCE"] = df["RESISTANCE"].str.strip()
counts = df["RESISTANCE"].value_counts()

# 6️⃣ Si aucun gène, crée graphique vide
if counts.empty:
    save_no_data("No data: no resistance classes found")
    sys.exit(0)

# 7️⃣ Création du graphique principal
plt.figure(figsize=(10, 6))
counts.sort_values().plot(kind="barh", color="tomato")
plt.xlabel("Nombre de gènes")
plt.title("Gènes de résistance par classe d’antibiotiques")
plt.tight_layout()
plt.savefig(output_file, dpi=300)
plt.close()

print(f"Plot saved to: {output_file}")
