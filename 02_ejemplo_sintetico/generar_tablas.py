from pathlib import Path
import numpy as np
import pandas as pd

# Datos artificiales para comprobar formatos y ejecutar modelos, no para simular GEMAS.
carpeta = Path(__file__).resolve().parent / "datos"
carpeta.mkdir(exist_ok=True)
rng = np.random.default_rng(20261024)
n_saliva, n_controles, n_asv = 206, 12, 320
ids = [f"SYN_{i:04d}" for i in range(1, n_saliva + 1)]
ids += [f"SYN_CTRL_{i:03d}" for i in range(1, n_controles + 1)]
asv = [f"ASV_SYN_{i:04d}" for i in range(1, n_asv + 1)]
n = len(ids)
fechas = pd.date_range("2018-03-01", "2020-06-30", freq="9D")
fecha_muestra = rng.choice(fechas, size=n).astype("datetime64[ns]")
fecha_muestra[n_saliva - 1:] = np.datetime64("NaT")
exac = rng.binomial(1, 0.35, n)
sexo = rng.permutation(np.resize(["Female", "Male"], n))
meta = pd.DataFrame({
    "sample.type": ["Saliva"] * n_saliva + ["Mock"] * 6 + ["NTC"] * 6,
    "Age": rng.integers(18, 81, n), "Sex.letter": sexo,
    "Pool": rng.choice(["1", "2", "3", "6"], n),
    "Antibiotic.2mo": rng.binomial(1, 0.25, n),
    "Cavities": rng.binomial(1, 0.30, n),
    "Chewing.gum.30min": rng.binomial(1, 0.12, n),
    "Exacerbations.6M": exac,
    "Exacerbations.12M": np.maximum(exac, rng.binomial(1, 0.40, n)),
    "Date": fecha_muestra, "Island": rng.choice(["Tenerife", "La Palma"], n),
    "synthetic": True, "Sex.number": (sexo == "Male").astype(int)
}, index=pd.Index(ids, name="sample_id"))

# Conteos multinomiales con heterogeneidad y ausencias; sin efectos clínicos insertados.
pesos = np.exp(rng.normal(-1, 1.3, n_asv))
presencia = rng.uniform(0.18, 0.95, n_asv)
counts = np.zeros((n_asv, n), dtype=np.int64)
for j in range(n):
    activo = rng.binomial(1, presencia).astype(bool)
    activo[:8] = True
    probabilidades = rng.gamma(0.6 + pesos)
    probabilidades[~activo] = 0
    profundidad = rng.integers(25000, 80001) if j < n_saliva else rng.integers(1000, 10001)
    counts[:, j] = rng.multinomial(profundidad, probabilidades / probabilidades.sum())
k = np.arange(n_asv)
tax = pd.DataFrame({"Kingdom": "Bacteria",
    "Phylum": [f"Phylum_SYN_{i // 80 + 1:02d}" for i in k],
    "Class": [f"Class_SYN_{i // 40 + 1:02d}" for i in k],
    "Order": [f"Order_SYN_{i // 20 + 1:02d}" for i in k],
    "Family": [f"Family_SYN_{i // 10 + 1:02d}" for i in k],
    "Genus": [f"Genus_SYN_{i // 5 + 1:03d}" for i in k],
    "Species": [f"Species_SYN_{i + 1:04d}" for i in k]},
    index=pd.Index(asv, name="taxon_id"))

# Árbol artificial bifurcante con longitudes positivas; no procede de secuencias reales.
nodos = list(asv)
while len(nodos) > 1:
    indices = sorted(rng.choice(len(nodos), 2, replace=False), reverse=True)
    a, b = [nodos.pop(int(i)) for i in indices]
    nodos.append(f"({a}:{rng.uniform(0.01, 1):.8f},{b}:{rng.uniform(0.01, 1):.8f})")
(carpeta / "arbol.nwk").write_text(nodos[0] + ";\n", encoding="utf-8")

# Ambiente también sintético: cadenas de episodios y PM10 generado de forma independiente.
calendario = pd.date_range("2016-01-01", "2020-12-31", freq="D")
ambientes = []
for isla in ["Tenerife", "La Palma"]:
    evento = np.zeros(len(calendario), dtype=int)
    for i in range(len(evento)):
        evento[i] = rng.binomial(1, 0.60 if i and evento[i-1] else 0.09)
    naf = evento * rng.gamma(3, 8, len(evento))
    pm10 = 8 + rng.gamma(3, 4, len(evento)) + naf
    inicio = (evento == 1) & np.r_[True, evento[:-1] == 0]
    numero = np.cumsum(inicio)
    ambientes.append(pd.DataFrame({"fecha": calendario, "isla": isla,
        "miteco_africano": evento, "aporte_naf_usable": naf, "PM10": pm10,
        "episodio_miteco_id": [f"{isla}_SYN_E{e}" if yes else "" for e, yes in zip(numero, evento)],
        "synthetic": True}))
meta.to_csv(carpeta / "metadatos.csv", date_format="%Y-%m-%d")
tax.to_csv(carpeta / "taxonomia.csv")
pd.DataFrame(counts, index=pd.Index(asv, name="taxon_id"), columns=ids).to_csv(carpeta / "conteos.csv")
pd.concat(ambientes).to_csv(carpeta / "ambiente_diario.csv", index=False, date_format="%Y-%m-%d")
assert len(meta.loc[meta["sample.type"] == "Saliva"]) == 206
assert meta["Date"].notna().sum() == 205
assert counts[:, :206].sum(axis=0).min() >= 10000
print("Datos artificiales escritos en", carpeta)
