# =====================================================================
# A GEE-ből exportált tél x hónap cellatáblák (cellak_tel_*.csv) összegzése
# rácscellánkénti ANOMÁLIÁVÁ:
#   Δ(x) = Σ_k w_k * [X_epizód,k(x) − B_k(x)] / Σ_k w_k
#   k = (tél, hónap) cella, w_k = a cella epizódnapjainak száma (n_epi)
#   B_k = az ugyanazon tél ugyanazon hónapjának nem epizódnapi átlaga;
#         ha a cellában kevés a nem epizódnap / érvényes megfigyelés, a hónap
#         telekre összevont nem epizódos átlaga (súlyozva) a tartalék-alapvonal.
# Futtatás: python kombinal_cellak.py
# =====================================================================
import glob
import numpy as np
import pandas as pd

# ---------------- BEÁLLÍTÁSOK ----------------------------------------
INPUT_GLOB   = 'cellak_tel_*.csv'
OUTPUT       = 'racs_anomalia.csv'
MIN_NON_DAYS = 5       # CAMS: legalább ennyi nem epizódnap a cellában (különben tartalék-alapvonal)
MIN_OBS_EPI  = 1.0     # S5P: pixelenkénti átlagos érvényes megfigyelések minimuma az epizódnapokon
MIN_OBS_NON  = 3.0     # S5P: ugyanez a nem epizódnapokon (az alapvonalhoz)
MIN_COV      = 0.5     # az epizódnapok ekkora hányadát kell lefednie az érvényes celláknak

# Mértékegység-átváltás a termékekre. A CAMS PM2,5-öt a GEE kg/m3-ben adja (pl. 1,6e-8 kg/m3 = 16 ug/m3),
# ezért 1e9-cel szorozzuk, hogy ug/m3 legyen. Az UVAI dimenziótlan, az NO2 mol/m2 marad.
UNIT_FACTOR = {'pm': 1e9, 'uvai': 1.0, 'no2': 1.0}

# termék: (epizód-oszlop, nem epizód-oszlop, epizód-count, nem epizód-count)
PRODUCTS = {
    'pm':   ('pm_epi',   'pm_non',   None,         None),
    'uvai': ('uvai_epi', 'uvai_non', 'n_uvai_epi', 'n_uvai_non'),
    'no2':  ('no2_epi',  'no2_non',  'n_no2_epi',  'n_no2_non'),
}

# ---------------- BEOLVASÁS --------------------------------------------
files = sorted(glob.glob(INPUT_GLOB))
assert files, 'Nincs a mintára illő fájl: ' + INPUT_GLOB
df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
print(f'{len(files)} fájl, {len(df)} sor, {df.cell_id.nunique()} rácscella, '
      f'{df.groupby(["tel", "honap"]).ngroups} tél x hónap cella')

# ---------------- ÖSSZEGZÉS TERMÉKENKÉNT -----------------------------------
out = pd.DataFrame({'cell_id': sorted(df.cell_id.unique())}).set_index('cell_id')
for name, (e, n, ce, cn) in PRODUCTS.items():
    d = df[['cell_id', 'tel', 'honap', 'n_epi', 'n_non', e, n] + [c for c in (ce, cn) if c]].copy()
    d[e] = d[e] * UNIT_FACTOR[name]          # mértékegység-átváltás (CAMS: kg/m3 -> ug/m3)
    d[n] = d[n] * UNIT_FACTOR[name]

    # érvényesség
    epi_ok = d[e].notna()
    non_ok = d[n].notna()
    if ce:
        epi_ok &= d[ce] >= MIN_OBS_EPI
        non_ok &= d[cn] >= MIN_OBS_NON
        d['w_non'] = d[cn].where(non_ok, 0.0)          # alapvonal-súly: érvényes megfigyelések
        local_ok = non_ok                               # a megfigyelésszám dönt
    else:
        d['w_non'] = d['n_non'].where(non_ok, 0.0)      # CAMS: napok száma
        local_ok = non_ok & (d['n_non'] >= MIN_NON_DAYS)

    # tartalék-alapvonal: hónap telekre összevont nem epizódos átlag (cellánként)
    tmp = d.loc[non_ok].assign(wx=lambda x: x['w_non'] * x[n])
    pool = tmp.groupby(['cell_id', 'honap']).agg(wx=('wx', 'sum'), ww=('w_non', 'sum'))
    pool['b_pool'] = pool.wx / pool.ww
    d = d.join(pool['b_pool'], on=['cell_id', 'honap'])

    d['baseline'] = np.where(local_ok, d[n], d['b_pool'])
    d['fallback'] = ~local_ok & d['b_pool'].notna()
    d['valid'] = epi_ok & d['baseline'].notna()
    d['delta'] = d[e] - d['baseline']

    # súlyozott összegzés rácscellánként
    g = d.groupby('cell_id')
    tot_w = g['n_epi'].sum()
    ok = d[d.valid]
    num = (ok['n_epi'] * ok['delta']).groupby(ok.cell_id).sum()
    den = ok.groupby('cell_id')['n_epi'].sum()
    out['d_' + name] = (num / den)
    out['cov_' + name] = (den / tot_w).reindex(out.index).fillna(0.0)
    out['fb_' + name] = d.groupby('cell_id')['fallback'].mean()   # tartalék-alapvonalat használó cellák aránya
    n_fb = int(d['fallback'].sum()); n_all = int(d['valid'].sum())
    print(f'[{name}] érvényes tél-hónap-rácscella: {n_all}, ebből tartalék-alapvonallal: {n_fb} ({n_fb / max(n_all, 1):.1%})')
    if n_fb / max(n_all, 1) > 0.05:
        print(f'  FIGYELEM [{name}]: a cellák több mint 5%-a tartalék-alapvonalat használ. A telekre összevont alapvonal visszahozza az évek közötti összetétel-torzítást; vizsgáld meg a küszöböket (MIN_NON_DAYS / MIN_OBS_NON).')
    if ce:
        print(f'[{name}] átlagos érvényes megfigyelés/pixel – epizód: median {d[ce].median():.1f}, nem epizód: median {d[cn].median():.1f}')

    # túl kevés lefedettségű rácscellák jelölése
    out.loc[out['cov_' + name] < MIN_COV, 'd_' + name] = np.nan

out.reset_index().to_csv(OUTPUT, index=False, encoding='utf-8-sig')
print('Mentve:', OUTPUT)
print(out.filter(like='d_').describe().round(3).to_string())
