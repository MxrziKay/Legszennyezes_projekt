# =====================================================================
# Fatüzelési kibocsátás (településszint) -> belterületi poligonok ->
#   (1) állomások körüli pufferek, (2) 5 km-es rács
# Futtatás: python emisszio_csatolas.py     (geopandas, pandas, numpy, openpyxl kell)
# =====================================================================
import re, os, unicodedata, difflib
import numpy as np, pandas as pd, geopandas as gpd
from shapely.geometry import box

# ---------------- 1. BEÁLLÍTÁSOK (ezeket írd át) ---------------------
XLSX        = 'fatuzeles_PM25_BC_forgatokonyvek.xlsx'
BELTERULET  = 'belterulet.gpkg'          # belterületi poligonok (gpkg/shp)
BELT_LAYER  = None                       # None = első réteg
BELT_NEV    = 'nev'                      # a település nevét tartalmazó oszlop
STATIONS    = 'allomasok_tisztitott.csv' # lon, lat, id oszlopokkal
GRID_FILE   = None                       # a GEE-ből exportált 5 km-es rács (shp/gpkg); None = saját EOV-rács
EMISSZIO    = 'PM2.5 (kg) – HU készülékmegoszlás szerint'   # az elemzett oszlop
BUFFERS_KM  = [1, 3, 5, 10, 20]
CRS         = 'EPSG:23700'               # EOV, méter
OUT         = 'kimenet'
os.makedirs(OUT, exist_ok=True)

# ---------------- 2. NÉVNORMALIZÁLÁS ---------------------------------
ROM = {'I':1,'V':5,'X':10,'L':50}
def roman(s):
    t = 0
    for i, c in enumerate(s):
        v = ROM[c]; t += -v if i+1 < len(s) and ROM[s[i+1]] > v else v
    return t

def strip_acc(s):
    return ''.join(c for c in unicodedata.normalize('NFKD', s) if not unicodedata.combining(c))

def key(s):
    s = re.sub(r'\s+', ' ', str(s)).strip()
    m = re.match(r'^budapest\s+(\d{1,2})\.?\s*(ker\.?|kerület)?$', s, re.I)       # "Budapest 01. ker."
    if m: return 'budapest_%d' % int(m.group(1))
    m = re.match(r'^budapest\s+([IVXL]+)\.?\s*(ker\.?|kerület)?$', s, re.I)        # "Budapest XVIII. kerület"
    if m: return 'budapest_%d' % roman(m.group(1).upper())
    m = re.match(r'^([IVXL]+)\.?\s*kerület$', s, re.I)                              # "XVIII. kerület"
    if m: return 'budapest_%d' % roman(m.group(1).upper())
    s = re.sub(r'\s*\(.*?\)', '', s)                  # zárójeles kiegészítés törlése
    s = re.sub(r'[\u2010-\u2015]', '-', s)            # kötőjel-változatok
    return s.lower()

# ---------------- 3. TÁBLA BEOLVASÁSA --------------------------------
t = pd.read_excel(XLSX, sheet_name='Adattábla', header=1)
t = t.rename(columns={t.columns[0]: 'nev'})
t['key'] = t.nev.map(key)
assert not t.key.duplicated().any(), 'duplikált kulcs a táblában'
t['pop_proxy'] = t['Fával'].fillna(0) + t['Hálózati (vezetékes) gázzal és fával'].fillna(0)

# ---------------- 4. POLIGONOK + ILLESZTÉS ----------------------------
b = gpd.read_file(BELTERULET, layer=BELT_LAYER).to_crs(CRS)
b['key'] = b[BELT_NEV].map(key)
b = b[~b.geometry.is_empty & b.geometry.notna()].copy()

# Ha a poligonréteg egyetlen "Budapest"-et tartalmaz, a kerületeket összevonjuk
if 'budapest' in set(b.key) and not any(k.startswith('budapest_') for k in b.key):
    bp = t.key.str.startswith('budapest_')
    s = t[bp].select_dtypes('number').sum()
    row = {c: s[c] for c in s.index}; row.update(nev='Budapest', key='budapest')
    t = pd.concat([t[~bp], pd.DataFrame([row])], ignore_index=True)
    print('Budapest kerületei összevonva egyetlen sorrá.')

bk, tk = set(b.key), set(t.key)
# második kör: ékezet nélküli egyezés
acc = {strip_acc(k): k for k in bk}
fix = {k: acc[strip_acc(k)] for k in tk - bk if strip_acc(k) in acc}
t['key'] = t.key.replace(fix); tk = set(t.key)

tabla_nelkul = sorted(tk - bk)      # van adat, nincs poligon
poli_nelkul  = sorted(bk - tk)      # van poligon, nincs adat
print(f'Illesztve: {len(tk & bk)} / {len(tk)} település | táblában poligon nélkül: {len(tabla_nelkul)} | poligon adat nélkül: {len(poli_nelkul)}')
lost = t[t.key.isin(tabla_nelkul)][EMISSZIO].sum() / t[EMISSZIO].sum() * 100
print(f'Nem illesztett települések kibocsátása: {lost:.2f}% (cél: < 1%)')
rows = []
for k in tabla_nelkul:
    javaslat = difflib.get_close_matches(k, poli_nelkul, n=3, cutoff=0.8)
    rows.append({'tabla_kulcs': k, 'kibocsatas_kg': t.loc[t.key == k, EMISSZIO].iloc[0], 'javaslat': ' | '.join(javaslat)})
pd.DataFrame(rows).to_csv(f'{OUT}/nem_illesztett_tabla.csv', index=False, encoding='utf-8-sig')
pd.DataFrame({'poligon_kulcs': poli_nelkul}).to_csv(f'{OUT}/nem_illesztett_poligon.csv', index=False, encoding='utf-8-sig')

# ---------------- 5. KIBOCSÁTÁS ELOSZTÁSA A POLIGONRÉSZEK KÖZÖTT -------
b['a_km2'] = b.geometry.area / 1e6
tot = b.groupby('key').a_km2.sum().rename('belter_km2')
b = b.join(tot, on='key').merge(t[['key', EMISSZIO, 'pop_proxy']], on='key', how='inner')
b['dens'] = b[EMISSZIO] / b.belter_km2                    # kg / km² belterület, településen belül egyenletes
print(f'Poligonrészek: {len(b)}; települések: {b.key.nunique()}; belterület összesen: {b.a_km2.sum():.0f} km²')
chk = b.assign(em=b.dens * b.a_km2).groupby('key').em.sum().sum()
print(f'Ellenőrzés: elosztott kibocsátás = {chk:,.0f} kg  (táblában az illesztetteké: {t[t.key.isin(set(b.key))][EMISSZIO].sum():,.0f} kg)')
b = b[['key', 'dens', 'geometry']]

def aggregate(windows):
    """windows: GeoDataFrame (poligonok). Visszaad: kibocsátás kg, belterület km² ablakonként (index szerint)."""
    w = windows[['geometry']].copy(); w['wid'] = w.index
    ov = gpd.overlay(w, b, how='intersection', keep_geom_type=True)
    ov['a'] = ov.geometry.area / 1e6
    ov['em'] = ov.dens * ov.a
    return ov.groupby('wid').agg(em_kg=('em', 'sum'), belter_km2=('a', 'sum')).reindex(windows.index).fillna(0)

# ---------------- 6. ÁLLOMÁSOK: PUFFEREK -------------------------------
try:
    st = pd.read_csv(STATIONS, encoding='utf-8-sig')
except UnicodeDecodeError:           # Excelből mentett (Windows-1250) CSV
    st = pd.read_csv(STATIONS, encoding='cp1250')
g = gpd.GeoDataFrame(st, geometry=gpd.points_from_xy(st.lon, st.lat), crs='EPSG:4326').to_crs(CRS)
res = st.copy()
for r in BUFFERS_KM:
    win = gpd.GeoDataFrame(geometry=g.geometry.buffer(r * 1000), index=g.index, crs=CRS)
    a = aggregate(win)
    res[f'em_kg_{r}km'] = a.em_kg.values
    res[f'em_dens_{r}km'] = a.em_kg.values / (np.pi * r * r)          # kg / km² (a teljes kör területére)
    res[f'belter_{r}km'] = a.belter_km2.values
res.to_csv(f'{OUT}/allomas_emisszio_pufferek.csv', index=False, encoding='utf-8-sig')
print('Mentve: allomas_emisszio_pufferek.csv')

# ---------------- 7. 5 km-ES RÁCS ---------------------------------------
if GRID_FILE:
    grid = gpd.read_file(GRID_FILE).to_crs(CRS)
else:
    x0, y0, x1, y1 = b.total_bounds
    xs = np.arange(np.floor(x0/5000)*5000, x1, 5000); ys = np.arange(np.floor(y0/5000)*5000, y1, 5000)
    grid = gpd.GeoDataFrame(geometry=[box(x, y, x+5000, y+5000) for x in xs for y in ys], crs=CRS)
a = aggregate(grid)
grid['em_kg'] = a.em_kg.values
grid['belter_km2'] = a.belter_km2.values
grid['em_dens'] = grid.em_kg / (grid.geometry.area / 1e6)
grid['belter_arany'] = grid.belter_km2 / (grid.geometry.area / 1e6)
grid = grid[grid.em_kg > 0] if not GRID_FILE else grid
grid.to_file(f'{OUT}/racs_emisszio.gpkg', driver='GPKG')
grid.drop(columns='geometry').to_csv(f'{OUT}/racs_emisszio.csv', index=False, encoding='utf-8-sig')
print(f'Rács: {len(grid)} cella; elosztott kibocsátás a rácson: {grid.em_kg.sum():,.0f} kg')
