# =====================================================================
# kalibracio_epizod.py
# ---------------------------------------------------------------------
# MIT CSINÁL?
#   Ellenőrzi (validálja), hogy az ERA5-ből definiált "epizódnapok"
#   (szélcsendes + hideg napok, legalább 3 egymás után) tényleg SZENNYEZETTEBB
#   napok-e a MÉRT PM2,5 alapján (OLM-állomások). Ugyanezzel hasonlítja össze
#   a különböző küszöbértékeket (szigorú / fő / laza definíció).
#
# MIÉRT KELL?
#   A szakirodalom nem ad számértéket az "országos átlagos szélre", ezért a
#   küszöböket (2 m/s, 5 °C, 3 nap) nem vehetjük készen. Ez a szkript a
#   független mérésekkel mutatja meg, hogy a választás értelmes-e.
#   Műholdas/modelladatot NEM használ, csak az ERA5-táblát és az állomásokat.
#
# A TELJES LOGIKA EGY MONDATBAN:
#   1) kiválasztjuk az epizódnapokat; 2) minden állomás napi értékéből levonjuk
#   az állomás "szokásos" szintjét (alapvonal) -> anomália; 3) az anomáliákat
#   átlagoljuk országos napi mutatóvá; 4) összevetjük az epizód- és a nem
#   epizódnapok átlagát; 5) a bizonytalanságot az EPIZÓDOK újramintavételével
#   (bootstrap) becsüljük, mert a szomszédos napok nem függetlenek.
#
# BEMENET (a szkripttel egy mappában):
#   ERA5_napi_szel_homerseklet.csv  - a GEE-ből exportált napi tábla (date, wind, temp)
#   pm25_napi_hosszu.csv            - mért napi PM2,5, állomás x nap (hosszú formátum)
#   allomasok_tisztitott.csv        - állomáslista (tipus_csoport oszlop a szűréshez)
# KIMENET:
#   a képernyőre táblázatok + kalibracio_eredmenyek.csv
# Futtatás: python kalibracio_epizod.py   (pandas, numpy, scipy kell)
# =====================================================================
import numpy as np
import pandas as pd
from scipy import stats

# ---------------- 1. BEÁLLÍTÁSOK (ezeket szabad átírni) -----------------
ERA_CSV      = 'ERA5_napi_szel_homerseklet.csv'
PM_CSV       = 'pm25_napi_hosszu.csv'
STATIONS_CSV = 'allomasok_tisztitott.csv'
OUT_CSV      = 'kalibracio_eredmenyek.csv'

# Mely állomástípusokat zárjuk ki? Pl. ['közlekedési'] ha csak háttérállomások kellenek.
# Az eddigi eredmények (6. táblázat) MINDEN állomással készültek, ezért üres.
EXCLUDE_GROUPS = []

# Melyik hónapokban értékelünk? A 6. táblázat november-februárra készült.
# None = minden elérhető hónap (okt. 15. - ápr. 15.). Az epizódok KIVÁLASZTÁSA
# mindig a teljes ERA5-táblán történik (hogy a hónaphatáron átnyúló epizódok se vesszenek el),
# csak az ÉRTÉKELÉS korlátozódik erre.
ANALYSIS_MONTHS = [11, 12, 1, 2]

# A vizsgált definíciók: név -> (szélküszöb m/s, minimális egymást követő napok, hőmérsékleti küszöb K)
# 278,15 K = 5 °C. A szél az országos átlagos napi szélsebesség (ERA5-Land).
DEFINITIONS = {
    'Szigorú (<1,2 m/s)': (1.2, 3, 278.15),
    'Fő (<2,0 m/s)':      (2.0, 3, 278.15),
    'Laza (<2,5 m/s)':    (2.5, 3, 278.15),
}
MAIN_DEFINITION = 'Fő (<2,0 m/s)'   # az alapvonal-összehasonlítás ezen készül

MIN_STATIONS = 10     # egy nap csak akkor számít, ha legalább ennyi állomás mért aznap
MIN_NON_DAYS = 5      # alapvonalhoz kell legalább ennyi nem epizódnap egy állomás-tél-hónapban
TOP_QUANTILE = 0.90   # a "legszennyezettebb napok" = az anomália felső 10%-a
N_BOOT       = 5000   # bootstrap-ismétlések száma
SEED         = 1      # a véletlen magja: így a futás ismételhető (a CI ~±0,1-et mozoghat más magnál)
RUN_GRID     = False  # True: feltáró rácskeresés sok küszöbkombinációra (lassabb)


# ---------------- 2. ADATBEOLVASÁS ---------------------------------------
def load_era(path):
    """Az ERA5-Land napi országos átlagok (szél m/s, hőmérséklet K)."""
    e = pd.read_csv(path, usecols=['date', 'wind', 'temp'])
    e['date'] = pd.to_datetime(e['date'])
    # dropna: ha egy napra hiányzik az érték, nem tehetjük epizódnappá.
    # sort_values: a "három egymást követő nap" logika sorbarendezett táblát igényel.
    return e.dropna().sort_values('date').reset_index(drop=True)


def load_pm(pm_csv, stations_csv):
    """Mért PM2,5 állomás x nap; csak a hozzárendelt állomások és az értékelt hónapok."""
    pm = pd.read_csv(pm_csv, parse_dates=['datum'])
    pm = pm.dropna(subset=['allomas_id', 'pm25'])      # nem hozzárendelt oszlopok / hiányzó értékek ki
    pm['allomas_id'] = pm['allomas_id'].astype(int)
    st = pd.read_csv(stations_csv, usecols=['id', 'tipus_csoport'], encoding='utf-8-sig')
    pm = pm.merge(st, left_on='allomas_id', right_on='id', how='left')
    if EXCLUDE_GROUPS:
        pm = pm[~pm['tipus_csoport'].isin(EXCLUDE_GROUPS)]
    pm['honap'] = pm['datum'].dt.month
    # "tel" = a tél KEZDŐÉVE (2021 = a 2021/22-es tél). A július előtti hónapok az előző évi télhez tartoznak.
    pm['tel'] = np.where(pm['honap'] >= 7, pm['datum'].dt.year, pm['datum'].dt.year - 1)
    if ANALYSIS_MONTHS:
        pm = pm[pm['honap'].isin(ANALYSIS_MONTHS)]
    return pm[['datum', 'allomas_id', 'pm25', 'honap', 'tel']].copy()


# ---------------- 3. EPIZÓDNAPOK KIVÁLASZTÁSA ----------------------------
def pick_days(era, wmax, window, tmax):
    """A GEE-szkript pickDays függvényének pontos mása.
    Epizódnap = olyan nap, amely része egy legalább `window` napos sorozatnak,
    amelynek MINDEN napján szél < wmax ÉS hőmérséklet < tmax."""
    # ok[i] igaz, ha az i-edik nap "szélcsendes és hideg"
    ok = ((era['wind'] < wmax) & (era['temp'] < tmax)).to_numpy()
    dates = era['date'].to_numpy()
    selected = set()
    for i in range(len(era) - window + 1):          # minden lehetséges `window` hosszú ablak
        if ok[i:i + window].all():                  # az ablak minden napja teljesíti-e a feltételt?
            # Naptári egymásutániság: az első és utolsó nap különbsége pontosan window-1 nap legyen.
            # Ez megakadályozza, hogy két tél (vagy egy hiányzó nap) között "összeragadjon" egy ablak.
            span = int((dates[i + window - 1] - dates[i]) / np.timedelta64(1, 'D'))
            if span == window - 1:
                selected.update(dates[i:i + window])   # az ablak MINDEN napja epizódnap (így a hosszabb epizódok is teljesen bekerülnek)
    return pd.DatetimeIndex(sorted(selected))


# ---------------- 4. ANOMÁLIA = napi érték - az állomás "szokásos" szintje ----
def add_anomaly(pm, episode_days, mode):
    """Az anomália az állomás napi PM2,5-értéke mínusz az alapvonal.
    mode:
      'tel_honap'      JAVASOLT: az állomás ugyanazon telének ugyanazon hónapjában mért
                       NEM epizódnapi átlag; ha abban kevesebb mint MIN_NON_DAYS nap van,
                       az állomás adott hónapjának telekre összevont nem epizódos átlaga (tartalék).
      'tel_honap_mind' RÉGI: ugyanaz, de az epizódnapok is benne vannak az alapvonalban
                       (az alapvonalat felhúzza -> a hatást alulbecsli).
      'honap_osszevont' a telek összevonva, minden nap (nem szűri ki az évek közötti
                       különbséget -> a hatást felülbecsli, mert az epizódok szennyezettebb telekre esnek).
    Visszaad: (DataFrame anomáliával, a tartalék-alapvonal aránya az epizódos állomásnapokon)."""
    d = pm.copy()
    d['epi'] = d['datum'].isin(episode_days)        # igaz, ha az adott nap epizódnap
    fb = 0.0
    if mode == 'tel_honap_mind':
        base = d.groupby(['allomas_id', 'tel', 'honap'])['pm25'].transform('mean')
    elif mode == 'honap_osszevont':
        base = d.groupby(['allomas_id', 'honap'])['pm25'].transform('mean')
    elif mode == 'tel_honap':
        non = d[~d['epi']]                           # csak a NEM epizódnapok adják az alapvonalat
        cell = non.groupby(['allomas_id', 'tel', 'honap'])['pm25'].agg(b_cell='mean', n_cell='count')
        pool = non.groupby(['allomas_id', 'honap'])['pm25'].mean().rename('b_pool')
        d = d.join(cell, on=['allomas_id', 'tel', 'honap']).join(pool, on=['allomas_id', 'honap'])
        use_cell = d['n_cell'] >= MIN_NON_DAYS       # elég nem epizódnap van a cellában? (NaN -> False -> tartalék)
        base = np.where(use_cell, d['b_cell'], d['b_pool'])
        fb = float((~use_cell)[d['epi']].mean())     # az epizódos állomásnapok mekkora részén kellett a tartalék
    else:
        raise ValueError(mode)
    d['anom'] = d['pm25'] - base
    return d, fb


def daily_index(d):
    """Állomásonkénti anomáliákból országos napi mutató: az aznap mérő állomások átlaga."""
    g = d.groupby('datum').agg(a=('anom', 'mean'), n=('anom', 'count'), epi=('epi', 'first')).reset_index()
    g = g[g['n'] >= MIN_STATIONS].copy()             # túl kevés állomásos napot eldobunk (zajos lenne)
    g['tel'] = np.where(g['datum'].dt.month >= 7, g['datum'].dt.year, g['datum'].dt.year - 1)
    return g


# ---------------- 5. KIÉRTÉKELÉS -----------------------------------------
def evaluate(daily, rng, daily_neutral):
    """Az epizód- és a nem epizódnapok anomáliájának összehasonlítása.
    daily_neutral: ugyanaz a napi mutató, de az epizód-definíciótól FÜGGETLEN alapvonallal
    (állomás x tél x hónap átlag, minden nap) - ez kell a 'lefedés' számításához."""
    non_mean = daily.loc[~daily['epi'], 'a'].mean()  # a nem epizódnapok átlagos anomáliája (összehasonlítási alap)
    e = daily[daily['epi']].sort_values('datum').copy()
    # Egymást követő epizódnapok = EGY epizód. A diff != 1 nap új epizódot indít; a cumsum számozza őket.
    e['grp'] = (e['datum'].diff().dt.days != 1).cumsum()
    per = e.groupby('grp')['a'].mean().to_numpy()    # epizódonkénti átlagos anomália (az epizód az egység, nem a nap)
    excess = per.mean() - non_mean                   # a "PM2,5-többlet": az epizódok átlaga mínusz a nem epizódnapok átlaga
    # Bootstrap: az epizódokat visszatevéssel újra húzzuk, minden húzásra újraszámoljuk a többletet.
    # A 2,5. és 97,5. percentilis a 95%-os konfidenciaintervallum. (Az epizódok szintjén dolgozunk,
    # mert egy epizód napjai egymástól nem függetlenek: a napokat külön húzni túl szűk intervallumot adna.)
    boot = [rng.choice(per, per.size).mean() - non_mean for _ in range(N_BOOT)]
    lo, hi = np.percentile(boot, [2.5, 97.5])
    npos = int((per > non_mean).sum())               # hány epizód szennyezettebb a nem epizódnapok átlagánál
    # Előjelpróba: ha az epizódoknak nem lenne hatása, mindegyik 50%-os eséllyel lenne pozitív.
    p_sign = stats.binomtest(npos, per.size, 0.5, alternative='greater').pvalue
    # Lefedés: a legszennyezettebb napok (az anomália felső 10%-a) mekkora részét fogja meg a definíció?
    # FONTOS: a "legszennyezettebb napokat" az epizód-definíciótól FÜGGETLEN anomáliából jelöljük ki
    # (daily_neutral). Ha ugyanazt az alapvonalat használnánk, mint a többlethez (amelyből az epizódnapok
    # ki vannak zárva), a mérték önmagát igazolná (körkörös lenne), és felfújt értékeket adna.
    thr = daily_neutral['a'].quantile(TOP_QUANTILE)
    top_days = set(daily_neutral.loc[daily_neutral['a'] >= thr, 'datum'])
    epi_days = set(daily.loc[daily['epi'], 'datum'])
    coverage = len(top_days & epi_days) / len(top_days)
    # Télenkénti megoszlás: mennyire koncentrálódnak az epizódnapok egyetlen télre?
    by_winter = e.groupby('tel').size()
    # Leave-one-winter-out: egy tél kihagyásával is megmarad-e a nap-súlyozott különbség?
    lowo = []
    for y in sorted(daily['tel'].unique()):
        s = daily[daily['tel'] != y]
        if s['epi'].sum() > 0:
            lowo.append(s.loc[s['epi'], 'a'].mean() - s.loc[~s['epi'], 'a'].mean())
    return dict(epizod=per.size, nap=len(e), tel=by_winter.size,
                max_telresz=by_winter.max() / len(e), tobblet=excess, ci_lo=lo, ci_hi=hi,
                pozitiv=f'{npos}/{per.size}', p_elojel=p_sign, lefedes=coverage,
                min_LOWO=min(lowo) if lowo else np.nan)


def run(pm, ep, rng, mode='tel_honap'):
    """Egy teljes futás egy epizódlistára: anomália (a kért alapvonallal) + definíciófüggetlen anomália
    a lefedéshez + kiértékelés."""
    d, fb = add_anomaly(pm, ep, mode)
    dn, _ = add_anomaly(pm, ep, 'tel_honap_mind')
    res = evaluate(daily_index(d), rng, daily_index(dn))
    res['tartalek_arany'] = fb
    return res


def continuous_relationship(pm, era):
    """Folytonos kapcsolat: a napi PM-anomália és a szél/hőmérséklet Spearman-korrelációja.
    Az anomália itt az epizód-definíciótól FÜGGETLEN (állomás x tél x hónap átlag, minden nap)."""
    d = pm.copy()
    d['anom'] = d['pm25'] - d.groupby(['allomas_id', 'tel', 'honap'])['pm25'].transform('mean')
    d['epi'] = False
    daily = daily_index(d).merge(era, left_on='datum', right_on='date')
    # 3 napos, középre igazított mozgóátlag-szél: a "tartós" szélcsendet méri
    w3 = era.set_index('date')['wind'].rolling(3, center=True).mean().rename('w3')
    daily = daily.merge(w3, left_on='datum', right_index=True)
    return {n: stats.spearmanr(daily[c], daily['a'], nan_policy='omit')[0]
            for n, c in [('szél (napi)', 'wind'), ('szél (3 napos átlag)', 'w3'), ('hőmérséklet', 'temp')]}


# ---------------- 6. FUTTATÁS --------------------------------------------
if __name__ == '__main__':
    pd.set_option('display.width', 200)
    era = load_era(ERA_CSV)
    pm = load_pm(PM_CSV, STATIONS_CSV)
    rng = np.random.default_rng(SEED)
    print(f'ERA5: {len(era)} nap | PM2,5: {pm.allomas_id.nunique()} állomás, {pm.datum.nunique()} nap, {len(pm)} állomásnap')

    # --- (a) Folytonos kapcsolat ---
    print('\n(a) A napi PM-anomália Spearman-korrelációja az időjárással:')
    for k, v in continuous_relationship(pm, era).items():
        print(f'    {k:22s} {v:+.2f}')

    # --- (b) A definíciók összehasonlítása (a javasolt alapvonallal) ---
    rows = []
    for name, (wmax, window, tmax) in DEFINITIONS.items():
        ep = pick_days(era, wmax, window, tmax)          # 1) epizódnapok az ERA5-ből
        res = run(pm, ep, rng)                            # 2-5) anomália, napi mutató, kiértékelés
        res.update(definicio=name, epizodnap_ERA5=len(ep))
        rows.append(res)
    R = pd.DataFrame(rows).set_index('definicio')
    show = R[['epizod', 'nap', 'tel', 'max_telresz', 'tobblet', 'ci_lo', 'ci_hi', 'pozitiv', 'p_elojel', 'lefedes', 'min_LOWO', 'tartalek_arany']]
    print('\n(b) Definíciók összehasonlítása (alapvonal: állomás x tél x hónap, nem epizódnapok):')
    print(show.round(3).to_string())
    R.to_csv(OUT_CSV, encoding='utf-8-sig')

    # --- (c) Az alapvonal megválasztásának hatása a fő definícióra ---
    wmax, window, tmax = DEFINITIONS[MAIN_DEFINITION]
    ep = pick_days(era, wmax, window, tmax)
    print(f'\n(c) Alapvonal-választás hatása ({MAIN_DEFINITION}), PM2,5-többlet µg/m3:')
    for mode, leiras in [('honap_osszevont', 'telek összevonva, minden nap (felülbecsül)'),
                         ('tel_honap_mind', 'tél x hónap, epizódnapokkal együtt (alulbecsül)'),
                         ('tel_honap', 'tél x hónap, csak nem epizódnapok (javasolt)')]:
        r = run(pm, ep, np.random.default_rng(SEED), mode)
        print(f'    {mode:16s} {r["tobblet"]:+.1f} ({r["ci_lo"]:+.1f}...{r["ci_hi"]:+.1f})   {leiras}')

    # --- (d) A tartalék-küszöb érzékenysége ---
    print(f'\n(d) Az alapvonal tartalék-küszöbének (MIN_NON_DAYS) hatása ({MAIN_DEFINITION}):')
    keep = MIN_NON_DAYS
    for n in [1, 3, 5, 7, 10, 14]:
        MIN_NON_DAYS = n
        r = run(pm, ep, np.random.default_rng(SEED))
        print(f'    küszöb {n:>2} nap: többlet {r["tobblet"]:+.2f}, tartalék az epizódos állomásnapok {r["tartalek_arany"]:.1%}-ában')
    MIN_NON_DAYS = keep

    # --- (e) Opcionális feltáró rácskeresés ---
    if RUN_GRID:
        print('\n(e) Feltáró rácskeresés (minden sor egy küszöbkombináció):')
        out = []
        for w in [1.0, 1.2, 1.5, 1.8, 2.0, 2.5]:
            for win in [2, 3, 4, 5]:
                for tk, t in [('nincs', 999.0), ('8 C', 281.15), ('5 C', 278.15), ('2 C', 275.15)]:
                    e_ = pick_days(era, w, win, t)
                    if len(e_) < 3:
                        continue
                    r = run(pm, e_, np.random.default_rng(SEED))
                    out.append(dict(szel=w, ablak=win, hom=tk, **{k: r[k] for k in ['epizod', 'nap', 'tel', 'tobblet', 'ci_lo', 'lefedes']}))
        print(pd.DataFrame(out).round(2).to_string(index=False))
        print('FIGYELEM: a rácskeresés feltáró; sok kombináció közül a "legjobb" kiválasztása torzít.')
