# Az adatfeldolgozás fájljainak leírása

Ez a dokumentum a szkriptek által létrehozott és a kutatásban használt táblák tartalmát írja le: mi egy sor, mit jelentenek az oszlopok, milyen mértékegységben szerepelnek az értékek.

## 1. Fájltérkép

| Fájl | Ki készíti | Egy sor = | Mire való |
| --- | --- | --- | --- |
| `ERA5_napi_szel_homerseklet.csv` | GEE-szkript | egy nap | az országos átlagos szél és hőmérséklet |
| `epizodnapok_lista.csv` | GEE-szkript | egy epizódnap | a kiválasztott epizódnapok |
| `cellak_tel_ÉÉÉÉ.csv` (télenként egy) | GEE-szkript | rácscella × hónap | epizód- és nem epizódnapi kompozitok |
| `HU_5km_racs_geometria` (shp) | GEE-szkript | rácscella | a cellák alakja (`cell_id` a kulcs) |
| `racs_anomalia.csv` | `kombinal_cellak.py` | rácscella | a tél × hónap cellákból összegzett anomália |
| `kalibracio_eredmenyek.csv` | `kalibracio_epizod.py` | epizód-definíció | az epizód-definíciók összehasonlítása mért PM2,5-tel |
| `allomasok_tisztitott.csv` | előkészített bemenet | állomás | az állomások adatai és a mérési lefedettség |
| `pm25_napi_hosszu.csv` | előkészített bemenet | állomás × nap | a mért napi PM2,5 |
| `allomas_emisszio_pufferek.csv` | `emisszio_csatolas.py` | állomás | a kibocsátás az állomások körüli pufferekben |
| `racs_emisszio.csv` / `.gpkg` | `emisszio_csatolas.py` | rácscella | a kibocsátás az 5 × 5 km-es rácson |
| `nem_illesztett_tabla.csv` | `emisszio_csatolas.py` | település | kibocsátási táblában szereplő, de poligon nélküli települések |
| `nem_illesztett_poligon.csv` | `emisszio_csatolas.py` | poligon-név | poligonréteg-beli, de adat nélküli nevek |

Az `emisszio_csatolas.py` műveleteit QGIS-ben végezném el (de feltöltöttem a repo-ba)
A `kombinal_cellak.py`-t még nem futtattam le
### A folyamat egy pillantásra

```
ERA5-tábla ──► epizódnapok ──► cellak_tel_*.csv ──► racs_anomalia.csv
     │                                                     ▲
     └──► kalibracio_epizod.py ◄── pm25_napi_hosszu.csv    │ (cell_id)
                 (a küszöbök ellenőrzése)         HU_5km_racs_geometria
KSH-adat ──► kibocsátás ──► emisszio_csatolas.py ──► racs_emisszio, allomas_emisszio_pufferek
```

A kalibráció (`kalibracio_epizod.py`) a GEE-szkripttől külön fut: csak az ERA5-táblát és a mért PM2,5-et használja, és az eredménye (a küszöbértékek) kerül vissza a GEE-szkript elejére.

## 2. GEE-kimenetek

### 2.1. `ERA5_napi_szel_homerseklet.csv`

Egy sor egy nap. Az exportált példában 911 sor van: az öt tél (2021/22–2025/26) október 15. és április 14. közötti napjai, hiány nélkül.

| Oszlop | Jelentés | Egység |
| --- | --- | --- |
| `date` | a nap dátuma | ÉÉÉÉ-HH-NN |
| `wind` | Magyarország területére átlagolt, 10 m-es szélsebesség a napi átlagos komponensekből (s = √(u² + v²)); a valódi napi átlagos sebességet kissé alábecsüli | m/s |
| `temp` | az országos átlagos 2 m-es léghőmérséklet | K (278,15 K = 5 °C) |

### 2.2. `epizodnapok_lista.csv`

Egy sor egy kiválasztott epizódnap (a 2,0 m/s, 3 nap, 5 °C alatti szabály szerint).

| Oszlop | Jelentés |
| --- | --- |
| `system:index` | a GEE sorszáma (0-tól indul, ezért az utolsó érték eggyel kisebb a sorok számánál) |
| `date` | az epizódnap dátuma |
| `.geo` | technikai oszlop (üres geometria), törölhető |

Az epizódnapok száma a fájl **adatsorainak** száma (fejléc nélkül), nem az utolsó `system:index`.

### 2.3. `cellak_tel_ÉÉÉÉ.csv`

Egy sor egy 5 × 5 km-es rácscella egy hónapra, egy télre. A 2025-ös fájlban 15 788 sor van: 3947 cella × 4 hónap (2025. november – 2026. február).

| Oszlop | Jelentés | Egység |
| --- | --- | --- |
| `cell_id` | a cella középpontja EPSG:3035-ben, méterben (`x_y`), pl. `4932500_2547500` | m |
| `tel` | a tél kezdőéve (2025 = 2025/26) | év |
| `honap` | a naptári hónap (a 2025-ös fájlban 1 = 2026. január) | 1–12 |
| `n_epi` | az epizódnapok száma az adott hónapban | nap |
| `n_non` | a nem epizódnapok száma az adott hónapban | nap |
| `pm_epi` | a CAMS PM2,5 átlaga az epizódnapokon | **kg/m³** (×10⁹ = µg/m³) |
| `pm_non` | a CAMS PM2,5 átlaga a nem epizódnapokon (alapvonal) | kg/m³ |
| `uvai_epi`, `uvai_non` | az UV aeroszol index átlaga epizód- / nem epizódnapokon | dimenziótlan |
| `n_uvai_epi`, `n_uvai_non` | a cellában az átlagba bekerült érvényes megfigyelések átlagos száma pixelenként | db |
| `no2_epi`, `no2_non` | a troposzférikus NO2-oszlop átlaga | mol/m² |
| `n_no2_epi`, `n_no2_non` | az NO2 érvényes megfigyeléseinek száma | db |
| `system:index`, `.geo` | technikai oszlopok, törölhetők | |

Értelmezés: a cella adott hónapjának anomáliája a `*_epi − *_non` különbség. A `n_epi` és az `n_non` az összes cellában azonos az adott hónapra, mert országos meteorológiai kritériumból származnak.

Példa (január, `cell_id = 4932500_2547500`): 17 epizód- és 14 nem epizódnap; CAMS 17,9 és 11,8 µg/m³ (különbség +6,1); UVAI −0,088 és −0,074; NO2 3,6·10⁻⁵ és 3,1·10⁻⁵ mol/m² (4, illetve 3 érvényes megfigyelésből).

## 3. A Python-szkriptek kimenetei

### 3.1. `racs_anomalia.csv` (`kombinal_cellak.py`)

Egy sor egy rácscella. A tél × hónap cellákból az epizódnapok számával súlyozott anomália.

| Oszlop | Jelentés | Egység |
| --- | --- | --- |
| `cell_id` | a cella kulcsa (a rácsgeometriához kapcsolható) | |
| `d_pm` | a CAMS PM2,5 anomália (epizód − alapvonal) | µg/m³ |
| `d_uvai` | az UVAI-anomália | dimenziótlan |
| `d_no2` | az NO2-anomália | mol/m² |
| `cov_pm`, `cov_uvai`, `cov_no2` | az epizódnapok mekkora hányadát fedték érvényes cellák (0–1); 0,5 alatt a `d_*` üres | arány |
| `fb_pm`, `fb_uvai`, `fb_no2` | a tél–hónap cellák mekkora hányadában kellett tartalék-alapvonal | arány |

Az alapvonal az ugyanazon tél ugyanazon hónapjának nem epizódnapi átlaga; ha abban kevés a nap vagy az érvényes megfigyelés (`MIN_NON_DAYS = 5`; S5P: `MIN_OBS_NON`), a hónap telekre összevont átlaga a tartalék.

### 3.2. `kalibracio_eredmenyek.csv` (`kalibracio_epizod.py`)

Egy sor egy epizód-definíció (szigorú, fő, laza).

| Oszlop | Jelentés |
| --- | --- |
| `epizod` | a független epizódok száma |
| `nap` | az epizódnapok száma az értékelésben |
| `tel` | hány télre esnek epizódok |
| `max_telresz` | a legnagyobb tél aránya az epizódnapokból |
| `tobblet` | a PM2,5-többlet: az epizódonkénti átlagos anomália mínusz a nem epizódnapok átlaga (µg/m³) |
| `ci_lo`, `ci_hi` | a többlet 95%-os konfidenciaintervalluma (epizódszintű bootstrap) |
| `pozitiv` | a pozitív epizódok aránya (pl. `26/29`) |
| `p_elojel` | az előjelpróba p-értéke |
| `lefedes` | a legszennyezettebb napok (az anomália felső 10%-a) mekkora részét fogja meg a definíció |
| `min_LOWO` | a nap-súlyozott különbség minimuma, ha egy telet kihagyunk |
| `tartalek_arany` | az epizódos állomásnapok hányada, ahol tartalék-alapvonal kellett |
| `epizodnap_ERA5` | az epizódnapok száma a teljes ERA5-táblán (az értékelt hónapok előtt) |

### 3.3. A kibocsátás-csatolás kimenetei (`emisszio_csatolas.py`)

**`allomas_emisszio_pufferek.csv`** (egy sor egy állomás): az állomástábla oszlopai, plusz minden `R` sugárra (1, 3, 5, 10, 20 km):

| Oszlop | Jelentés | Egység |
| --- | --- | --- |
| `em_kg_Rkm` | az `R` km sugarú pufferbe eső éves kibocsátás | kg |
| `em_dens_Rkm` | az előző a teljes kör területére vetítve (π·R²) | kg/km² |
| `belter_Rkm` | a puffer belterületi poligonjainak területe | km² |

**`racs_emisszio.csv` / `.gpkg`** (egy sor egy rácscella): `em_kg` (a cellára eső kibocsátás, kg), `belter_km2` (a belterület a cellában), `em_dens` (kg/km² a cella területére), `belter_arany` (a cella belterületi aránya).

**`nem_illesztett_tabla.csv`**: `tabla_kulcs`, `kibocsatas_kg` (a nem illesztett település kibocsátása), `javaslat` (hasonló nevű poligon). **`nem_illesztett_poligon.csv`**: `poligon_kulcs`.

A számolás kibocsátásként az `EMISSZIO` beállítás szerinti oszlopot használja (alapértelmezés: PM2,5, HU készülékmegoszlás szerint).

## 4. Előkészített bemenetek

### 4.1. `pm25_napi_hosszu.csv`

Egy sor egy állomás egy napi PM2,5-értéke (csak a nem hiányzó mérések).

| Oszlop | Jelentés |
| --- | --- |
| `datum` | a nap dátuma |
| `tel` | a fűtési szezon (pl. `2024/2025`) |
| `allomas_id` | az állomás azonosítója (`allomasok_tisztitott.csv` `id` oszlopa) |
| `pm_oszlop` | az állomás oszlopneve az eredeti OLM-fájlokban |
| `pm25` | a napi PM2,5-koncentráció (µg/m³) |

### 4.2. `allomasok_tisztitott.csv`

Egy sor egy állomás. A feltöltött változatban 37 sor van: 36 PM2,5-öt mérő állomás és egy mérés nélküli (Budapest, Kosztolányi Dezső tér).

| Oszlop | Jelentés |
| --- | --- |
| `id` | állomásazonosító (az eredeti 52 soros lista sorszáma) |
| `varos`, `cím` | a település és a cím |
| `Állomás típusa` | az OLM szerinti besorolás (városi/külvárosi/vidéki háttér, közlekedési, ipari) |
| `lon`, `lat` | WGS84-koordináták |
| `geokodolt` | `True`: az eredeti geokódolásból (17 állomás); `False`: utólag kiegészített koordináta (20) |
| `pm_oszlop` | az oszlopnév a PM-fájlokban (az illesztés kulcsa) |
| `tipus_csoport` | egyszerűsített csoport: háttér (23), közlekedési (12), ipari (2) |
| `budapest` | 1, ha budapesti (10 sor, ebből 9 mér PM2,5-öt) |
| `lefed_2021_2022` … `lefed_2025_2026` | a napok %-a, amikor volt mérés az adott télen (183 napból) |
| `napok_adattal` | a mérési napok száma összesen (legfeljebb 915 = 5 × 183) |
| `pm25_atlag` | az összes mért nap átlaga (µg/m³); az eltérő hiányzó időszakok miatt állomások között nem közvetlenül összehasonlítható |
| `pm25_van` | 1, ha van PM2,5-mérés |

Néhány tudnivaló: 21 állomás éri el mind az öt télen a 75%-os lefedettséget; Szentgotthárd csak a 2025/26-os télben mér, Szolnok mérése a 2025/26-os télre megszűnt.

## 5. Ismert tudnivalók és ellenőrizendők

- **Mértékegység:** a GEE-ből a CAMS PM2,5 kg/m³-ben érkezik; µg/m³-hez 10⁹-cel kell szorozni (a `kombinal_cellak.py` ezt elvégzi).
- **Hiányzó NO2:** a 2025/26-os téli fájlban az epizódnapi NO2 a cellák 43%-ában (november) és 31%-ában (december) hiányzik (felhőzet); január és február teljes.
- **CAMS-felbontás:** a CAMS értékei havonta csak 98 különböző értéket vesznek fel, vagyis ennyi független CAMS-pixel fedi az országot.
- **Epizódnapok száma:** az ERA5-exportból számolva 137 nap adódik
