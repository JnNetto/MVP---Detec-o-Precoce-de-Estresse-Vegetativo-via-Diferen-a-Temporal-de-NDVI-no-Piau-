"""
ndvi_sentinel2.py
=================
1. Busca as datas disponiveis na API SATVeg/Embrapa para o bbox do ponto
2. Para cada data MODIS disponivel, verifica se existe cena Sentinel-2
   com baixa cobertura de nuvens
3. Baixa a cena, calcula NDVI, gera imagem RGB de cor natural
4. Salva tudo em output_ndvi/serie_sentinel2.json
5. Para ao atingir meta_pares entradas

Credenciais (credentials.json):
    {
        "client_id":       "email Copernicus",
        "client_secret":   "senha Copernicus",
        "consumer_key":    "AgroAPI consumer key",
        "consumer_secret": "AgroAPI consumer secret"
    }
"""

import json, os, sys, zipfile
from datetime import datetime, timedelta

if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

import numpy as np
import requests
import rasterio
from rasterio.enums import Resampling
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

CREDENTIALS_FILE  = "credentials.json"

class TokenExpiredError(Exception):
    pass
TOKEN_URL         = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
CATALOG_URL       = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
DOWNLOAD_URL      = "https://zipper.dataspace.copernicus.eu/odata/v1/Products"
AGROAPI_TOKEN_URL = "https://api.cnptia.embrapa.br/token"
SATVEG_SERIE_URL  = "https://api.cnptia.embrapa.br/satveg/v2/seriespoligono"


# ---------------------------------------------------------------------------
# Credenciais
# ---------------------------------------------------------------------------

def load_credentials(path=CREDENTIALS_FILE):
    if not os.path.exists(path):
        raise FileNotFoundError(f"'{path}' nao encontrado.")
    with open(path, encoding="utf-8") as f:
        return json.load(f)

def get_copernicus_token(client_id, client_secret):
    print("[AUTH-COPERNICUS] Obtendo token...")
    r = requests.post(TOKEN_URL, data={
        "grant_type": "password",
        "username": client_id,
        "password": client_secret,
        "client_id": "cdse-public",
    }, headers={"Content-Type": "application/x-www-form-urlencoded"})
    if r.status_code != 200:
        raise RuntimeError(f"Falha: {r.status_code} - {r.text[:200]}")
    print("[AUTH-COPERNICUS] OK.")
    return r.json()["access_token"]

def get_agroapi_token(consumer_key, consumer_secret):
    print("[AUTH-AGROAPI] Obtendo token...")
    r = requests.post(
        AGROAPI_TOKEN_URL,
        data={"grant_type": "client_credentials"},
        auth=(consumer_key, consumer_secret),
        headers={"Content-Type": "application/x-www-form-urlencoded"},
    )
    if r.status_code != 200:
        raise RuntimeError(f"Falha AgroAPI: {r.status_code} - {r.text[:200]}")
    print("[AUTH-AGROAPI] OK.")
    return r.json()["access_token"]


# ---------------------------------------------------------------------------
# API SATVeg — busca datas disponíveis para o bbox
# ---------------------------------------------------------------------------

def bbox_to_poligono(bbox):
    """[minLon, minLat, maxLon, maxLat] -> string 'lon lat,...' para API SATVeg."""
    min_lon, min_lat, max_lon, max_lat = bbox
    pts = [
        (min_lon, min_lat),
        (max_lon, min_lat),
        (max_lon, max_lat),
        (min_lon, max_lat),
        (min_lon, min_lat),
    ]
    return ",".join(f"{lon} {lat}" for lon, lat in pts)

def fetch_modis_dates(agroapi_token, bbox):
    """
    Busca a serie MODIS completa para o bbox via API SATVeg.
    Retorna lista de datas ISO disponíveis ordenadas (mais recente primeiro).
    """
    poligono = bbox_to_poligono(bbox)
    payload = {
        "tipoPerfil": "ndvi",
        "poligono":   poligono,
        "satelite":   "comb",
        "preFiltro":  3,
    }
    print(f"[MODIS] Buscando datas disponiveis na API SATVeg...")
    r = requests.post(
        SATVEG_SERIE_URL,
        json=payload,
        headers={
            "Authorization": f"Bearer {agroapi_token}",
            "Content-Type":  "application/json",
        },
    )
    if r.status_code != 200:
        raise RuntimeError(f"Erro SATVeg: {r.status_code} - {r.text[:300]}")

    data    = r.json()
    datas   = data.get("listaDatas", [])
    valores = data.get("listaSerie", [])

    # Filtra datas com valor valido (nao None, nao negativo) e a partir de 2023
    ano_minimo = 2023
    pares = []
    for d, v in zip(datas, valores):
        if v is not None and float(v) > 0 and d[:4] >= str(ano_minimo):
            pares.append((d, float(v)))

    # Ordena mais recente primeiro
    pares.sort(key=lambda x: x[0], reverse=True)
    print(f"[MODIS] {len(pares)} datas validas encontradas.")
    return pares  # lista de (data_iso, ndvi_modis)


# ---------------------------------------------------------------------------
# Catálogo Sentinel-2
# ---------------------------------------------------------------------------

def build_bbox_wkt(bbox):
    a, b, c, d = bbox
    return f"POLYGON(({a} {b},{c} {b},{c} {d},{a} {d},{a} {b}))"

def search_sentinel_date(cop_token, bbox, date, max_cloud):
    """Busca cena Sentinel-2 L2A para uma data exata."""
    wkt = build_bbox_wkt(bbox)
    params = {
        "$filter": (
            f"Collection/Name eq 'SENTINEL-2' "
            f"and OData.CSC.Intersects(area=geography'SRID=4326;{wkt}') "
            f"and ContentDate/Start gt {date}T00:00:00Z "
            f"and ContentDate/Start lt {date}T23:59:59Z "
            f"and Attributes/OData.CSC.DoubleAttribute/any(att:att/Name eq 'cloudCover' "
            f"and att/OData.CSC.DoubleAttribute/Value le {max_cloud}) "
            f"and contains(Name,'L2A')"
        ),
        "$orderby": "ContentDate/Start asc",
        "$top": "3",
        "$expand": "Attributes",
    }
    r = requests.get(CATALOG_URL, params=params, headers={"Authorization": f"Bearer {cop_token}"})
    if r.status_code == 403:
        raise TokenExpiredError("Token Copernicus expirado (403)")
    if r.status_code != 200:
        raise RuntimeError(f"Erro catalogo: {r.status_code}")
    return r.json().get("value", [])

def search_sentinel_window(cop_token, bbox, date_iso, window_days, max_cloud):
    """
    Busca cena Sentinel-2 numa janela de +-window_days em torno da data MODIS.
    Retorna o produto mais próximo encontrado, ou None.
    """
    target = datetime.strptime(date_iso, "%Y-%m-%d")
    for delta in range(0, window_days + 1):
        for sign in ([0] if delta == 0 else [1, -1]):
            candidate = (target + timedelta(days=delta * sign)).strftime("%Y-%m-%d")
            products = search_sentinel_date(cop_token, bbox, candidate, max_cloud)
            if products:
                return products[0], candidate
    return None, None


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

def download_product(cop_token, product_id, product_name, out_dir):
    zip_path = os.path.join(out_dir, f"{product_name}.zip")
    if os.path.exists(zip_path):
        print(f"  [DOWNLOAD] Cache encontrado.")
        return zip_path
    url = f"{DOWNLOAD_URL}({product_id})/$value"
    print(f"  [DOWNLOAD] Baixando (aguarde)...")
    with requests.get(url, headers={"Authorization": f"Bearer {cop_token}"}, stream=True) as r:
        if r.status_code == 403:
            raise TokenExpiredError("Token Copernicus expirado durante download (403)")
        if r.status_code != 200:
            raise RuntimeError(f"Erro download: {r.status_code}")
        with open(zip_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=4*1024*1024):
                f.write(chunk)
    print(f"  [DOWNLOAD] Concluido.")
    return zip_path


# ---------------------------------------------------------------------------
# Extração de bandas (NDVI + RGB)
# ---------------------------------------------------------------------------

def extract_bands(zip_path, out_dir, resolution="medium"):
    """Extrai B02, B03, B04 (RGB) e B8A/B08 (NIR) do ZIP."""
    res_map = {"high": "R10m", "medium": "R20m", "low": "R60m"}
    target_res = res_map.get(resolution, "R20m")
    bands = {}
    with zipfile.ZipFile(zip_path, "r") as zf:
        all_files = zf.namelist()

        def find_band(codes, res):
            for name in all_files:
                if f"IMG_DATA/{res}/" in name and name.endswith(".jp2"):
                    for code in codes:
                        if f"_{code}_" in name or name.endswith(f"_{code}.jp2"):
                            return name
            return None

        red_file  = find_band(["B04"], target_res) or find_band(["B04"], "R10m")
        nir_file  = find_band(["B08"], "R10m") if target_res == "R10m" else (find_band(["B8A"], target_res) or find_band(["B08"], "R10m"))
        blue_file = find_band(["B02"], target_res) or find_band(["B02"], "R10m")
        green_file= find_band(["B03"], target_res) or find_band(["B03"], "R10m")
        # SCL: Scene Classification Layer — classifica nuvens pixel a pixel
        scl_file  = find_band(["SCL"], "R20m") or find_band(["SCL"], "R60m")

        if not red_file or not nir_file:
            raise RuntimeError("Bandas B04/B08 nao encontradas no ZIP.")

        for key, fname in [
            ("red",   red_file),
            ("nir",   nir_file),
            ("blue",  blue_file),
            ("green", green_file),
            ("scl",   scl_file),
        ]:
            if fname is None:
                continue
            dest = os.path.join(out_dir, os.path.basename(fname))
            if not os.path.exists(dest):
                with zf.open(fname) as src, open(dest, "wb") as dst:
                    dst.write(src.read())
            bands[key] = dest
    return bands


# ---------------------------------------------------------------------------
# Verificação de nuvens no bbox via SCL
# ---------------------------------------------------------------------------

def check_cloud_bbox(scl_path, bbox, max_cloud_pct=30):
    """
    Verifica a porcentagem de pixels nublados DENTRO do bbox usando a SCL.
    Classes SCL consideradas nuvem/invalido:
      0  = sem dados
      1  = pixels saturados
      3  = sombra de nuvem
      8  = nuvem media probabilidade
      9  = nuvem alta probabilidade
      10 = cirrus
    Retorna (cloud_pct, aceita) onde aceita=True se cloud_pct <= max_cloud_pct.
    """
    if scl_path is None or not os.path.exists(scl_path):
        return None, True  # sem SCL, aceita sem verificar

    from rasterio.windows import from_bounds
    from pyproj import Transformer

    CLASSES_NUVEM = {0, 1, 3, 8, 9, 10}

    with rasterio.open(scl_path) as src:
        min_lon, min_lat, max_lon, max_lat = bbox
        try:
            transformer = Transformer.from_crs("EPSG:4326", src.crs.to_epsg(), always_xy=True)
            x_min, y_min = transformer.transform(min_lon, min_lat)
            x_max, y_max = transformer.transform(max_lon, max_lat)
            win = from_bounds(x_min, y_min, x_max, y_max, src.transform)
            col_off = max(0, int(win.col_off))
            row_off = max(0, int(win.row_off))
            width   = min(int(win.width),  src.width  - col_off)
            height  = min(int(win.height), src.height - row_off)
            if width <= 0 or height <= 0:
                return None, True
            from rasterio.windows import Window
            safe_win = Window(col_off, row_off, width, height)
            scl = src.read(1, window=safe_win)
        except Exception as e:
            return None, True

    total = scl.size
    nuvem = sum(np.sum(scl == c) for c in CLASSES_NUVEM)
    cloud_pct = nuvem / max(total, 1) * 100
    aceita = cloud_pct <= max_cloud_pct
    return cloud_pct, aceita


# ---------------------------------------------------------------------------
# NDVI
# ---------------------------------------------------------------------------

def calculate_ndvi(red_path, nir_path, out_dir, date_label):
    with rasterio.open(red_path) as rs:
        red = rs.read(1).astype(np.float32)
        profile = rs.profile.copy()
    with rasterio.open(nir_path) as ns:
        if ns.shape != (profile["height"], profile["width"]):
            nir = ns.read(1, out_shape=(profile["height"], profile["width"]),
                          resampling=Resampling.bilinear).astype(np.float32)
        else:
            nir = ns.read(1).astype(np.float32)
    np.seterr(divide="ignore", invalid="ignore")
    # Pixels com red=0 E nir=0 sao fora do swath (sem dados), nao solo exposto
    fora_do_swath = (red == 0) & (nir == 0)
    ndvi = np.where((nir + red) == 0, np.nan, (nir - red) / (nir + red))
    ndvi[fora_do_swath] = np.nan  # garante que fora do swath seja NaN
    profile.update(driver="GTiff", dtype=rasterio.float32, count=1, compress="lzw", nodata=np.nan)
    tif_path = os.path.join(out_dir, f"ndvi_{date_label}.tif")
    with rasterio.open(tif_path, "w", **profile) as dst:
        dst.write(ndvi.astype(np.float32), 1)
    return ndvi, tif_path

def compute_stats(ndvi):
    valid = ndvi[~np.isnan(ndvi)]
    return {
        "mean":   float(np.nanmean(valid)),
        "median": float(np.nanmedian(valid)),
        "min":    float(np.nanmin(valid)),
        "max":    float(np.nanmax(valid)),
        "std":    float(np.nanstd(valid)),
        "pixels_valid": int(valid.size),
        "pixels_total": int(ndvi.size),
    }


# ---------------------------------------------------------------------------
# Imagem RGB de cor natural
# ---------------------------------------------------------------------------

def save_rgb(bands, out_dir, date_label, bbox=None):
    """
    Gera imagem RGB de cor natural (B04=R, B03=G, B02=B) recortada para o bbox.
    Salva como PNG em output_ndvi/rgb_YYYYMMDD.png
    Retorna None se a area for toda fora do swath.
    """
    if not all(k in bands for k in ("red", "green", "blue")):
        print("  [RGB] Bandas RGB nao disponiveis, pulando.")
        return None

    from rasterio.windows import from_bounds, Window

    def read_cropped(path, win):
        with rasterio.open(path) as src:
            # Intersecta a window com os limites reais do arquivo
            col_off = max(0, int(win.col_off))
            row_off = max(0, int(win.row_off))
            width   = min(int(win.width),  src.width  - col_off)
            height  = min(int(win.height), src.height - row_off)
            if width <= 0 or height <= 0:
                return np.zeros((10, 10), dtype=np.float32)
            safe_win = Window(col_off, row_off, width, height)
            return src.read(1, window=safe_win).astype(np.float32)

    # Calcula window a partir do bbox — reprojetando para o CRS do arquivo (UTM)
    with rasterio.open(bands["red"]) as rs:
        if bbox:
            min_lon, min_lat, max_lon, max_lat = bbox
            # Sentinel-2 L2A JP2 esta em UTM — precisa reprojetar o bbox de WGS84 para UTM
            from pyproj import Transformer
            try:
                transformer = Transformer.from_crs("EPSG:4326", rs.crs.to_epsg(), always_xy=True)
                x_min, y_min = transformer.transform(min_lon, min_lat)
                x_max, y_max = transformer.transform(max_lon, max_lat)
                win = from_bounds(x_min, y_min, x_max, y_max, rs.transform)
                print(f"  [RGB] CRS arquivo: {rs.crs}, bbox UTM: ({x_min:.0f},{y_min:.0f})->({x_max:.0f},{y_max:.0f})")
            except Exception as e:
                print(f"  [RGB] Falha na reprojecao ({e}), tentando WGS84 direto...")
                win = from_bounds(min_lon, min_lat, max_lon, max_lat, rs.transform)
        else:
            win = Window(0, 0, rs.width, rs.height)

    r = read_cropped(bands["red"],   win)
    g = read_cropped(bands["green"], win)
    b = read_cropped(bands["blue"],  win)

    # Reamostrar g e b para o shape de r se divergirem (sem dependencias extras)
    def resample(arr, target_shape):
        if arr.shape == target_shape:
            return arr
        from rasterio.transform import from_bounds as tfrom_bounds
        zoom_r = target_shape[0] / arr.shape[0]
        zoom_c = target_shape[1] / arr.shape[1]
        row_idx = np.clip((np.arange(target_shape[0]) / zoom_r).astype(int), 0, arr.shape[0]-1)
        col_idx = np.clip((np.arange(target_shape[1]) / zoom_c).astype(int), 0, arr.shape[1]-1)
        return arr[np.ix_(row_idx, col_idx)]

    g = resample(g, r.shape)
    b = resample(b, r.shape)

    # Detecta pixels fora do swath (as tres bandas == 0)
    fora = (r == 0) & (g == 0) & (b == 0)
    pct_valido = (~fora).sum() / max(fora.size, 1) * 100
    print(f"  [RGB] {pct_valido:.1f}% pixels validos no recorte.")

    if pct_valido < 5:
        print(f"  [RGB] Area totalmente fora do swath — imagem nao gerada.")
        return None

    def normalize(arr, mascara_invalida):
        vals = arr[~mascara_invalida]
        if len(vals) == 0:
            return np.zeros_like(arr, dtype=np.float32)
        p2  = np.percentile(vals, 2)
        p98 = np.percentile(vals, 98)
        out = np.clip((arr - p2) / max(p98 - p2, 1e-6), 0, 1)
        out[mascara_invalida] = 0.82  # cinza claro para fora do swath
        return out.astype(np.float32)

    rgb = np.dstack([normalize(r, fora), normalize(g, fora), normalize(b, fora)])

    date_fmt = f"{date_label[:4]}-{date_label[4:6]}-{date_label[6:]}"
    png_path = os.path.join(out_dir, f"rgb_{date_label}.png")
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.imshow(rgb, interpolation="nearest")
    ax.set_title(f"Cor natural (RGB) Sentinel-2 {date_fmt}\n{pct_valido:.0f}% pixels validos no bbox")
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  [RGB] Imagem salva: {png_path}")
    return png_path


# ---------------------------------------------------------------------------
# Série JSON
# ---------------------------------------------------------------------------

def load_serie(out_dir):
    path = os.path.join(out_dir, "serie_sentinel2.json")
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    return {}

def save_serie(out_dir, serie):
    path = os.path.join(out_dir, "serie_sentinel2.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(serie, f, indent=2, ensure_ascii=False)

def datas_ja_coletadas(serie, ponto_id):
    return {e["data"] for e in serie.get(ponto_id, [])}

def datas_modis_ja_coletadas(serie, ponto_id):
    """Retorna set de datas MODIS que ja tem par no JSON."""
    return {e["data_modis"] for e in serie.get(ponto_id, []) if "data_modis" in e}

def total_coletado(serie, ponto_id):
    return len(serie.get(ponto_id, []))


# ---------------------------------------------------------------------------
# Loop principal
# ---------------------------------------------------------------------------

def run(config, credentials_path=CREDENTIALS_FILE, out_dir="output_ndvi"):
    os.makedirs(out_dir, exist_ok=True)

    ponto_id   = config.get("ponto_id", "ponto_1")
    bbox       = config["bbox"]
    max_cloud  = config.get("maxCloudCover", 25)
    resolution = config.get("resolution", "medium")
    meta       = config.get("meta_pares", 12)
    # Janela em dias para buscar Sentinel-2 em torno da data MODIS
    # 0 = data exata apenas; 1 = +/-1 dia; etc.
    window     = config.get("janela_dias", 1)
    # Gerar RGB apenas para a primeira cena nova (evita processar todas)
    rgb_apenas_primeira = config.get("rgb_apenas_primeira", True)

    creds = load_credentials(credentials_path)
    cop_token    = get_copernicus_token(creds["client_id"], creds["client_secret"])
    agroapi_token = get_agroapi_token(creds["consumer_key"], creds["consumer_secret"])

    # Busca datas MODIS disponíveis para o bbox
    modis_pares = fetch_modis_dates(agroapi_token, bbox)

    serie = load_serie(out_dir)
    ja_coletadas       = datas_ja_coletadas(serie, ponto_id)
    ja_coletadas_modis = datas_modis_ja_coletadas(serie, ponto_id)
    atual = total_coletado(serie, ponto_id)

    print(f"\n[LOOP] Ponto         : {ponto_id}")
    print(f"[LOOP] Meta          : {meta} entradas")
    print(f"[LOOP] Ja coletadas  : {atual}")
    print(f"[LOOP] Faltam        : {max(0, meta - atual)}")
    print(f"[LOOP] Datas MODIS   : {len(modis_pares)} disponiveis")
    print(f"[LOOP] Janela busca  : +-{window} dias em torno de cada data MODIS\n")

    if atual >= meta:
        print("[LOOP] Meta ja atingida. Rode analise_correlacao.py")
        return

    coletados_agora = 0
    rgb_gerado = False

    # --- Balanceamento sazonal ---
    # Meses obrigatorios: pelo menos 1 entrada por mes
    # Se algum mes nao tiver dados, prioriza meses chuvosos como fallback
    MESES_CHUVOSOS  = config.get("meses_chuvosos",  [2, 3, 4])
    MESES_SECOS     = config.get("meses_secos",     [8, 9, 10])
    TODOS_MESES     = list(range(1, 13))

    def meses_ja_cobertos(serie, ponto_id):
        return {datetime.strptime(e["data"], "%Y-%m-%d").month
                for e in serie.get(ponto_id, [])}

    def mes_prioritario(data_modis_str, serie, ponto_id):
        """
        Retorna True se a data MODIS deve ser tentada agora.
        Logica:
          1. Se o mes ainda nao tem nenhuma entrada -> prioritario
          2. Se todos os meses ja tem pelo menos 1 entrada -> qualquer mes serve
          3. Se ainda faltam meses sem cobertura -> so tenta meses sem cobertura
             (priorizando chuvosos primeiro)
        """
        mes = datetime.strptime(data_modis_str, "%Y-%m-%d").month
        cobertos = meses_ja_cobertos(serie, ponto_id)
        faltando = [m for m in TODOS_MESES if m not in cobertos]

        if not faltando:
            # Todos os meses cobertos, aceita qualquer data
            return True
        if mes in faltando:
            # Esse mes ainda nao tem entrada, prioritario
            return True
        # Esse mes ja tem entrada, mas ainda faltam outros — pula por enquanto
        return False

    # Reordena modis_pares: primeiro datas de meses sem cobertura
    # (chuvosos antes de secos), depois o resto
    def prioridade_mes(par):
        mes = datetime.strptime(par[0], "%Y-%m-%d").month
        serie_atual = load_serie(out_dir)
        cobertos = meses_ja_cobertos(serie_atual, ponto_id)
        faltando = [m for m in TODOS_MESES if m not in cobertos]
        if mes in faltando:
            if mes in MESES_CHUVOSOS:
                return 0  # maxima prioridade
            elif mes in MESES_SECOS:
                return 1
            else:
                return 2
        return 3  # mes ja coberto

    modis_pares_ord = sorted(modis_pares, key=prioridade_mes)

    for data_modis, ndvi_modis in modis_pares_ord:
        if total_coletado(load_serie(out_dir), ponto_id) >= meta:
            break

        # Pula se a data MODIS ja tem par no JSON (verificacao rapida, sem chamar API)
        if data_modis in ja_coletadas_modis:
            continue

        # Pula datas Sentinel ja coletadas
        if data_modis in ja_coletadas:
            continue

        # Verifica balanceamento sazonal
        serie_atual = load_serie(out_dir)
        if not mes_prioritario(data_modis, serie_atual, ponto_id):
            continue

        mes_nome = datetime.strptime(data_modis, "%Y-%m-%d").strftime("%b")
        print(f"[LOOP] Data MODIS {data_modis} ({mes_nome}, NDVI={ndvi_modis:.4f})...", end=" ", flush=True)

        # Busca Sentinel-2 na janela — renova token se expirado
        try:
            product, data_sentinel = search_sentinel_window(
                cop_token, bbox, data_modis, window, max_cloud
            )
        except TokenExpiredError:
            print("token expirado, renovando...")
            cop_token = get_copernicus_token(creds["client_id"], creds["client_secret"])
            product, data_sentinel = search_sentinel_window(
                cop_token, bbox, data_modis, window, max_cloud
            )

        if product is None:
            print(f"sem cena Sentinel-2 em +-{window} dias.")
            continue

        product_id   = product["Id"]
        product_name = product["Name"]
        product_date = product["ContentDate"]["Start"][:10]
        cloud_attr   = next((a["Value"] for a in product.get("Attributes", []) if a["Name"] == "cloudCover"), "N/D")
        cloud_str    = f"{cloud_attr:.2f}%" if isinstance(cloud_attr, float) else str(cloud_attr)
        lag_dias     = abs((datetime.strptime(product_date, "%Y-%m-%d") - datetime.strptime(data_modis, "%Y-%m-%d")).days)

        print(f"ENCONTRADO: {product_date} (lag={lag_dias}d, nuvens={cloud_str})")

        try:
            zip_path = download_product(cop_token, product_id, product_name, out_dir)
            bands    = extract_bands(zip_path, out_dir, resolution)
            date_label = product_date.replace("-", "")

            # Verifica nuvens NO BBOX via SCL antes de calcular NDVI
            max_cloud_bbox = config.get("maxCloudBbox", 20)
            cloud_bbox_pct, bbox_ok = check_cloud_bbox(
                bands.get("scl"), bbox, max_cloud_pct=max_cloud_bbox
            )
            if cloud_bbox_pct is not None:
                print(f"  [SCL] Nuvens no bbox: {cloud_bbox_pct:.1f}%", end="")
                if not bbox_ok:
                    print(f" — acima de {max_cloud_bbox}%, pulando.")
                    continue
                else:
                    print(f" — OK.")

            ndvi, _  = calculate_ndvi(bands["red"], bands["nir"], out_dir, date_label)
            stats    = compute_stats(ndvi)

            # Valida se o bbox tem pixels validos (nao esta fora do swath)
            pixels_validos_pct = stats["pixels_valid"] / max(stats["pixels_total"], 1) * 100
            if stats["pixels_valid"] == 0 or pixels_validos_pct < 20:
                print(f"  [SKIP] {pixels_validos_pct:.1f}% pixels validos — bbox fora do swath, pulando.")
                continue
            if np.isnan(stats["mean"]):
                print(f"  [SKIP] NDVI medio = NaN — bbox sem dados validos, pulando.")
                continue

            # RGB — gera para todas as cenas, ou só primeira conforme config
            if not rgb_apenas_primeira or not rgb_gerado:
                save_rgb(bands, out_dir, date_label, bbox=bbox)
                rgb_gerado = True

        except TokenExpiredError:
            print("  [TOKEN] Token expirado durante processamento, renovando...")
            cop_token = get_copernicus_token(creds["client_id"], creds["client_secret"])
            print("  [TOKEN] Token renovado. Tente rodar novamente — esta data sera reprocessada.")
            continue
        except Exception as e:
            print(f"  [ERRO] {e}")
            continue

        # Salva no JSON com ndvi_modis já preenchido
        serie = load_serie(out_dir)
        if ponto_id not in serie:
            serie[ponto_id] = []
        serie[ponto_id] = [e for e in serie[ponto_id] if e["data"] != product_date]
        serie[ponto_id].append({
            "data":            product_date,
            "data_modis":      data_modis,
            "lag_dias":        lag_dias,
            "bbox":            bbox,
            "ndvi_sentinel2":  round(stats["mean"], 6),
            "nuvens_pct":      round(cloud_attr, 2) if isinstance(cloud_attr, float) else cloud_attr,
            "produto":         product_name,
            "ndvi_modis":      round(ndvi_modis, 6),  # preenchido automaticamente
        })
        serie[ponto_id].sort(key=lambda e: e["data"])
        save_serie(out_dir, serie)

        ja_coletadas.add(product_date)
        ja_coletadas_modis.add(data_modis)
        coletados_agora += 1
        total_agora = total_coletado(serie, ponto_id)

        print(f"  [OK] S2={stats['mean']:.4f} | MODIS={ndvi_modis:.4f} | Total: {total_agora}/{meta}")

    # Resumo
    total_final = total_coletado(load_serie(out_dir), ponto_id)
    print(f"\n{'='*60}")
    print(f"LOOP ENCERRADO")
    print(f"  Coletados agora : {coletados_agora}")
    print(f"  Total no JSON   : {total_final}/{meta}")
    if total_final >= meta:
        print(f"\n  Meta atingida! Rode: python analise_correlacao.py")
    else:
        print(f"\n  Faltam {meta - total_final} entradas.")
        print(f"  Tente aumentar 'janela_dias' ou 'maxCloudCover' no CONFIG.")
    print(f"{'='*60}")


# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------

CONFIG = {
    "ponto_id": "ponto_10",

    # Meta de pares a coletar
    "meta_pares": 10,

    # Bbox do ponto: [minLon, minLat, maxLon, maxLat]
    "bbox": [
-40.75080,
-7.23072,
-40.74565,
-7.22547],

    # Janela em dias para buscar Sentinel-2 em torno de cada data MODIS
    # 0 = data exata; 1 = +/-1 dia para compensar fuso horario
    # Mantenha baixo (0 ou 1) para garantir comparacao temporal valida
    "janela_dias": 3,

    # Se True, gera RGB apenas para a primeira cena nova encontrada
    # Se False, gera RGB para todas as cenas (mais lento)
    "rgb_apenas_primeira": False,

    "maxCloudCover": 35,       # nuvens maximas no TILE inteiro (filtro do catalogo)
    "maxCloudBbox": 30,        # nuvens maximas NO BBOX especificamente (via SCL)
    "resolution": "medium",

    # Balanceamento sazonal
    # O script garante pelo menos 1 entrada por mes antes de repetir meses
    # Quando faltam meses, prioriza chuvosos > secos > transicao
    "meses_chuvosos": [2, 3, 4],       # fev, mar, abr
    "meses_secos":    [8, 9, 10],      # ago, set, out
}

if __name__ == "__main__":
    if len(sys.argv) > 1:
        with open(sys.argv[1], encoding="utf-8") as f:
            CONFIG = json.load(f)
    run(CONFIG)