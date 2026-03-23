"""
analise_correlacao.py
=====================
1. Busca automaticamente os valores NDVI MODIS via API SATVeg/Embrapa
   para todas as datas coletadas no serie_sentinel2.json
2. Preenche ndvi_modis no JSON automaticamente (sem trabalho manual)
3. Calcula r², MAE, Pearson, ΔNDVI sazonal
4. Gera gráficos comparativos

Credenciais (credentials.json):
    {
        "client_id": "...",        (Copernicus)
        "client_secret": "...",    (Copernicus)
        "consumer_key": "...",     (AgroAPI/SATVeg)
        "consumer_secret": "..."   (AgroAPI/SATVeg)
    }

Rode:  python analise_correlacao.py
"""

import json
import os
import sys
from datetime import datetime, timedelta

if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

import numpy as np
import requests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SERIE_FILE   = os.path.join("output_ndvi", "serie_sentinel2.json")
CREDS_FILE   = "credentials.json"
OUT_DIR      = "output_ndvi"

AGROAPI_TOKEN_URL = "https://api.cnptia.embrapa.br/token"
SATVEG_URL        = "https://api.cnptia.embrapa.br/satveg/v2/seriespoligono"

# Períodos sazonais — ajuste se necessário
PERIODO_CHUVOSO = (2, 4)   # fev–abr
PERIODO_SECO    = (8, 10)  # ago–out


# ---------------------------------------------------------------------------
# Credenciais
# ---------------------------------------------------------------------------

def load_credentials():
    if not os.path.exists(CREDS_FILE):
        raise FileNotFoundError(f"'{CREDS_FILE}' nao encontrado.")
    with open(CREDS_FILE, encoding="utf-8") as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Auth AgroAPI
# ---------------------------------------------------------------------------

def get_agroapi_token(consumer_key, consumer_secret):
    print("[AGROAPI] Obtendo token...")
    r = requests.post(
        AGROAPI_TOKEN_URL,
        data={"grant_type": "client_credentials"},
        auth=(consumer_key, consumer_secret),
        headers={"Content-Type": "application/x-www-form-urlencoded"},
    )
    if r.status_code != 200:
        raise RuntimeError(f"Falha no token AgroAPI: {r.status_code} - {r.text[:300]}")
    print("[AGROAPI] Token obtido.")
    return r.json()["access_token"]


# ---------------------------------------------------------------------------
# Busca série MODIS via API SATVeg — polígono
# ---------------------------------------------------------------------------

def bbox_to_poligono(bbox):
    """
    Converte [minLon, minLat, maxLon, maxLat] para string de polígono
    no formato da API SATVeg: "lat lon,lat lon,..." fechando no primeiro ponto.
    """
    min_lon, min_lat, max_lon, max_lat = bbox
    pts = [
        (min_lat, min_lon),
        (min_lat, max_lon),
        (max_lat, max_lon),
        (max_lat, min_lon),
        (min_lat, min_lon),  # fecha o polígono
    ]
    return ",".join(f"{lon} {lat}" for lat, lon in pts)


def fetch_modis_series(token, bbox, satelite="comb", pre_filtro=3):
    """
    Busca a série temporal NDVI MODIS completa para o polígono do bbox.
    Retorna dict {data_str: ndvi_value} com todas as datas disponíveis.
    """
    poligono = bbox_to_poligono(bbox)
    payload = {
        "tipoPerfil": "ndvi",
        "poligono":   poligono,
        "satelite":   satelite,
        "preFiltro":  pre_filtro,
    }
    print(f"[MODIS] Poligono enviado: {poligono}")
    print(f"[MODIS] Buscando serie NDVI MODIS via API SATVeg...")
    r = requests.post(
        SATVEG_URL,
        json=payload,
        headers={
            "Authorization": f"Bearer {token}",
            "Content-Type":  "application/json",
        },
    )
    if r.status_code != 200:
        raise RuntimeError(f"Erro API SATVeg: {r.status_code} - {r.text[:500]}")

    data = r.json()
    print(f"[DEBUG] Chaves: {list(data.keys())}")

    valores = data.get("listaSerie") or []
    datas   = data.get("listaDatas") or []

    print(f"[DEBUG] len listaSerie={len(valores)}, len listaDatas={len(datas)}")
    if datas:
        print(f"[DEBUG] primeiras datas: {datas[:3]}")
    else:
        print(f"[DEBUG] listaDatas vazia ou ausente — todas as chaves: {list(data.keys())}")
        # Tenta chaves alternativas
        for key in data.keys():
            if "data" in key.lower() or "date" in key.lower():
                print(f"[DEBUG] chave candidata: {key} -> {data[key][:3]}")

    if not valores:
        raise RuntimeError(f"listaSerie vazia na resposta: {str(data)[:300]}")

    # Se nao veio listaDatas, gera datas sinteticas a partir de 2000-02-18
    # com passo de 16 dias (composicao MODIS MOD13Q1)
    if not datas:
        print(f"[DEBUG] Gerando datas sinteticas para {len(valores)} valores MODIS...")
        from datetime import date
        inicio = date(2000, 2, 18)
        datas = [(inicio + timedelta(days=16*i)).strftime("%d/%m/%Y") for i in range(len(valores))]

    serie = {}
    for d, v in zip(datas, valores):
        try:
            # API retorna datas no formato ISO (2000-02-18)
            dt = d if len(d) == 10 and d[4] == "-" else datetime.strptime(d, "%d/%m/%Y").strftime("%Y-%m-%d")
            serie[dt] = float(v) if v is not None else None
        except Exception:
            pass

    print(f"[MODIS] {len(serie)} datas processadas.")
    return serie


def find_closest_modis(modis_serie, target_date_str, max_days=16):
    """
    Encontra o valor MODIS mais próximo de target_date dentro de uma janela
    de max_days dias (padrão 16 = resolução temporal do MODIS).
    """
    target = datetime.strptime(target_date_str, "%Y-%m-%d")
    best_date = None
    best_delta = None

    for d, v in modis_serie.items():
        if v is None:
            continue
        dt = datetime.strptime(d, "%Y-%m-%d")
        delta = abs((dt - target).days)
        if delta <= max_days:
            if best_delta is None or delta < best_delta:
                best_delta = delta
                best_date = d

    if best_date:
        return modis_serie[best_date], best_date
    return None, None


# ---------------------------------------------------------------------------
# Série — leitura e escrita
# ---------------------------------------------------------------------------

def load_serie():
    if not os.path.exists(SERIE_FILE):
        raise FileNotFoundError(
            f"'{SERIE_FILE}' nao encontrado.\n"
            "Rode ndvi_sentinel2.py primeiro."
        )
    with open(SERIE_FILE, encoding="utf-8") as f:
        return json.load(f)

def save_serie(serie):
    with open(SERIE_FILE, "w", encoding="utf-8") as f:
        json.dump(serie, f, indent=2, ensure_ascii=False)


# ---------------------------------------------------------------------------
# Preenchimento automático do ndvi_modis
# ---------------------------------------------------------------------------

def fill_modis_values(serie, modis_serie, ponto_id):
    """
    Preenche automaticamente ndvi_modis em todas as entradas do ponto
    que ainda estão com null, usando o valor MODIS mais próximo.
    """
    entries = serie.get(ponto_id, [])
    preenchidos = 0
    nao_encontrados = []

    for entry in entries:
        if entry.get("ndvi_modis") is not None:
            continue  # já preenchido, não sobrescreve

        v, data_modis = find_closest_modis(modis_serie, entry["data"])
        if v is not None:
            entry["ndvi_modis"] = round(v, 6)
            entry["data_modis_referencia"] = data_modis
            preenchidos += 1
        else:
            nao_encontrados.append(entry["data"])

    serie[ponto_id] = entries

    if preenchidos:
        print(f"[MODIS] {preenchidos} entradas preenchidas automaticamente.")
    if nao_encontrados:
        print(f"[MODIS] Sem correspondencia MODIS para: {nao_encontrados}")

    return serie


# ---------------------------------------------------------------------------
# Métricas
# ---------------------------------------------------------------------------

def r_squared(y_true, y_pred):
    y_true = np.array(y_true, dtype=float)
    y_pred = np.array(y_pred, dtype=float)
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    return float(1 - ss_res / ss_tot) if ss_tot != 0 else float("nan")

def mae(y_true, y_pred):
    return float(np.mean(np.abs(np.array(y_true) - np.array(y_pred))))

def pearson(y_true, y_pred):
    return float(np.corrcoef(y_true, y_pred)[0, 1])

def delta_ndvi(entries, fonte):
    def mes(e):
        return datetime.strptime(e["data"], "%Y-%m-%d").month
    chuvosos = [e[fonte] for e in entries if PERIODO_CHUVOSO[0] <= mes(e) <= PERIODO_CHUVOSO[1] and e[fonte] is not None]
    secos    = [e[fonte] for e in entries if PERIODO_SECO[0]    <= mes(e) <= PERIODO_SECO[1]    and e[fonte] is not None]
    if not chuvosos or not secos:
        return None, None, None
    return (
        round(float(np.mean(chuvosos)) - float(np.mean(secos)), 4),
        round(float(np.mean(chuvosos)), 4),
        round(float(np.mean(secos)),    4),
    )


# ---------------------------------------------------------------------------
# Visualizações
# ---------------------------------------------------------------------------

def plot_series(ponto_id, entries, out_dir):
    datas   = [datetime.strptime(e["data"], "%Y-%m-%d") for e in entries]
    s2      = [e["ndvi_sentinel2"] for e in entries]
    modis   = [e["ndvi_modis"]     for e in entries]

    # Normalizado
    s2_arr = np.array(s2, dtype=float)
    m_arr  = np.array(modis, dtype=float)
    s2_norm = (s2_arr - np.nanmin(s2_arr)) / (np.nanmax(s2_arr) - np.nanmin(s2_arr) or 1)
    m_norm  = (m_arr  - np.nanmin(m_arr))  / (np.nanmax(m_arr)  - np.nanmin(m_arr)  or 1)

    fig, axes = plt.subplots(2, 1, figsize=(13, 9))
    fig.suptitle(f"Comparacao NDVI - {ponto_id}", fontsize=12)

    # Faixas sazonais
    for ax in axes:
        for e, dt in zip(entries, datas):
            m = dt.month
            if PERIODO_CHUVOSO[0] <= m <= PERIODO_CHUVOSO[1]:
                ax.axvspan(dt - timedelta(days=4), dt + timedelta(days=4), alpha=0.08, color="blue")
            elif PERIODO_SECO[0] <= m <= PERIODO_SECO[1]:
                ax.axvspan(dt - timedelta(days=4), dt + timedelta(days=4), alpha=0.08, color="orange")

    axes[0].plot(datas, s2,    "o-",  color="#1a9641", label="Sentinel-2 (mediana filtrada)")
    axes[0].plot(datas, modis, "s--", color="#d73027", label="MODIS/SATVeg")
    axes[0].set_ylabel("NDVI")
    axes[0].set_title("Valores absolutos  |  azul=chuvoso  laranja=seco")
    axes[0].legend(); axes[0].grid(True, alpha=0.3)

    axes[1].plot(datas, s2_norm, "o-",  color="#1a9641", label="Sentinel-2 (norm.)")
    axes[1].plot(datas, m_norm,  "s--", color="#d73027", label="MODIS (norm.)")
    axes[1].set_ylabel("NDVI normalizado (0-1)")
    axes[1].set_title("Normalizado - compara padrao temporal")
    axes[1].legend(); axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    png = os.path.join(out_dir, f"correlacao_{ponto_id}.png")
    plt.savefig(png, dpi=150, bbox_inches="tight"); plt.close()
    print(f"[VIZ] Serie: {png}")


def plot_scatter(ponto_id, entries, out_dir):
    s2    = [e["ndvi_sentinel2"] for e in entries]
    modis = [e["ndvi_modis"]     for e in entries]
    r2    = r_squared(modis, s2)
    pcc   = pearson(modis, s2)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(modis, s2, color="#1a9641", edgecolors="#0a5c2a", zorder=3)
    lim = [min(min(modis), min(s2)) - 0.05, max(max(modis), max(s2)) + 0.05]
    ax.plot(lim, lim, "k--", alpha=0.4, label="y = x (ideal)")
    m, b = np.polyfit(modis, s2, 1)
    x_fit = np.linspace(lim[0], lim[1], 100)
    ax.plot(x_fit, m * x_fit + b, color="#d73027", label=f"Regressao (r²={r2:.3f}, r={pcc:.3f})")
    ax.set_xlabel("NDVI MODIS/SATVeg"); ax.set_ylabel("NDVI Sentinel-2")
    ax.set_title(f"Dispersao - {ponto_id}")
    ax.legend(); ax.grid(True, alpha=0.3)
    plt.tight_layout()
    png = os.path.join(out_dir, f"scatter_{ponto_id}.png")
    plt.savefig(png, dpi=150, bbox_inches="tight"); plt.close()
    print(f"[VIZ] Scatter: {png}")


# ---------------------------------------------------------------------------
# Análise por ponto
# ---------------------------------------------------------------------------

def analyze_ponto(ponto_id, entries, out_dir):
    complete = [e for e in entries if e.get("ndvi_modis") is not None]
    print(f"\n{'='*60}")
    print(f"PONTO: {ponto_id}")
    print(f"{'='*60}")
    print(f"  Entradas totais      : {len(entries)}")
    print(f"  Com MODIS preenchido : {len(complete)}")

    if len(complete) < 3:
        print(f"  [AVISO] Menos de 3 pares - r² nao calculavel.")
        return None

    s2    = [e["ndvi_sentinel2"] for e in complete]
    modis = [e["ndvi_modis"]     for e in complete]

    r2    = r_squared(modis, s2)
    mae_v = mae(modis, s2)
    pcc   = pearson(modis, s2)

    print(f"\n  --- Metricas de correlacao ---")
    print(f"  r²      : {r2:.4f}  {'OK (>= 0.75)' if r2 >= 0.75 else 'abaixo do criterio (< 0.75)'}")
    print(f"  Pearson : {pcc:.4f}")
    print(f"  MAE     : {mae_v:.4f}")

    for fonte, label in [("ndvi_sentinel2", "Sentinel-2"), ("ndvi_modis", "MODIS/SATVeg")]:
        delta, med_c, med_s = delta_ndvi(complete, fonte)
        print(f"\n  --- DELTA NDVI {label} ---")
        if delta is not None:
            print(f"  Media chuvoso (fev-abr) : {med_c:.4f}")
            print(f"  Media seco    (ago-out) : {med_s:.4f}")
            print(f"  DELTA NDVI              : {delta:.4f}")
        else:
            print(f"  Sem dados sazonais suficientes.")

    plot_series(ponto_id, complete, out_dir)
    plot_scatter(ponto_id, complete, out_dir)
    return {"r2": r2, "mae": mae_v, "pearson": pcc}


# ---------------------------------------------------------------------------
# Pipeline principal
# ---------------------------------------------------------------------------

def run():
    os.makedirs(OUT_DIR, exist_ok=True)

    creds = load_credentials()
    consumer_key    = creds.get("consumer_key")
    consumer_secret = creds.get("consumer_secret")
    if not consumer_key or not consumer_secret:
        raise ValueError("credentials.json precisa ter 'consumer_key' e 'consumer_secret' da AgroAPI.")

    serie = load_serie()
    token = get_agroapi_token(consumer_key, consumer_secret)

    resultados = {}

    for ponto_id, entries in serie.items():
        if not entries:
            continue

        # Pega o bbox do config do ponto (guardado nas entradas) ou pede manualmente
        # Como o bbox não está no JSON da série, vamos usar o centro das datas
        # para buscar — o usuário precisa ter o bbox no CONFIG do ndvi_sentinel2.py
        # Alternativa: extrair lat/lon do nome do produto (não confiável)
        # SOLUÇÃO: lemos o CONFIG embutido no ndvi_sentinel2.py se existir,
        # senão pedimos que o usuário adicione bbox ao serie_sentinel2.json

        bbox = entries[0].get("bbox")
        if not bbox:
            print(f"\n[AVISO] Ponto '{ponto_id}' nao tem 'bbox' nas entradas.")
            print(f"  Adicione o campo 'bbox' manualmente no serie_sentinel2.json")
            print(f"  Exemplo: adicione  \"bbox\": [-42.73924, -4.48541, -42.73705, -4.48330]")
            print(f"  em cada entrada do ponto '{ponto_id}', ou so na primeira entrada.")
            continue

        # Busca série MODIS completa uma vez por ponto
        try:
            modis_serie = fetch_modis_series(token, bbox)
        except Exception as e:
            print(f"[ERRO] Falha ao buscar MODIS para '{ponto_id}': {e}")
            continue

        # Preenche ndvi_modis automaticamente
        serie = fill_modis_values(serie, modis_serie, ponto_id)
        save_serie(serie)

        # Analisa
        res = analyze_ponto(ponto_id, serie[ponto_id], OUT_DIR)
        if res:
            resultados[ponto_id] = res

    # Resumo geral
    if len(resultados) > 1:
        print(f"\n{'='*60}")
        print("RESUMO GERAL")
        print(f"{'='*60}")
        r2_vals = [v["r2"] for v in resultados.values()]
        print(f"  r² medio : {np.mean(r2_vals):.4f}")
        print(f"  r² min   : {np.min(r2_vals):.4f}")
        print(f"  r² max   : {np.max(r2_vals):.4f}")
        aprov = sum(1 for v in r2_vals if v >= 0.75)
        print(f"  Pontos com r² >= 0.75: {aprov}/{len(r2_vals)}")

    out_path = os.path.join(OUT_DIR, "resultado_correlacao.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(resultados, f, indent=2, ensure_ascii=False)
    print(f"\n[OK] Resultado salvo em: {out_path}")


if __name__ == "__main__":
    run()