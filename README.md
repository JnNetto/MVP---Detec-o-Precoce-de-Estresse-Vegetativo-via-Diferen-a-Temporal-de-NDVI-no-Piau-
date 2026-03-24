# NDVI Sentinel-2 × MODIS/SATVeg — Pipeline de Correlação

Pipeline Python para coleta automática de NDVI via **Sentinel-2 (Copernicus)** e **MODIS (SATVeg/Embrapa)**, com cálculo de correlação temporal entre as duas fontes. Desenvolvido como suporte empírico ao trabalho de pesquisa sobre detecção de estresse vegetativo no semiárido piauiense.

---

## Contexto

O estudo compara o NDVI gerado pelo Sentinel-2 (10–20 m/pixel) com a série histórica do MODIS disponível no SATVeg/Embrapa (250 m/pixel). A hipótese é que os dois sensores apresentam **correlação temporal significativa (Pearson r ≥ 0,75)**, mesmo com offset sistemático nos valores absolutos — esperado pela diferença de resolução espacial entre eles.

---

## Arquivos

```
ndvi_sentinel2.py       # Script 1 — coleta de dados
analise_correlacao.py   # Script 2 — análise e geração de gráficos
credentials.json        # Credenciais (não versionar)
output_ndvi/            # Pasta gerada automaticamente com todos os resultados
```

---

## Pré-requisitos

### Instalação das dependências

```bash
pip install requests rasterio numpy matplotlib pyproj
```

### Contas necessárias

| Serviço | Para que serve | Cadastro |
|---|---|---|
| Copernicus Data Space | Download das cenas Sentinel-2 | https://dataspace.copernicus.eu |
| AgroAPI Embrapa | Série MODIS via API SATVeg | https://www.agroapi.cnptia.embrapa.br/store |

### Arquivo `credentials.json`

Crie na mesma pasta dos scripts:

```json
{
  "client_id":       "seu_email@exemplo.com",
  "client_secret":   "sua_senha_copernicus",
  "consumer_key":    "sua_consumer_key_agroapi",
  "consumer_secret": "sua_consumer_secret_agroapi"
}
```

> **Atenção:** Como esses valores não sçao commitados, você pode se logar nas plataformas e arranja-las ou pedir as credenciais diretamente.

---

## Script 1 — `ndvi_sentinel2.py`

### O que faz

1. Autentica no Copernicus e na AgroAPI simultaneamente
2. Consulta a API SATVeg para obter todas as datas MODIS disponíveis para o bbox do ponto (filtradas a partir de 2023)
3. Para cada data MODIS, busca uma cena Sentinel-2 L2A dentro de uma janela de ±N dias
4. Baixa a cena, extrai as bandas necessárias e verifica nuvens pixel a pixel no bbox via SCL (Scene Classification Layer)
5. Calcula o NDVI médio da área e gera uma imagem RGB de cor natural recortada para o bbox
6. Salva o par (NDVI Sentinel-2 + NDVI MODIS) no arquivo `output_ndvi/serie_sentinel2.json`
7. Para automaticamente ao atingir a meta de pares definida no CONFIG

### Balanceamento sazonal

O script garante cobertura temporal equilibrada ao priorizar a busca por meses ainda sem dados, na ordem: meses chuvosos (fev–abr) → meses secos (ago–out) → demais meses. Só repete um mês já coberto quando todos os outros já têm pelo menos uma entrada.

### Resiliência

- **Token expirado (403):** renova o token do Copernicus automaticamente e continua sem interromper
- **Queda de rede:** tenta novamente até 3 vezes com intervalo de 30 segundos antes de pular a data
- **Área fora do swath:** detecta pixels sem dados (red=0 e nir=0) e descarta a cena
- **Cache de downloads:** ZIPs e bandas JP2 já baixados são reutilizados automaticamente
- **Cache de pares:** ao reiniciar, pula imediatamente datas MODIS que já têm par salvo no JSON

### Configuração (CONFIG)

```python
CONFIG = {
    "ponto_id":   "ponto_1",   # identificador do ponto — use nomes diferentes para pontos distintos
    "meta_pares": 12,           # quantidade de pares a coletar antes de encerrar

    # Bounding box: [minLon, minLat, maxLon, maxLat] em graus decimais WGS84
    # Recomendado: área homogênea de caatinga com pelo menos 25 ha (500×500 m)
    "bbox": [-42.77136, -4.43949, -42.76626, -4.43483],

    "janela_dias": 1,          # janela de busca do Sentinel-2 em torno da data MODIS
                               # 0 = data exata; 1 = ±1 dia (recomendado)

    "rgb_apenas_primeira": False,  # True = gera RGB só da 1ª cena nova (mais rápido)

    "maxCloudCover": 25,       # % máxima de nuvens no tile inteiro (filtro do catálogo)
    "maxCloudBbox":  20,       # % máxima de nuvens dentro do bbox (via SCL, mais preciso)
    "resolution":    "medium", # "high" = 10m | "medium" = 20m | "low" = 60m

    "meses_chuvosos": [2, 3, 4],   # meses de prioridade máxima no balanceamento
    "meses_secos":    [8, 9, 10],  # meses de prioridade secundária
}
```

### Como usar

```bash
python ndvi_sentinel2.py
```

Para usar um CONFIG externo em JSON:

```bash
python ndvi_sentinel2.py meu_ponto.json
```

### Saídas geradas

```
output_ndvi/
├── serie_sentinel2.json          # série temporal acumulada de todos os pares
├── ndvi_YYYYMMDD.tif             # GeoTIFF do NDVI para cada data
└── rgb_YYYYMMDD.png              # imagem RGB de cor natural recortada para o bbox
```

### Estrutura do `serie_sentinel2.json`

```json
{
  "ponto_1": [
    {
      "data":           "2025-08-21",   // data da cena Sentinel-2
      "data_modis":     "2025-08-21",   // data MODIS correspondente
      "lag_dias":       0,              // diferença em dias entre as duas datas
      "bbox":           [-42.77, ...],  // bbox do ponto
      "ndvi_sentinel2": 0.422541,       // NDVI médio do bbox no Sentinel-2
      "nuvens_pct":     0.03,           // cobertura de nuvens no tile inteiro (%)
      "produto":        "S2C_MSIL2A_...",
      "ndvi_modis":     0.7942          // NDVI MODIS da mesma data (preenchido automaticamente)
    }
  ]
}
```

---

## Script 2 — `analise_correlacao.py`

### O que faz

1. Lê o `serie_sentinel2.json` gerado pelo Script 1
2. Para cada ponto, consulta a API SATVeg com o bbox da primeira entrada e busca a série MODIS completa
3. Preenche automaticamente qualquer `ndvi_modis` ainda nulo, encontrando o valor MODIS mais próximo (janela de ±16 dias)
4. Calcula as métricas de correlação: Pearson r, r², MAE e ΔNDVI sazonal
5. Gera dois gráficos por ponto: série temporal comparativa e diagrama de dispersão

### Como usar

```bash
python analise_correlacao.py
```

> Rode somente após ter coletado pelo menos 3 pares com o Script 1. O ideal é ter cobertura dos dois períodos sazonais (chuvoso e seco) antes de interpretar os resultados.

### Saídas geradas

```
output_ndvi/
├── serie_sentinel2.json          # atualizado com ndvi_modis preenchidos
├── correlacao_ponto_1.png        # série temporal absoluta e normalizada
├── scatter_ponto_1.png           # diagrama de dispersão com linha de regressão
└── resultado_correlacao.json     # métricas consolidadas de todos os pontos
```

### Métricas calculadas

| Métrica | Descrição |
|---|---|
| **Pearson r** | Correlação temporal entre as séries — métrica primária |
| **r²** | Concordância absoluta — métrica secundária (penaliza offset entre escalas) |
| **MAE** | Erro médio absoluto em relação ao MODIS |
| **ΔNDVI chuvoso–seco** | Diferença média entre período chuvoso (fev–abr) e seco (ago–out) |

---

## Fluxo completo

```
1. Configurar credentials.json
2. Definir bbox e ponto_id no CONFIG do ndvi_sentinel2.py
3. Rodar ndvi_sentinel2.py  →  coleta pares até atingir meta_pares
4. (opcional) Rodar de novo para outros pontos com ponto_id diferente
5. Rodar analise_correlacao.py  →  gera métricas e gráficos
```

---

## Notas metodológicas

**Offset sistemático esperado:** o MODIS (250 m) suaviza a mistura espectral entre vegetação e solo, produzindo valores de NDVI absolutos mais altos que o Sentinel-2 (10–20 m), que resolve essa mistura pixel a pixel. Por isso o r² tende a ser negativo mesmo quando a correlação temporal é forte. O Pearson r é a métrica adequada para comparar o padrão temporal entre os dois sensores.

**Escolha do bbox:** recomenda-se área de pelo menos 25 ha (500×500 m) com vegetação homogênea e sem estradas, construções ou corpos d'água. Áreas menores que um pixel MODIS (250×250 m = 6,25 ha) tornam a comparação instável.

**Cobertura de nuvens:** o filtro `maxCloudCover` atua no tile inteiro (100×100 km), enquanto `maxCloudBbox` usa a SCL para verificar pixel a pixel dentro do bbox. Use ambos para garantir que a área de interesse esteja efetivamente livre de nuvens.

**Janela temporal:** o MODIS usa composição de 16 dias (melhor pixel da janela), portanto a data reportada não corresponde necessariamente à foto do dia exato. `janela_dias: 1` é suficiente para a maioria dos casos.
