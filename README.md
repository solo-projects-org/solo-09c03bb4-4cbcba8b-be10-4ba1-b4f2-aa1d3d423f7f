# 农情遥感补采分析 API

用于在安排补采前，从地块编号图和各景无云掩膜中找出覆盖率不足的地块。
服务提供 FastAPI + Swagger UI，可直接上传 GeoTIFF 试算。

## 启动

```bash
docker compose up --build -d
```

打开：

- Swagger：<http://localhost:8000/docs>
- 健康检查：<http://localhost:8000/health>

可用 `APP_PORT` 修改宿主机端口，例如：

```bash
APP_PORT=8080 docker compose up --build -d
```

## 接口

### `POST /analyze`

以 `multipart/form-data` 上传：

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `parcel` | 文件，单个 | 地块编号图，单波段 GeoTIFF |
| `masks` | 文件，可重复 | 各景无云掩膜，单波段 GeoTIFF，至少一景 |
| `target_coverage` | 数值 | 目标覆盖率，范围 0 到 1 |

输入约束：

- 全部文件必须是单波段 GeoTIFF，且北向上（无旋转、无剪切）。
- 全部影像使用相同坐标系、相同分辨率。
- 掩膜原点可与地块图原点不同，但偏移必须是整数个像元，保证像元格网重合。
- 地块图中只有有效的正整数是地块编号；0、负数、小数和 nodata 会被忽略。
- 掩膜中有效值 `1` 表示该像元无云；同一像元被多景覆盖时只计一次。

不兼容时返回 HTTP 400，并在 `detail.incompatible_files` 中逐个列出文件名、
角色（`parcel` / `mask`）和原因，例如：

```json
{
  "detail": {
    "message": "输入影像不兼容，未执行分析。",
    "incompatible_files": [
      {
        "filename": "bad_scene.tif",
        "role": "mask",
        "reason": "坐标系与地块图不一致（EPSG:4326 vs 地块图 EPSG:32650）"
      }
    ]
  }
}
```

成功时返回覆盖率**严格低于** `target_coverage` 的地块，按编号升序排列：

```json
{
  "target_coverage": 0.9,
  "total_parcels_analyzed": 4,
  "below_target_count": 2,
  "parcels": [
    {
      "parcel_id": 102,
      "total_pixels": 2200,
      "covered_pixels": 600,
      "coverage_ratio": 0.2727272727272727,
      "below_target": true
    },
    {
      "parcel_id": 103,
      "total_pixels": 2750,
      "covered_pixels": 500,
      "coverage_ratio": 0.18181818181818182,
      "below_target": true
    }
  ]
}
```

## 样例

仓库自带一组可直接在 Swagger 中上传的样例：

- `samples/parcel_map.tif`：120×100，EPSG:32650，像元 10 m，地块编号 101–104
- `samples/scene_01_cloud_free.tif`：相对地块图原点偏移 `(-7, -5)` 像元
- `samples/scene_02_cloud_free.tif`：相对地块图原点偏移 `(+15, +10)` 像元
- `samples/boundary_update_sample.geojson`：WGS84 边界更正样例，互改 101/102 两块地

两景的幅面尺寸不同、原点错开且空间范围互相重叠，像元格网与地块图对齐，
用于验证窗口映射和多景覆盖合并。

同样可以用 curl 试算：

```bash
curl -X POST http://localhost:8000/analyze \
  -F "parcel=@samples/parcel_map.tif" \
  -F "masks=@samples/scene_01_cloud_free.tif" \
  -F "masks=@samples/scene_02_cloud_free.tif" \
  -F "target_coverage=0.9"
```

样例中地块 101 完全覆盖，地块 104 覆盖率为 2750/3025 ≈ 90.91%，
不低于 0.9；因此只返回地块 102 和 103。

重新生成样例：

```bash
python scripts/generate_samples.py
```

### `POST /update-parcel-boundaries`

县里复测后，用新的地块边界批量更正地块编号图。以 `multipart/form-data` 上传：

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `parcel` | 文件，单个 | 底图地块编号图：北向上的单波段 UInt16/UInt32 GeoTIFF，NoData 为 0 |
| `boundaries` | 文件，单个 | WGS84（EPSG:4326）GeoJSON，FeatureCollection；每个原地块编号对应一个 Polygon，最多 32 个要素 |

GeoJSON 要素的地块编号取 `properties.parcel_id`（也兼容要素顶层 `id`），
必须是图上已存在的正整数。边界会从 WGS84 重投影到底图坐标系，
**像元中心落在多边形内或边上都算命中**（边的判定带 1e-6 像元容差）。

更正规则：

- 本批地块的旧像元视为可复用资源，新边界可以占用同批其他地块的旧址；
  未被任何新边界命中的旧像元置 0。
- 不在本批中的地块像元保持不变。
- 两个新边界命中同一像元（`new_new` 冲突），或新边界占用了未更正地块
  的像元（`new_old` 冲突），都会在错误中给出两个编号和一个冲突像元的
  行列号（从 0 开始）。
- 旧编号在图上不存在，或某条新边界在图内没有命中任何像元，也会分别
  在 `unknown_parcel_ids`、`no_hit_parcel_ids` 中指出编号。
- **有错不出图**：任一校验不通过都返回 HTTP 400，不写出任何像元。
- 整批核对通过后才分块写出新 GeoTIFF，保留原图的尺寸、坐标系、格网、
  数据类型和 NoData，成功时以 `image/tiff` 下载返回。

错误响应示例：

```json
{
  "detail": {
    "message": "边界核对未通过，未生成新的地块编号图。",
    "geojson_errors": [],
    "unknown_parcel_ids": [99],
    "no_hit_parcel_ids": [],
    "conflicts": [
      {
        "type": "new_new",
        "parcel_ids": [101, 102],
        "row": 12,
        "column": 57
      }
    ]
  }
}
```

随附样例 `samples/boundary_update_sample.geojson` 同时更正地块 101、102，
两条新边界各呈 L 形、互占对方一条旧址且互不重叠（生成脚本：
`scripts/generate_boundary_sample.py`）。Swagger 上传示例：

```bash
curl -X POST http://localhost:8000/update-parcel-boundaries \
  -F "parcel=@samples/parcel_map.tif" \
  -F "boundaries=@samples/boundary_update_sample.geojson" \
  --output parcel_map_updated.tif
```

更正后：101 共 2100 个像元、102 共 2300 个像元，各自有 100 个像元落在
对方旧址上；103、104 及行间空隙保持不变。把返回的图再传给 `/analyze`
即可直接用于补采分析。

本地运行与测试（可选）：

```bash
python -m pip install -r requirements-dev.txt
uvicorn app.main:app --reload
pytest -q
```

## 计算与内存

- 以地块图为基准，按其自身分块遍历；异常大的条带/分块会再切到
  ≤1024×1024 像元的窗口（`app.analysis`、`app.boundaries`）。不会整幅载入地块图或掩膜。
- 每个窗口只读取与地块窗口相交的掩膜窗口。越界部分按未覆盖处理。
- 每个地块的总像元数和覆盖像元数按窗口增量累计，多景掩膜用布尔 OR 合并，
  重叠不重复计数。
- 边界更正分两遍走块：第一遍只做命中核对并汇总全部错误，第二遍才向新文件
  写像元；每块内先用栅格化缩小候选像元，再逐候选点做精确的像元中心判定。
- 上传文件先保存到私有临时目录，请求结束后自动删除；不使用上传文件名作为
  本地路径。
