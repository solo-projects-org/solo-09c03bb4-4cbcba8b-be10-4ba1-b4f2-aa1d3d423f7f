# 农情遥感补采分析 API

用于在安排补采前，从地块编号图和各景无云掩膜中找出覆盖率不足的地块；
也支持把县里复测得到的地块边界（WGS84 GeoJSON）批量写回地块编号图。
服务提供 FastAPI + Swagger UI，可直接上传 GeoTIFF / GeoJSON 试算。

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

### `POST /correct-boundaries`

以 `multipart/form-data` 上传：

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `parcel` | 文件，单个 | 原地块编号图：北向上、单波段 UInt16/UInt32 GeoTIFF，NoData 为 0 |
| `boundaries` | 文件，单个 | WGS84（EPSG:4326）GeoJSON `FeatureCollection`，每个原地块编号对应一个 `Polygon` |

GeoJSON 中每个 Feature 的编号取自 `properties.parcel_id`（也兼容数字类型的
Feature `id`），一批最多 32 个，编号不可重复。成功时返回更正后的 GeoTIFF
（`image/tiff`），保留原图的尺寸、坐标系、像元格网和数据类型。

更正规则：

- 多边形先重投影到地块图坐标系（长边自动加密，避免投影切角），再按**像元中心**
  判定归属；像元中心落在多边形内或边上都算命中（`covers` 语义）。
- 新边界命中的像元写入对应编号；本批更正地块的旧像元若未被任何新边界重新命中，
  置为 NoData `0`；不在本批中的地块原样保留。
- 整批核对全部通过后才出图。任一检查失败返回 HTTP 400，`detail.errors` 列出
  全部编号级问题，且不返回影像：

| `kind` | 含义 |
| --- | --- |
| `new_new_conflict` | 两个新边界命中同一像元，给出两个编号和一个冲突像元的 `row`/`column` |
| `new_old_conflict` | 新边界占用了未更正地块的像元，给出新旧两个编号和冲突像元行列号 |
| `unknown_parcel_id` | 编号在原地块图中不存在 |
| `no_pixels_hit` | 新边界在地块图范围内一个像元都没命中 |
| `invalid_raster` / `invalid_geojson` / `invalid_geometry` | 输入结构不合法 |
| `duplicate_parcel_id` / `too_many_parcels` / `id_out_of_range` | 批次结构问题 |

冲突响应示例：

```json
{
  "detail": {
    "message": "边界更正核对未通过，未输出影像。",
    "errors": [
      {
        "kind": "new_new_conflict",
        "message": "地块 101 与 102 的新边界同时命中像元（行 0，列 40）。",
        "parcel_ids": [101, 102],
        "column": 40,
        "row": 0
      }
    ]
  }
}
```

curl 调用：

```bash
curl -X POST http://localhost:8000/correct-boundaries \
  -F "parcel=@samples/parcel_map.tif" \
  -F "boundaries=@samples/boundary_corrections_101_102.geojson" \
  --output parcel_corrected.tif
```

## 样例

仓库自带一组可直接在 Swagger 中上传的样例：

- `samples/parcel_map.tif`：120×100，EPSG:32650，像元 10 m，地块编号 101–104
- `samples/scene_01_cloud_free.tif`：相对地块图原点偏移 `(-7, -5)` 像元
- `samples/scene_02_cloud_free.tif`：相对地块图原点偏移 `(+15, +10)` 像元

两景的幅面尺寸不同、原点错开且空间范围互相重叠，像元格网与地块图对齐，
用于验证窗口映射和多景覆盖合并。

- `samples/boundary_corrections_101_102.geojson`：WGS84 边界更正样例，
  重定义地块 101、102，使两块地各占用对方的部分旧址（101 向下扩到第 30 行、
  102 向上扩到第 30 行）。可在 Swagger 的 `/correct-boundaries` 中与
  `parcel_map.tif` 一起上传；更正后 103、104 不变，101/102 未被重新命中的
  旧像元置 0。

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

本地运行与测试（可选）：

```bash
python -m pip install -r requirements-dev.txt
uvicorn app.main:app --reload
pytest -q
```

## 计算与内存

- 以地块图为基准，按其自身分块遍历；异常大的条带/分块会再切到
  ≤1024×1024 像元的窗口（`app.analysis`、`app.boundaries`）。不会整幅载入
  地块图或掩膜。
- `/analyze` 每个窗口只读取与地块窗口相交的掩膜窗口，越界部分按未覆盖处理；
  每个地块的像元数按窗口增量累计，多景掩膜用布尔 OR 合并，重叠不重复计数。
- `/correct-boundaries` 在每个窗口内只对落在多边形包围盒内的候选像元做点-多边形
  判断（shapely 向量化谓词），冲突、越界与编号核对全部在分块扫描中完成；
  只有全部检查通过，写出的 GeoTIFF 才会作为响应返回，否则临时文件直接丢弃。
- 上传文件先保存到私有临时目录，请求结束后自动删除；不使用上传文件名作为
  本地路径。
