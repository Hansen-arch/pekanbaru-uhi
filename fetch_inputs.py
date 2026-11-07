"""Build pipeline inputs for the Pekanbaru AOI from Google Earth Engine.

Writes to ./inputs:
  lst.tif, qa.tif      Landsat 9 C02 L2 ST_B10 and QA_PIXEL, EPSG:32647, 30 m
  neigh.gpkg           2 km grid cells used as neighbourhood stand-ins (name = cell id)
  green.gpkg           per-cell NDVI change (2022-2024 minus 2016-2018, dry season)
"""
import os

import ee
import geopandas as gpd
import numpy as np
import rasterio
import requests
from rasterio.io import MemoryFile
from shapely.geometry import box

CRS = "EPSG:32647"
OUT = "inputs"
AOI_BOUNDS = [101.35, 0.40, 101.55, 0.60]
SCENE = "LANDSAT/LC09/C02/T1_L2/LC09_127059_20230613"
CELL_M = 2000


def download_tif(img, region, path):
    url = img.getDownloadURL({"scale": 30, "crs": CRS, "region": region,
                              "format": "GEO_TIFF", "filePerBand": False})
    r = requests.get(url, timeout=600)
    r.raise_for_status()
    with MemoryFile(r.content) as mf, mf.open() as src:
        profile = src.profile.copy()
        arr = src.read()
        names = list(src.descriptions)
    return arr, profile, names


def fetch_scene():
    img = ee.Image(SCENE)
    aoi = ee.Geometry.Rectangle(AOI_BOUNDS)
    arr, prof, names = download_tif(img.select(["ST_B10", "QA_PIXEL"]), aoi, None)
    st_idx = names.index("ST_B10") if "ST_B10" in names else 0
    qa_idx = names.index("QA_PIXEL") if "QA_PIXEL" in names else 1
    prof.update(count=1, compress="deflate")
    for fname, idx in [("lst.tif", st_idx), ("qa.tif", qa_idx)]:
        with rasterio.open(os.path.join(OUT, fname), "w", **prof) as d:
            d.write(arr[idx], 1)
    return prof["crs"], prof["transform"]


def grid_geojson(crs_str):
    """2 km grid over the AOI, built in UTM so cells are square."""
    ll = gpd.GeoDataFrame(geometry=[box(101.35, 0.40, 101.55, 0.60)], crs="EPSG:4326")
    utm = ll.to_crs(CRS).total_bounds
    xs = np.arange(utm[0], utm[2], CELL_M)
    ys = np.arange(utm[1], utm[3], CELL_M)
    cells, names = [], []
    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            cells.append(box(x, y, x + CELL_M, y + CELL_M))
            names.append(f"C{i:02d}{j:02d}")
    g = gpd.GeoDataFrame({"name": names}, geometry=cells, crs=CRS)
    return g


def ndvi_change_by_cell(grid: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    def dry_composite(start, end):
        cols = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2").merge(
            ee.ImageCollection("LANDSAT/LC09/C02/T1_L2"))
        cols = cols.filterBounds(ee.Geometry.Rectangle(AOI_BOUNDS)).filterDate(start, end) \
            .filter(ee.Filter.lt("CLOUD_COVER", 30))

        def mask_and_ndvi(i):
            qa = i.select("QA_PIXEL")
            clear = qa.bitwiseAnd(1 << 3).eq(0).And(qa.bitwiseAnd(1 << 4).eq(0))
            sr = i.select(["SR_B5", "SR_B4"]).multiply(0.0000275).add(-0.2)
            ndvi = sr.normalizedDifference(["SR_B5", "SR_B4"]).rename("ndvi")
            return ndvi.updateMask(clear)
        return ee.ImageCollection(cols.map(mask_and_ndvi)).median()

    before = dry_composite("2016-01-01", "2018-12-31")
    after = dry_composite("2022-01-01", "2024-12-31")
    change = after.subtract(before).rename("ndvi_change")

    grid_ll = grid.to_crs("EPSG:4326")
    fc = ee.FeatureCollection([
        ee.Feature(ee.Geometry(row.geometry.__geo_interface__, "EPSG:4326"), {"name": row["name"]})
        for _, row in grid_ll.iterrows()
    ])
    stats = change.reduceRegions(collection=fc, reducer=ee.Reducer.mean(),
                                 scale=30, tileScale=2).getInfo()
    recs = {f["properties"]["name"]: f["properties"].get("mean") for f in stats["features"]}
    out = grid.copy()
    out["ndvi_change"] = out["name"].map(recs)
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    ee.Initialize()
    crs, transform = fetch_scene()
    print("scene written", crs, transform)
    grid = grid_geojson(CRS)
    grid.to_file(os.path.join(OUT, "neigh.gpkg"), driver="GPKG")
    green = ndvi_change_by_cell(grid)
    green[["ndvi_change", "geometry"]].to_file(os.path.join(OUT, "green.gpkg"), driver="GPKG")
    print("cells", len(grid), "with change", int(green["ndvi_change"].notna().sum()))
    print("ndvi_change mean", round(float(green["ndvi_change"].mean()), 3))


if __name__ == "__main__":
    main()
