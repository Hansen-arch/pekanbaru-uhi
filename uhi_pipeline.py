"""Urban heat island mapping: Landsat 8/9 thermal band -> LST -> neighbourhood stats.

Inputs (paths set via CLI):
  --lst-tif     Landsat 8/9 ST_B10 GeoTIFF (Collection 2 Level-2, uint16)
  --qa-tif      QA_PIXEL GeoTIFF from the same scene
  --neigh       Neighbourhood polygons (GeoPackage/GeoJSON) with a 'name' column
  --green       Green-space change polygons or raster folder (GeoPackage/GeoJSON,
                with a numeric 'ndvi_change' or 'green_loss' column)
Outputs go to --out (default ./outputs).
"""
import argparse
import os

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask
from rasterio.transform import array_bounds
from scipy import ndimage

# Landsat Collection 2 Level-2 thermal scaling (USGS):
# LST[K] = DN * 0.00341802 + 149.0 ; Celsius = K - 273.15
ST_SCALE = 0.00341802
ST_OFFSET = 149.0
KELVIN_TO_C = 273.15

# QA_PIXEL bits (Collection 2): bit 1 dilated cloud, 2 cirrus, 3 cloud, 4 cloud shadow
CLOUD_BITS = (1 << 1) | (1 << 2) | (1 << 3) | (1 << 4)


def st_to_celsius(dn: np.ndarray) -> np.ndarray:
    """Convert ST_B10 digital numbers to degrees Celsius (float32)."""
    kelvin = dn.astype(np.float32) * ST_SCALE + ST_OFFSET
    return kelvin - KELVIN_TO_C


def cloud_mask(qa: np.ndarray) -> np.ndarray:
    """True where pixel is clear (no dilated cloud, cirrus, cloud or shadow)."""
    return (qa.astype(np.uint16) & CLOUD_BITS) == 0


def load_lst(lst_path: str, qa_path: str):
    with rasterio.open(lst_path) as src, rasterio.open(qa_path) as qsrc:
        dn = src.read(1)
        qa = qsrc.read(1)
        nodata = src.nodata
        lst = st_to_celsius(dn)
        valid = cloud_mask(qa) & (dn != 0)
        if nodata is not None:
            valid &= dn != nodata
        # Outside the scene = zeros connected to the raster border; zeros inside
        # the scene are missing data and are drawn grey like cloud.
        zeros = dn == 0
        labels, _ = ndimage.label(zeros)
        edge = np.unique(np.concatenate([labels[0], labels[-1], labels[:, 0], labels[:, -1]]))
        footprint = ~np.isin(labels, edge[edge != 0])
        lst = np.where(valid, lst, np.nan).astype(np.float32)
        return lst, src.transform, src.crs, src.res, footprint


def zonal_stats(lst, transform, crs, res, polys: gpd.GeoDataFrame,
                unc: np.ndarray = None, src: np.ndarray = None,
                footprint: np.ndarray = None) -> pd.DataFrame:
    polys = polys.to_crs(crs)
    rows = []
    for _, row in polys.iterrows():
        mask = geometry_mask(
            [row.geometry], out_shape=lst.shape, transform=transform,
            invert=True, all_touched=False,
        )
        vals = lst[mask & np.isfinite(lst)]
        rec = {"name": row["name"], "lst_mean": np.nan, "lst_min": np.nan,
               "lst_max": np.nan, "n_pixels": int(vals.size)}
        # Coverage is measured against in-scene area only, so the white swath
        # gaps do not make every cell look under-observed.
        if footprint is not None:
            n_scene = int((mask & footprint).sum())
            rec["pct_clear"] = 100.0 * vals.size / n_scene if n_scene else np.nan
        if unc is not None:
            u = unc[mask & np.isfinite(unc)]
            rec["lst_unc_mean"] = float(u.mean()) if u.size else np.nan
        if src is not None:
            s = src[mask & np.isfinite(lst)]
            rec["pct_regression"] = 100.0 * float((s == 2).mean()) if s.size else np.nan
            rec["pct_composite"] = 100.0 * float((s == 1).mean()) if s.size else np.nan
        if vals.size == 0:
            rows.append(rec)
            continue
        rec.update({"lst_mean": float(vals.mean()), "lst_min": float(vals.min()),
                    "lst_max": float(vals.max())})
        rows.append(rec)
    return pd.DataFrame(rows)


def green_change_by_zone(neigh: gpd.GeoDataFrame, green: gpd.GeoDataFrame) -> pd.Series:
    """Mean green-loss value per neighbourhood. Uses the first numeric column
    named like change/loss, otherwise the first numeric column."""
    green = green.to_crs(neigh.crs)
    num_cols = [c for c in green.select_dtypes("number").columns]
    if not num_cols:
        raise ValueError("green-space layer has no numeric change column")
    pref = [c for c in num_cols if any(k in c.lower() for k in ("loss", "change", "ndvi"))]
    col = (pref or num_cols)[0]
    joined = gpd.sjoin(green[[col, "geometry"]], neigh[["name", "geometry"]],
                       how="inner", predicate="intersects")
    return joined.groupby("name")[col].mean().rename("green_change")


def zone_intervention(df: pd.DataFrame, n: int = 3) -> pd.DataFrame:
    """Rank neighbourhoods: hottest first, weighted by greatest green loss."""
    d = df.copy()
    z_lst = (d["lst_mean"] - d["lst_mean"].mean()) / d["lst_mean"].std(ddof=0)
    if "green_change" in d and d["green_change"].notna().any():
        loss = -d["green_change"].fillna(0.0)  # more negative change = more loss; no polygon = no measured loss
        z_loss = (loss - loss.mean()) / loss.std(ddof=0)
        d["priority"] = 0.6 * z_lst + 0.4 * z_loss
    else:
        d["priority"] = z_lst
    return d.sort_values("priority", ascending=False).head(n)


def heatmap(lst, footprint, transform, neigh: gpd.GeoDataFrame, zones: pd.DataFrame, path: str):
    fig, ax = plt.subplots(figsize=(9, 8))
    extent = array_bounds(lst.shape[0], lst.shape[1], transform)
    left, bottom, right, top = extent[0], extent[1], extent[2], extent[3]
    vmin, vmax = np.nanpercentile(lst, [2, 98])
    # One RGBA image: data colours, grey for cloud/shadow inside the scene,
    # transparent (white) outside the scene footprint.
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    cmap = plt.get_cmap("inferno")
    rgba = cmap(norm(np.where(np.isfinite(lst), lst, vmin)))
    rgba[footprint & ~np.isfinite(lst)] = matplotlib.colors.to_rgba("#6e6e6e")
    rgba[~footprint] = (1, 1, 1, 0)
    ax.imshow(rgba, extent=(left, right, bottom, top))
    im = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    neigh.to_crs(neigh.crs).boundary.plot(ax=ax, color="white", linewidth=0.5)
    order = list(zones["name"])
    for _, r in neigh.iterrows():
        if r["name"] in order:
            c = r.geometry.representative_point()
            # Stack labels by rank so nearby zones do not overprint each other.
            k = order.index(r["name"])
            ax.annotate(str(r["name"]), (c.x, c.y), xytext=(0, 6 + 14 * k),
                        textcoords="offset points", color="cyan", fontsize=9,
                        ha="center", bbox=dict(boxstyle="round,pad=0.2", fc="black",
                                                alpha=0.6, ec="none"))
    fig.colorbar(im, ax=ax, label="Land surface temperature (°C)", extend="both")
    ax.set_title("Surface urban heat island, Pekanbaru\n"
                 "dry-season composite 2022-2024 | cyan = intervention zones, "
                 "grey = cloud not filled, white = outside swath", fontsize=11)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def uncertainty_map(unc, footprint, transform, neigh, path):
    fig, ax = plt.subplots(figsize=(9, 8))
    left, bottom, right, top = array_bounds(unc.shape[0], unc.shape[1], transform)
    rgba = np.zeros(unc.shape + (4,))
    fin = np.isfinite(unc)
    if fin.any():
        cmap = plt.get_cmap("viridis")
        norm = matplotlib.colors.Normalize(0, np.nanpercentile(unc[fin], 99), clip=True)
        rgba[fin] = cmap(norm(unc[fin]))
    else:
        norm, cmap = matplotlib.colors.Normalize(0, 1), plt.get_cmap("viridis")
    rgba[footprint & ~fin] = matplotlib.colors.to_rgba("#e8e8e8")
    rgba[~footprint] = (1, 1, 1, 0)
    ax.imshow(rgba, extent=(left, right, bottom, top))
    neigh.boundary.plot(ax=ax, color="black", linewidth=0.4, alpha=0.5)
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                 label="1σ regression residual (°C)", extend="max")
    ax.set_title("Uncertainty of regression-filled pixels\n"
                 "light grey = observed (no fill needed), white = outside scene")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def provenance_map(src, footprint, transform, neigh, path):
    fig, ax = plt.subplots(figsize=(9, 8))
    left, bottom, right, top = array_bounds(src.shape[0], src.shape[1], transform)
    colours = {0: "#d73027", 1: "#fdae61", 2: "#4575b4"}
    rgba = np.zeros(src.shape + (4,))
    for k, c in colours.items():
        rgba[src == k] = matplotlib.colors.to_rgba(c)
    rgba[~footprint] = (1, 1, 1, 0)
    ax.imshow(rgba, extent=(left, right, bottom, top))
    neigh.boundary.plot(ax=ax, color="black", linewidth=0.4, alpha=0.6)
    handles = [matplotlib.patches.Patch(color=colours[0], label="primary scene (13 Jun 2023)"),
               matplotlib.patches.Patch(color=colours[1], label="medoid composite (other dates)"),
               matplotlib.patches.Patch(color=colours[2], label="regression fill (modelled)")]
    ax.legend(handles=handles, loc="lower left", fontsize=8, framealpha=0.9)
    ax.set_title("Source of each LST pixel")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lst-tif", required=True)
    ap.add_argument("--qa-tif", required=True)
    ap.add_argument("--neigh", required=True)
    ap.add_argument("--green", required=False)
    ap.add_argument("--unc-tif", required=False)
    ap.add_argument("--source-tif", required=False)
    ap.add_argument("--footprint-tif", required=False)
    ap.add_argument("--out", default="outputs")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    lst, transform, crs, res, footprint = load_lst(a.lst_tif, a.qa_tif)
    if a.footprint_tif and os.path.exists(a.footprint_tif):
        # Explicit swath mask from the composite step, so interior data gaps
        # are not mistaken for the edge of the scene.
        with rasterio.open(a.footprint_tif) as f:
            footprint = f.read(1) > 0
    neigh = gpd.read_file(a.neigh)
    if "name" not in neigh.columns:
        raise ValueError("neighbourhood layer needs a 'name' column")

    unc = None
    if a.unc_tif and os.path.exists(a.unc_tif):
        with rasterio.open(a.unc_tif) as u:
            unc = u.read(1)
    src = None
    if a.source_tif and os.path.exists(a.source_tif):
        with rasterio.open(a.source_tif) as s:
            src = s.read(1)
    stats = zonal_stats(lst, transform, crs, res, neigh, unc=unc, src=src,
                        footprint=footprint)
    if a.green:
        green = gpd.read_file(a.green)
        gc = green_change_by_zone(neigh.to_crs(crs), green)
        stats = stats.merge(gc, left_on="name", right_index=True, how="left")
    else:
        stats["green_change"] = np.nan

    stats = stats.sort_values("lst_mean", ascending=False)
    if stats["green_change"].notna().any():
        print("Correlation LST mean vs green change:",
              round(stats["lst_mean"].corr(stats["green_change"]), 3))
    stats.to_csv(os.path.join(a.out, "zonal_stats.csv"), index=False)

    zones = zone_intervention(stats)
    zones.to_csv(os.path.join(a.out, "intervention_zones.csv"), index=False)
    heatmap(lst, footprint, transform, neigh.to_crs(crs), zones,
            os.path.join(a.out, "uhi_heatmap.png"))
    if unc is not None:
        uncertainty_map(unc, footprint, transform, neigh.to_crs(crs),
                        os.path.join(a.out, "lst_uncertainty_map.png"))
    if src is not None:
        provenance_map(src, footprint, transform, neigh.to_crs(crs),
                       os.path.join(a.out, "lst_provenance_map.png"))
    print(stats.to_string(index=False))
    print("\nTop 3 intervention zones:")
    print(zones[["name", "lst_mean", "green_change", "priority"]].to_string(index=False))


if __name__ == "__main__":
    main()
