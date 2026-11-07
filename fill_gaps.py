"""Build a per-pixel clear-sky thermal composite and fill gaps in inputs/lst.tif.

For each pixel, the composite uses the medoid of the clear observations across
several 2022-2024 Landsat 8/9 scenes on path 127 row 59 (medoid: one real
reading per pixel, the clear observation nearest the per-pixel median). The
primary scene (LC09_127059_20230613) is kept where it is clear. Filled values
come from different dates, so they are not strictly comparable.

Writes:
  inputs/lst_primary.tif, inputs/qa_primary.tif   backup of the original inputs
  inputs/lst.tif, inputs/qa.tif                   composite used by uhi_pipeline.py
"""
import os
import shutil

import ee
import numpy as np
import rasterio
from affine import Affine

from fetch_inputs import AOI_BOUNDS, OUT, download_tif

PRIMARY = "LANDSAT/LC09/C02/T1_L2/LC09_127059_20230613"
N_SCENES = 15


def clear_mask(qa, dn):
    cloud_bits = (1 << 1) | (1 << 2) | (1 << 3) | (1 << 4)
    return ((qa & cloud_bits) == 0) & (dn != 0)


def box_mean(a, w):
    pad = w // 2
    p = np.pad(a, pad, mode="edge")
    c = np.pad(np.cumsum(np.cumsum(p, 0), 1), ((1, 0), (1, 0)))
    s = c[w:, w:] - c[:-w, w:] - c[w:, :-w] + c[:-w, :-w]
    return s / (w * w)


def regress_fill(dn, valid, win=9):
    """Fill invalid pixels by a linear regression on the local neighbourhood mean.

    The neighbourhood mean excludes the pixel itself, so the fit does not
    use the value it predicts. Returns the filled DN, the 1-sigma residual
    in DN for each filled pixel, and the mask of filled pixels.
    """
    v = np.where(valid, dn, 0).astype(np.float64)
    ok = valid.astype(np.float64)
    vsum = box_mean(v, win) * win * win - v
    osum = box_mean(ok, win) * win * win - ok
    nbr = np.where(osum > 0, vsum / np.maximum(osum, 1e-9), np.nan)
    fit = valid & np.isfinite(nbr)
    b, a = np.polyfit(nbr[fit], v[fit], 1)
    # Local 1-sigma: residuals pooled in a moving window, so heterogeneous
    # (urban) areas carry a wider error bar than uniform ones.
    resid2 = np.where(fit, (v - (a + b * np.nan_to_num(nbr))) ** 2, 0.0)
    num = box_mean(resid2, 25)
    den = box_mean(fit.astype(np.float64), 25)
    sigma_local = np.sqrt(num / np.maximum(den, 1e-6))
    gsig = float(np.std(v[fit] - (a + b * nbr[fit])))
    sigma_local = np.where(np.isfinite(sigma_local) & (sigma_local > 0.1), sigma_local, gsig)
    gap = ~valid & np.isfinite(nbr)
    filled = dn.copy()
    filled[gap] = np.round(a + b * nbr[gap]).astype(dn.dtype)
    sig = np.full(dn.shape, np.nan, dtype=np.float32)
    sig[gap] = np.maximum(sigma_local[gap] * 0.00341802, 0.3).astype(np.float32)  # DN to kelvin
    return filled, sig, gap


def candidate_ids():
    aoi = ee.Geometry.Rectangle(AOI_BOUNDS)
    coll = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2").merge(
        ee.ImageCollection("LANDSAT/LC08/C02/T1_L2"))
    coll = coll.filterBounds(aoi).filterDate("2022-01-01", "2024-12-31") \
        .filter(ee.Filter.eq("WRS_PATH", 127)).filter(ee.Filter.eq("WRS_ROW", 59)) \
        .filter(ee.Filter.lt("CLOUD_COVER", 60))
    return coll.aggregate_array("system:index").getInfo()


def fetch(scene_id, tries=3):
    """Download ST_B10 and QA_PIXEL, caching to disk so retries are cheap."""
    cache = os.path.join(OUT, "scene_cache")
    os.makedirs(cache, exist_ok=True)
    path = os.path.join(cache, f"{scene_id}.npz")
    if os.path.exists(path):
        z = np.load(path)
        return z["dn"], z["qa"]
    img = ee.Image(_asset(scene_id))
    for k in range(tries):
        try:
            arr, _, _ = download_tif(img.select(["ST_B10", "QA_PIXEL"]),
                                     ee.Geometry.Rectangle(AOI_BOUNDS), None)
            dn, qa = arr[0], arr[1].astype(np.uint16)
            np.savez_compressed(path, dn=dn, qa=qa)
            return dn, qa
        except Exception as e:
            if k == tries - 1:
                raise
            print(f"retry {scene_id} after {type(e).__name__}")
            import time
            time.sleep(5 * (k + 1))


def _asset(index):
    # index looks like "2_LC08_127059_20220704"; the collection is the part after the prefix
    _, prod, _, _ = index.split("_")
    coll = "LANDSAT/LC09/C02/T1_L2" if prod == "LC09" else "LANDSAT/LC08/C02/T1_L2"
    return f"{coll}/{index.split('_', 1)[1]}"


def write_uncertainty(sig, profile):
    p = profile.copy()
    p.update(dtype="float32", nodata=np.nan)
    with rasterio.open(os.path.join(OUT, "lst_uncertainty.tif"), "w", **p) as d:
        d.write(sig, 1)


def write_source(src, profile):
    """0 = primary scene, 1 = medoid composite, 2 = regression fill."""
    p = profile.copy()
    p.update(dtype="uint8", nodata=None)
    with rasterio.open(os.path.join(OUT, "lst_source.tif"), "w", **p) as d:
        d.write(src.astype(np.uint8), 1)


def main():
    ee.Initialize()
    primary_lst = os.path.join(OUT, "lst.tif")
    primary_qa = os.path.join(OUT, "qa.tif")
    backup_lst = os.path.join(OUT, "lst_primary.tif")
    backup_qa = os.path.join(OUT, "qa_primary.tif")
    if not os.path.exists(backup_lst):
        shutil.copy(primary_lst, backup_lst)
        shutil.copy(primary_qa, backup_qa)

    with rasterio.open(backup_lst) as s, rasterio.open(backup_qa) as q:
        dn0 = s.read(1)
        qa0 = q.read(1).astype(np.uint16)
        profile = s.profile.copy()
        qa_profile = q.profile.copy()

    ids = candidate_ids()
    print("candidates", len(ids))

    # score each candidate by clear share, keep the best N
    scored = []
    for sid in ids:
        dn, qa = fetch(sid)
        scored.append((clear_mask(qa, dn).mean(), sid, dn, qa))
    scored.sort(key=lambda t: -t[0])
    chosen = scored[:N_SCENES]
    print("using", len(chosen), "scenes, best clear share", round(100 * chosen[0][0], 1))

    stack_dn = np.stack([np.where(clear_mask(qa, dn), dn, 0) for _, _, dn, qa in chosen]).astype(np.float64)
    stack_ok = np.stack([clear_mask(qa, dn) for _, _, dn, qa in chosen])
    counts = stack_ok.sum(axis=0)
    # medoid: per pixel, the clear observation closest to the per-pixel median.
    # Unlike a median, it is always one real reading from one date.
    med = np.nanmedian(np.where(stack_ok, stack_dn, np.nan), axis=0)
    med = np.nan_to_num(med, nan=0.0)
    dist = np.where(stack_ok, np.abs(stack_dn - med[None]), np.inf)
    pick = np.argmin(dist, axis=0)
    comp_dn = np.take_along_axis(stack_dn, pick[None], axis=0)[0].astype(dn0.dtype)
    comp_ok = counts > 0

    primary_ok = clear_mask(qa0, dn0)
    use_comp = ~primary_ok & comp_ok
    out_dn = np.where(use_comp, comp_dn, dn0).astype(dn0.dtype)
    out_qa = np.where(use_comp, 0, qa0).astype(qa0.dtype)  # 0 = clear in the output

    # regression fill for pixels still masked after the medoid step.
    # Restricted to the scene footprint so nothing is invented outside the swath.
    # The swath is the convex hull of every pixel any scene ever observed.
    # Zeros beyond it are the real scene edge, so nothing is invented outside
    # the satellite swath, while never-observed cloud inside it can be filled.
    from rasterio.features import rasterize
    from shapely.geometry import MultiPoint
    obs = (dn0 != 0) | (counts > 0)
    ys, xs = np.nonzero(obs)
    hull = MultiPoint(np.column_stack([xs, ys])).convex_hull
    inside = rasterize([(hull, 1)], out_shape=dn0.shape,
                       transform=Affine.identity(), fill=0).astype(bool)
    # Repeated passes: each one fills a ring of pixels further into the cloud
    # blobs, because newly filled pixels become neighbours for the next pass.
    valid = inside & clear_mask(out_qa, out_dn)
    sig = np.full(out_dn.shape, np.nan, dtype=np.float32)
    reg_gap = np.zeros(out_dn.shape, dtype=bool)
    work = out_dn.copy()
    for it in range(4):
        filled, sig_it, gap = regress_fill(work, valid)
        gap &= inside & ~valid
        if not gap.any():
            break
        work[gap] = filled[gap]
        sig[gap] = sig_it[gap]
        valid |= gap
        reg_gap |= gap
    out_dn = np.where(reg_gap, work, out_dn).astype(dn0.dtype)
    out_qa = np.where(reg_gap, 0, out_qa).astype(qa0.dtype)
    print(f"regression passes stopped after {it + 1}")
    write_uncertainty(sig, profile)
    source = np.zeros(dn0.shape, dtype=np.uint8)
    source[use_comp] = 1
    source[reg_gap] = 2
    write_source(source, profile)
    p = profile.copy()
    p.update(dtype="uint8", nodata=None)
    with rasterio.open(os.path.join(OUT, "lst_footprint.tif"), "w", **p) as d:
        d.write(inside.astype(np.uint8), 1)

    with rasterio.open(primary_lst, "w", **profile) as d:
        d.write(out_dn, 1)
    with rasterio.open(primary_qa, "w", **qa_profile) as d:
        d.write(out_qa, 1)

    print(f"primary clear   {100 * primary_ok.mean():.1f}%")
    print(f"filled          {100 * use_comp.sum() / dn0.size:.1f}% of raster")
    print(f"regression fill {100 * reg_gap.sum() / dn0.size:.1f}% of raster")
    if reg_gap.any():
        print(f"regression 1-sigma  median {np.median(sig[reg_gap]):.2f} K, "
              f"p95 {np.percentile(sig[reg_gap], 95):.2f} K")
    print(f"composite clear {100 * clear_mask(out_qa, out_dn).mean():.1f}%")


if __name__ == "__main__":
    main()
