"""Check saved evidence and regenerate the small tables without numerical dependencies."""

import argparse
import csv
import io
import json
import math
import sys

from evidence import (ARCHIVE_SHA256, ARCHIVE_URL, CONFIGURATIONS, DATA, EVIDENCE,
                      MAIN, REVISION, ROOT, SOLVER_URL, one, records, relaxation,
                      resolution, sha256)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def csv_bytes(rows):
    """Use round-trip float representations and a stable column order."""
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return out.getvalue().encode()


def json_bytes(value):
    return (json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode()


def derive():
    tables, claims, sources = {}, {}, {}

    def source(prov):
        key = f"{prov['file']}:{prov['line']}"
        sources[key] = prov
        return key

    def ratio(name, numerator, denominator, nsource, dsource, nfield, dfield, reference):
        value = numerator / denominator
        claims[name] = {"value": value, "numerator": numerator, "denominator": denominator,
                        "numerator_source": source(nsource), "denominator_source": source(dsource),
                        "numerator_field": nfield, "denominator_field": dfield,
                        "reference": reference}
        return value

    require(len(records("tc6_week2_relax.jsonl")) == 12, "Expected 12 archived configurations")
    baseline, pb = relaxation("none", 0, 0)
    rr = []
    selected = {}
    for region, par, rotation in CONFIGURATIONS:
        r, p = relaxation(region, par, rotation)
        selected[region, par, rotation] = r, p
        require(set(r["l2_vs_T85_by_day"]) == {"1", "3", "5", "7", "10"},
                "The archived relaxation time grid has changed")
        for day, error in sorted(r["l2_vs_T85_by_day"].items(), key=lambda x: int(x[0])):
            rr.append(dict(region=region, par_deg=par, control_rotation_deg=rotation, M=128,
                           day=int(day), reference="T85dt45", weighted_area=r["weighted_area"],
                           l2=error, baseline_l2=baseline["l2_vs_T85_by_day"][day],
                           reduction_vs_none=baseline["l2_vs_T85_by_day"][day] / error,
                           source=source(p)))
    tables["relaxation.csv"] = csv_bytes(rr)
    for key, label in [(("bp", 20, 0), "cones_day10_reduction"),
                       (("bp", 20, 60), "controls_day10_reduction"),
                       (("polar", 60, 0), "polar_day10_reduction")]:
        r, p = selected[key]
        ratio(label, baseline["l2_vs_T85_by_day"]["10"], r["l2_vs_T85_by_day"]["10"],
              pb, p, 'l2_vs_T85_by_day.10', 'l2_vs_T85_by_day.10', "T85dt45")
    pairs = [("bp", 20), ("disc", 20), ("annulus", 10), ("disc", 10)]
    claims["paired_shapes_favouring_cones"] = {
        "value": sum(selected[r, p, 0][0]["l2_vs_T85_by_day"]["10"] <
                     selected[r, p, 60][0]["l2_vs_T85_by_day"]["10"] for r, p in pairs),
        "total": 4, "pairs": [list(x) for x in pairs], "reference": "T85dt45"}
    for rotation, name in [(0, "cones"), (60, "controls")]:
        r, p = selected["bp", 20, rotation]
        claims[f"{name}_weighted_area"] = {"value": r["weighted_area"], "source": source(p),
                                           "field": "weighted_area"}

    res, context = [], []

    def series(r, p, method, mesh, size, unique, computed, group, reference="T213dt22"):
        for day, error in sorted(r["l2_by_day"].items(), key=lambda x: int(x[0])):
            context.append(dict(comparison=group, method=method, grid=mesh, grid_size=size,
                                unique_points=unique, computed_height_sites=computed,
                                day=int(day), reference=reference, l2=error, source=source(p)))

    for M, trunc, dt in [(128, 42, 900.0), (256, 85, 450.0), (512, 170, 225.0)]:
        endpoints = resolution(M, trunc, dt)
        spectral, ps = endpoints[2]
        for method, (r, p) in zip(("main", "trsk", "spectral"), endpoints):
            unique = M * M // 2
            computed = r["ncells"] if method == "trsk" else r["npts"]
            value = ratio(f"{method}_over_spectral_M{M}", r["l2"], spectral["l2"], p, ps,
                          "l2", "l2", "T213dt22")
            res.append(dict(method=method, M=M if method != "spectral" else "",
                            truncation=trunc if method == "spectral" else "", day=14,
                            reference="T213dt22", unique_points=unique,
                            computed_height_sites=computed, record_npts=r["npts"],
                            velocity_edges=r.get("nedges", ""), dt_s=r["dt"],
                            l2=r["l2"], spectral_l2=spectral["l2"], error_over_spectral=value,
                            source=source(p)))
            if method == "main":
                r, p = one("tc6_series.jsonl", {"M": M, "scheme": MAIN})
            series(r, p, method, "gaussian" if method == "spectral" else "torus",
                   trunc if method == "spectral" else M, unique, computed, "equal_unique_points")
    tables["resolution.csv"] = csv_bytes(res)
    fine = claims["trsk_over_spectral_M512"]["value"]
    claims["fine_trsk_percent_above_spectral"] = {"value": 100 * (fine - 1),
        "formula": "100 * (trsk_over_spectral_M512 - 1)", "reference": "T213dt22"}

    early = []
    for M, level in [(144, 5), (288, 6), (576, 7)]:
        main, pm = one("tc6_series.jsonl", {"M": M, "scheme": MAIN})
        trsk, pt = one("mpas_ref.jsonl", {"solver": "trsk", "case": "tc6",
                                           "grid": ["icos", level], "days": 14.0})
        series(main, pm, "main", "torus", M, main["nunique"], M * M, "near_equal_icos_cells")
        series(trsk, pt, "trsk", "icos", level, trsk["nunique"], trsk["ncells"], "near_equal_icos_cells")
        for day in (1, 3, 5):
            early.append(ratio(f"icos_over_main_M{M}_day{day}", trsk["l2_by_day"][str(day)],
                               main["l2_by_day"][str(day)], pt, pm, f"l2_by_day.{day}",
                               f"l2_by_day.{day}", "T213dt22"))
    claims["early_icos_over_main_range"] = {"minimum": min(early), "maximum": max(early),
                                          "days": [1, 3, 5], "M": [144, 288, 576]}
    tables["context.csv"] = csv_bytes(context)

    restart = []
    for _, r in records("tc6_week2_restart.jsonl"):
        _, p = one("tc6_week2_restart.jsonl", {k: r[k] for k in ("M", "label", "start_day", "scheme")})
        for day, error in sorted(r["l2_by_day"].items(), key=lambda x: int(x[0])):
            extra = r["after_3_days"] if int(day) == r["after_3_days"]["day"] else {}
            row = dict(M=r["M"], label=r["label"], scheme=r["scheme"], start_day=r["start_day"],
                       day=int(day), elapsed_day=int(day) - r["start_day"], reference="T213dt22",
                       l2=error, zonal_mean_share=extra.get("zonal_mean_share", ""))
            for lo in range(-90, 90, 30):
                row[f"rms_lat_{lo}_to_{lo+30}_m"] = extra.get("rms_by_lat_m", {}).get(f"{lo}..{lo+30}", "")
            row["source"] = source(p)
            restart.append(row)
    tables["restart.csv"] = csv_bytes(restart)
    start, p0 = one("tc6_week2_restart.jsonl", {"M": 128, "label": "main", "start_day": 0, "scheme": MAIN})
    late, p7 = one("tc6_week2_restart.jsonl", {"M": 128, "label": "main", "start_day": 7, "scheme": MAIN})
    ratio("restart_elapsed3_ratio", late["l2_by_day"]["10"], start["l2_by_day"]["3"],
          p7, p0, "l2_by_day.10", "l2_by_day.3", "T213dt22")
    twins = []
    for shape in ("random", "zonal_P2", "tetrahedral_l3"):
        r, p = one("tc6_week2_twins7.jsonl", {"M": 128, "shape": shape, "scheme": MAIN})
        for day, value in [(7, r["initial_rms_m"])] + sorted(
                [(int(d), v) for d, v in r["rms_difference_m_by_day"].items()]):
            twins.append(dict(M=128, shape=shape, start_day=7, day=day, elapsed_day=day-7,
                              initial_condition=day == 7, rms_difference_m=value,
                              reference="paired_solver_runs", source=source(p)))
        claims[f"twin_{shape}_day14_rms_m"] = {"value": r["rms_difference_m_by_day"]["14"],
            "initial_rms_m": r["initial_rms_m"], "source": source(p), "field": "rms_difference_m_by_day.14"}
    tables["twins.csv"] = csv_bytes(twins)
    for name, selector in [
        ("target_T85dt45_day14_l2", {"run": "T85_tc6_2h", "case": "tc6", "against": "T213dt22"}),
        ("reference_truncation_difference", {"run": "T213dt90_tc6", "against": "T341dt90"}),
        ("reference_timestep_difference", {"run": "T213dt45_tc6", "against": "T213dt22"}),
    ]:
        r, p = one("stswm_errors.jsonl", selector)
        claims[name] = {"value": r["l2_by_day"]["14"], "source": source(p),
                        "field": "l2_by_day.14", "reference": r["against"]}
    for day in (0, 3, 7, 10, 14):
        r, p = one("tc6_week2_sharpness.jsonl", {"mode": "sharpness", "day": day})
        claims[f"high_latitude_length_day{day}_km"] = {
            "value": r["bands"]["|lat| 60-80"]["length_m"] / 1000,
            "source": source(p), "field": "bands.|lat| 60-80.length_m", "divisor": 1000}
    expected = {"cones_day10_reduction": 5.33265185658, "controls_day10_reduction": 2.81015510331,
                "restart_elapsed3_ratio": 12.7711668561, "trsk_over_spectral_M512": 1.04582187688,
                "target_T85dt45_day14_l2": 1.3166923113e-3}
    for name, value in expected.items():
        require(math.isclose(claims[name]["value"], value, rel_tol=5e-11, abs_tol=0),
                f"Headline does not reproduce: {name}")
    require(claims["paired_shapes_favouring_cones"]["value"] == 4, "A paired comparison has changed")
    require([f'{claims[f"main_over_spectral_M{m}"]["value"]:.2f}' for m in (128, 256, 512)] ==
            ["1.86", "1.94", "1.92"], "Main/spectral ratios do not reproduce")
    require([f'{claims[f"trsk_over_spectral_M{m}"]["value"]:.2f}' for m in (128, 256)] ==
            ["1.91", "1.54"], "TRSK/spectral ratios do not reproduce")
    tables["claims.json"] = json_bytes({"schema": 1, "experiment": "archived",
        "metric": "sqrt(sum(w * (h - h_ref)^2) / sum(w * h_ref^2))",
        "claims": claims, "sources": sources})
    return tables, sources, claims


def check_fields():
    """Reconcile the bundled field measurements with their summary and run sidecars."""
    summary = json.loads((DATA / "visual-runs.json").read_text())
    validation = json.loads((DATA / "validation.json").read_text())
    with (DATA / "field-metrics.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    require(len(rows) == 3 * 2 * 169 * 3, "Field table must contain 3,042 measurements")
    indexed = {}
    for row in rows:
        key = row["run"], row["reference"], float(row["t_hours"]), row["region"]
        require(key not in indexed, f"Duplicate field measurement: {key}")
        indexed[key] = row
        numerator = float(row["weighted_error_sq_m2"])
        denominator = float(row["reference_height_sq_m2"])
        require(math.isclose(float(row["l2"]), math.sqrt(numerator / denominator), rel_tol=1e-9),
                f"Field norm does not reconcile: {key}")
        require(float(row["pair_rms_difference_m"]) == float(row["pair_max_difference_m"]) == 0,
                f"Unequal double-cover errors: {key}")
    for name, run in summary["runs"].items():
        sidecar = json.loads((ROOT / run["sidecar"]).read_text())
        require(run["l2_vs_T85dt45_by_day"] == sidecar["l2_vs_T85_by_day"],
                f"T85 sidecar scores differ: {name}")
        times = sorted(k[2] for k in indexed if k[0] == name and k[1] == "T213dt22" and k[3] == "global")
        require(len(times) == 169 and times[0] == 0 and times[-1] == 336,
                f"Unexpected snapshot coverage: {name}")
        for hour in times:
            for reference in ("T85dt45", "T213dt22"):
                for region in ("global", "cones20", "controls20"):
                    require((name, reference, hour, region) in indexed,
                            f"Missing field measurement: {name}, {reference}, {hour}, {region}")
        for day, value in run["l2_vs_T213dt22_by_day"].items():
            row = indexed[name, "T213dt22", 24 * int(day), "global"]
            require(value == float(row["l2"]), f"T213 summary differs: {name}, day {day}")
        require(validation["runs"][name]["n_snapshots"] == len(times),
                f"Snapshot validation differs: {name}")
    for reference, treatments in summary["reductions_vs_main"].items():
        key = f"l2_vs_{reference}_by_day"
        for name, days in treatments.items():
            for day, ratio in days.items():
                expected = summary["runs"]["main"][key][day] / summary["runs"][name][key][day]
                require(ratio == expected, f"Visual reduction differs: {reference}, {name}, day {day}")
    comparison = summary["main_vs_archived_T213"]
    archived, _ = one("tc6_series.jsonl", comparison["selector"])
    for day, item in comparison["days"].items():
        require(item["archived_l2"] == archived["l2_by_day"][day], "Archived T213 score differs")
        require(item["field_l2"] == summary["runs"]["main"]["l2_vs_T213dt22_by_day"][day],
                "New T213 comparison score differs")
        difference = item["field_l2"] / item["archived_l2"] - 1
        require(item["relative_difference"] == difference and abs(difference) < 6e-6,
                "T213 comparison exceeds its reported discrepancy")
    assets = json.loads((DATA / "field-assets.json").read_text())
    require(len(assets["assets"]) == 6 and
            sum(a["bytes"] for a in assets["assets"]) == assets["total_npz_bytes"],
            "Field asset count or total size differs")
    return len(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true",
                        help="Regenerate derived tables after verifying the evidence hashes.")
    args = parser.parse_args()
    try:
        path = DATA / "manifest.json"
        manifest = json.loads(path.read_text())
        for name, info in manifest["evidence"].items():
            require(sha256(EVIDENCE / name) == info["sha256"], f"Evidence hash differs: {name}")
            require((EVIDENCE / name).stat().st_size == info["bytes"], f"Evidence size differs: {name}")
            require(len(records(name)) == info["records"], f"Record count differs: {name}")
        require(set(manifest["evidence"]) == {p.name for p in EVIDENCE.glob("*.jsonl")},
                "The evidence file set differs from the manifest")
        for name, info in manifest.get("supplemental", {}).items():
            require(sha256(DATA / name) == info["sha256"], f"Supplemental input hash differs: {name}")
            require((DATA / name).stat().st_size == info["bytes"], f"Supplemental input size differs: {name}")
        figures = json.loads((ROOT / "figures/manifest.json").read_text())
        listed_figures = [item for visual in figures["visuals"].values() for item in visual["outputs"]]
        names = [item["file"] for item in listed_figures]
        require(len(names) == len(set(names)), "Duplicate figure output in the manifest")
        require("relaxation-fields.gif" in names, "The inline correction GIF is missing from the manifest")
        bundled = {p.name for p in (ROOT / "figures").iterdir()
                   if p.suffix.lower() in {".png", ".svg", ".mp4", ".gif"}}
        require(set(names) == bundled, "The figure file set differs from the manifest")
        for visual in figures["visuals"].values():
            for item in visual["outputs"]:
                require("sha256" in item and "bytes" in item, f"Figure hash or size missing: {item['file']}")
                require(sha256(ROOT / "figures" / item["file"]) == item["sha256"],
                        f"Figure hash differs: {item['file']}")
                require((ROOT / "figures" / item["file"]).stat().st_size == item["bytes"],
                        f"Figure size differs: {item['file']}")
        for name, info in figures.get("qa_files", {}).items():
            require(sha256(ROOT / name) == info["sha256"], f"Figure check hash differs: {name}")
        field_count = check_fields()
        tables, sources, claims = derive()
        if args.write:
            for name, content in tables.items():
                (DATA / name).write_bytes(content)
            manifest["derived"] = {name: {"sha256": sha256(DATA / name), "bytes": len(content)}
                                   for name, content in tables.items()}
            manifest["selectors"] = sources
            path.write_bytes(json_bytes(manifest))
        else:
            require(manifest["selectors"] == sources, "Manifest selectors differ from the selected records")
            for name, content in tables.items():
                require((DATA / name).read_bytes() == content, f"Derived table differs: {name}; run --write")
                require(sha256(DATA / name) == manifest["derived"][name]["sha256"],
                        f"Derived hash differs: {name}")
        print(f"PASS: {len(manifest['evidence'])} evidence files, {len(sources)} unique selectors, {len(tables)} derived files")
        print(f"PASS: {field_count} field measurements, run summaries, independent-reference ratios and figure hashes")
        print(f"Day 10 vs T85: cones {claims['cones_day10_reduction']['value']:.11f}x; "
              f"controls {claims['controls_day10_reduction']['value']:.11f}x")
        print("Main/spectral: " + ", ".join(f"{claims[f'main_over_spectral_M{m}']['value']:.9f}" for m in (128,256,512)))
        print("TRSK/spectral: " + ", ".join(f"{claims[f'trsk_over_spectral_M{m}']['value']:.9f}" for m in (128,256,512)))
        print(f"Restart, three elapsed days: {claims['restart_elapsed3_ratio']['value']:.10f}x")
        print("This check uses saved records and field measurements; it does not run a simulation.")
    except (ValueError, KeyError, OSError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
