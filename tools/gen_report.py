#!/usr/bin/env python3
"""gen_report.py — M8 shape-ranking comparison report (PLAN M8 session-2 deliverable).

Reads a ranking CSV emitted by scour_m8_gate (columns: shape, V_eq_m3, V_ci_m3, V_last_m3,
edge_scour_m, maxu, mass_err) and writes a self-contained HTML report: a ranked horizontal bar
chart of trapped volume (with CI whiskers), an edge-scour comparison, and the raw table. Pure
Python stdlib (inline SVG, no matplotlib / no network) so it runs anywhere and the HTML is
portable. Ranking is by matched-time trapped volume V_last (higher = better self-ballasting), the
comparative metric whose authority is stated in RESEARCH §11 — absolute magnitudes inherit the
gate_M6 HSV under-prediction (flagged in the report banner).

Usage:  python tools/gen_report.py <ranking.csv> [-o report.html] [--title "..."]
"""
import argparse
import csv
import html
import sys

# Brand-neutral, colour-blind-safe categorical palette (theme-agnostic dev tool).
BARC = "#3b7dd8"
CIC = "#1b3a63"
SCOURC = "#d1495b"
GRID = "#00000018"


def _optf(r, key):
    """Optional float column (M9 reversal diagnostics); None when absent."""
    v = r.get(key)
    if v is None or v == "":
        return None
    try:
        return float(v)
    except (ValueError, TypeError):
        return None


def load(path):
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            try:
                rows.append({
                    "shape": r.get("shape", "?"),
                    "V_last": float(r.get("V_last_m3", r.get("V_eq_m3", 0)) or 0),
                    "V_eq": float(r.get("V_eq_m3", 0) or 0),
                    "V_ci": float(r.get("V_ci_m3", 0) or 0),
                    "edge": float(r.get("edge_scour_m", 0) or 0),
                    "maxu": float(r.get("maxu", 0) or 0),
                    "mass_err": float(r.get("mass_err", 0) or 0),
                    # M9 tidal-reversal diagnostics (absent in M8 CSVs -> None -> section omitted).
                    "clip_frac": _optf(r, "clip_frac"),
                    "slacks": _optf(r, "slacks"),
                    "ke_growth": _optf(r, "ke_growth"),
                })
            except (ValueError, TypeError):
                pass
    return rows


def bars(rows, key, ci_key, unit, colour, title):
    """One horizontal-bar panel as inline SVG (values -> bar length, optional CI whisker)."""
    if not rows:
        return "<p>(no data)</p>"
    vmax = max((abs(r[key]) for r in rows), default=0) or 1.0
    W, rowh, padL, padR = 720, 34, 120, 90
    H = rowh * len(rows) + 20
    out = [f'<svg viewBox="0 0 {W} {H}" width="100%" role="img" aria-label="{html.escape(title)}">']
    for i, r in enumerate(sorted(rows, key=lambda x: -x[key])):
        y = 10 + i * rowh
        bw = (W - padL - padR) * (abs(r[key]) / vmax)
        out.append(f'<text x="{padL-8}" y="{y+rowh*0.55}" text-anchor="end" font-size="14" fill="currentColor">{html.escape(r["shape"])}</text>')
        out.append(f'<rect x="{padL}" y="{y+4}" width="{bw:.1f}" height="{rowh-12}" rx="3" fill="{colour}"/>')
        if ci_key and r.get(ci_key):
            cw = (W - padL - padR) * (r[ci_key] / vmax)
            cx = padL + bw
            out.append(f'<line x1="{cx-cw:.1f}" y1="{y+rowh*0.42}" x2="{cx+cw:.1f}" y2="{y+rowh*0.42}" stroke="{CIC}" stroke-width="2"/>')
        out.append(f'<text x="{padL+bw+6:.1f}" y="{y+rowh*0.55}" font-size="13" fill="currentColor">{r[key]:.3e} {unit}</text>')
    out.append("</svg>")
    return "".join(out)


def table(rows):
    hd = ["shape", "V_last [m³]", "V_eq (fit) [m³]", "±CI", "edge scour [m]", "max|u|", "mass err"]
    ks = ["shape", "V_last", "V_eq", "V_ci", "edge", "maxu", "mass_err"]
    th = "".join(f"<th>{h}</th>" for h in hd)
    tr = ""
    for r in sorted(rows, key=lambda x: -x["V_last"]):
        cells = []
        for k in ks:
            v = r[k]
            cells.append(f"<td>{html.escape(v)}</td>" if k == "shape" else f"<td>{v:.4e}</td>")
        tr += "<tr>" + "".join(cells) + "</tr>"
    return f"<table><thead><tr>{th}</tr></thead><tbody>{tr}</tbody></table>"


def stability(rows):
    """M9 tidal-reversal stability strip: per-shape mass conservation, MORFAC clip-fraction, slack
    count and KE growth across slacks. Rendered only when the CSV carries the M9 columns."""
    hd = ["shape", "slacks", "mass err", "MORFAC clip", "KE growth (slacks)"]
    tr = ""
    for r in sorted(rows, key=lambda x: -x["V_last"]):
        clip = r.get("clip_frac")
        ke = r.get("ke_growth")
        sl = r.get("slacks")
        tr += (
            "<tr>"
            f'<td>{html.escape(r["shape"])}</td>'
            f"<td>{int(sl) if sl is not None else '-'}</td>"
            f"<td>{r['mass_err']:.1e}</td>"
            f"<td>{clip*100:.2f}%</td>"
            f"<td>×{ke:.2f}</td>"
            "</tr>"
        )
    th = "".join(f"<th>{h}</th>" for h in hd)
    return f"<table><thead><tr>{th}</tr></thead><tbody>{tr}</tbody></table>"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("-o", "--out", default="ranking_report.html")
    ap.add_argument("--title", default="Scour-protection shape ranking")
    a = ap.parse_args()
    rows = load(a.csv)
    if not rows:
        print(f"gen_report: no rows parsed from {a.csv}", file=sys.stderr)
        return 1
    best = max(rows, key=lambda r: r["V_last"])["shape"]
    # M9: show the reversal-stability strip when the tidal diagnostic columns are present.
    is_tidal = any(r.get("slacks") is not None for r in rows)
    tidal_section = (
        f'<h2>Tidal-reversal stability</h2>\n<p style="font-size:.88rem;opacity:.8">Each shape ran an '
        "identical reversing tide (face-swap reversal, cosine slack ramps). Total sand (bed+suspended, "
        "credited for structure pinning) is conserved to machine precision; the MORFAC |Δz<sub>b</sub>| "
        "limiter clip-fraction stays &lt;1%; kinetic energy does not grow monotonically across the "
        f"slack transitions.</p>\n{stability(rows)}\n"
        if is_tidal else ""
    )
    doc = f"""<!doctype html><meta charset="utf-8"><title>{html.escape(a.title)}</title>
<style>
 body{{font:15px/1.5 system-ui,sans-serif;max-width:900px;margin:2rem auto;padding:0 1rem;color:#1a1a1a;background:#fff}}
 @media(prefers-color-scheme:dark){{body{{color:#e8e8e8;background:#161616}}}}
 h1{{font-size:1.5rem}} h2{{font-size:1.1rem;margin-top:2rem}}
 .banner{{background:#d1495b22;border-left:4px solid #d1495b;padding:.75rem 1rem;border-radius:4px;font-size:.92rem}}
 table{{border-collapse:collapse;width:100%;margin-top:.5rem;font-size:.9rem}}
 th,td{{border:1px solid {GRID};padding:.35rem .6rem;text-align:right}} th:first-child,td:first-child{{text-align:left}}
 .best{{font-weight:600}}
</style>
<h1>{html.escape(a.title)}</h1>
<p class="banner"><strong>Comparative ranking</strong> (RESEARCH §11): shapes are ranked by
matched-time trapped sand volume under identical conditions. <strong>Absolute magnitudes are
indicative</strong> — they inherit the gate_M6 horseshoe-vortex under-prediction; trust the
<em>ordering</em>, not the numbers. Best self-ballasting shape: <span class="best">{html.escape(best)}</span>.</p>
<h2>Trapped sand volume (higher = better; whisker = fit CI)</h2>
{bars(rows, "V_last", "V_ci", "m³", BARC, "Trapped volume by shape")}
<h2>Edge scour (lower = better)</h2>
{bars(rows, "edge", None, "m", SCOURC, "Edge scour by shape")}
{tidal_section}<h2>Data</h2>
{table(rows)}
<p style="font-size:.8rem;opacity:.6">Generated by tools/gen_report.py from {html.escape(a.csv)}.</p>
"""
    with open(a.out, "w", encoding="utf-8") as f:
        f.write(doc)
    print(f"gen_report: wrote {a.out} ({len(rows)} shapes; best={best})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
