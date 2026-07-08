#!/usr/bin/env python3
"""Generate aggregate-only public figures from aggregate_estimates.json."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np


COLORS = {
    "ink": "#26364A",
    "main_light": "#F5F8FB",
    "line": "#6B7280",
    "muted": "#697386",
    "muted_light": "#AEB6C2",
    "green": "#009E73",
    "tiva": "#0072B2",
    "volatile": "#D55E00",
    "raw": "#5A5A5A",
}


def load_estimates(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def add_card(ax, x, y, w, h, title, body, color, fill="white", title_size=6.9, body_size=6.05):
    rect = patches.FancyBboxPatch(
        (x, y - h / 2),
        w,
        h,
        boxstyle="round,pad=0.006,rounding_size=0.012",
        facecolor=fill,
        edgecolor=color,
        linewidth=0.9,
        transform=ax.transAxes,
    )
    ax.add_patch(rect)
    ax.add_patch(
        patches.Rectangle(
            (x + 0.008, y - h / 2 + 0.010),
            0.008,
            h - 0.020,
            facecolor=color,
            edgecolor=color,
            linewidth=0,
            transform=ax.transAxes,
        )
    )
    ax.text(x + 0.034, y + 0.18 * h, title, ha="left", va="center",
            fontsize=title_size, fontweight="semibold", color=color, transform=ax.transAxes)
    ax.text(x + 0.034, y - 0.20 * h, body, ha="left", va="center",
            fontsize=body_size, color="#1F2937", linespacing=1.16, transform=ax.transAxes)


def arrow(ax, start, end, lw=0.75):
    ax.annotate("", xy=end, xytext=start, xycoords=ax.transAxes,
                arrowprops=dict(arrowstyle="-|>", color=COLORS["line"], lw=lw, mutation_scale=6.5))


def add_note(ax, x, y, title, body):
    ax.plot([x - 0.030, x - 0.010], [y, y], transform=ax.transAxes,
            color=COLORS["muted_light"], lw=0.55, solid_capstyle="round")
    ax.text(x, y + 0.034, title, ha="left", va="center", fontsize=5.45,
            fontweight="semibold", color=COLORS["muted"], transform=ax.transAxes)
    ax.text(x, y - 0.018, body, ha="left", va="top", fontsize=4.95,
            color=COLORS["muted"], linespacing=1.12, transform=ax.transAxes)


def draw_flow(data: dict, outdir: Path) -> None:
    flow = data["cohort_flow"]
    fig = plt.figure(figsize=(84 / 25.4, 116 / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    main_x, main_w, note_x, node_h = 0.055, 0.590, 0.705, 0.083
    ys = [0.930, 0.772, 0.614, 0.456, 0.298]
    labels = [
        ("source_cohort", "n = {n:,} noncardiac cases\n2016-2017", COLORS["ink"], "white"),
        ("initially_eligible", "n = {n:,}", COLORS["ink"], COLORS["main_light"]),
        ("adequate_bis_variability", "n = {n:,}", COLORS["ink"], COLORS["main_light"]),
        ("primary_complete_case", "n = {n:,}\nOME >0 with non-missing ASA", COLORS["green"], "white"),
        ("mac_available_subset", "n = {n:,}\nTechnique-stratified analyses", COLORS["green"], "white"),
    ]
    for y, (key, body_template, color, fill) in zip(ys, labels):
        node = flow[key]
        add_card(ax, main_x, y, main_w, node_h, node["label"], body_template.format(n=node["n"]), color, fill=fill)
    center = main_x + main_w / 2
    for y1, y2 in zip(ys[:-1], ys[1:]):
        arrow(ax, (center, y1 - node_h / 2 - 0.004), (center, y2 + node_h / 2 + 0.004))

    exc = flow["exclusions"]
    add_note(ax, note_x, ys[1], "Excluded before eligibility",
             "Non-general anesthesia, n=345\nNo BIS monitoring, n=425\nAge <18 yr, n=48\nDuration/ASA criteria, n=13")
    add_note(ax, note_x, ys[2], "Excluded", "Insufficient BIS data for\nvariability computation,\nn=121")
    add_note(ax, note_x, ys[3], "Primary exclusions",
             "No opioid recorded, n=255\nASA missing among opioid-\npositive cases, n=96")
    add_note(ax, note_x, ys[4], "Excluded", f"Missing TWA-MAC, n={exc['missing_twa_mac']['n']}")

    fork_y, split_y, split_h = 0.212, 0.142, 0.088
    tiva_x, tiva_w = 0.055, 0.280
    vol_x, vol_w = 0.372, 0.348
    tiva_center, vol_center = tiva_x + tiva_w / 2, vol_x + vol_w / 2
    ax.plot([tiva_center, vol_center], [fork_y, fork_y], transform=ax.transAxes,
            color=COLORS["line"], lw=0.75, solid_capstyle="round")
    arrow(ax, (center, ys[4] - node_h / 2 - 0.004), (center, fork_y))
    arrow(ax, (tiva_center, fork_y), (tiva_center, split_y + split_h / 2 + 0.004))
    arrow(ax, (vol_center, fork_y), (vol_center, split_y + split_h / 2 + 0.004))
    add_card(ax, tiva_x, split_y, tiva_w, split_h, "TIVA", "TWA-MAC <0.3\nn = 3,086",
             COLORS["tiva"], fill="#F7FBFF", title_size=6.45, body_size=5.75)
    add_card(ax, vol_x, split_y, vol_w, split_h, "Volatile-supplemented", "TWA-MAC >=0.3\nn = 1,979",
             COLORS["volatile"], fill="#FFF8F2", title_size=5.85, body_size=5.65)
    ax.text(0.055, 0.040, "Strata were defined using time-weighted average MAC (TWA-MAC).",
            ha="left", va="center", fontsize=5.15, color=COLORS["muted"], transform=ax.transAxes)

    for ext in ["png", "pdf", "svg"]:
        fig.savefig(outdir / f"figure1_flow_public.{ext}", dpi=600, facecolor="white", edgecolor="none")
    plt.close(fig)


def p_label(p_value: float) -> str:
    if p_value < 0.001:
        return "P < 0.001"
    return f"P = {p_value:.3f}"


def draw_forest(data: dict, outdir: Path) -> None:
    rows = data["core_models"]
    fig, ax = plt.subplots(figsize=(6.8, 3.9))
    y = np.arange(len(rows))[::-1]
    for idx, row in enumerate(rows):
        color = COLORS["raw"] if row["outcome"].startswith("Raw") else (
            COLORS["tiva"] if "TIVA" in row["label"] else COLORS["volatile"] if "Volatile" in row["label"] else COLORS["ink"]
        )
        marker = "^" if row["outcome"].startswith("Raw") else "o"
        if "Pooled" in row["label"] or "pooled" in row["label"]:
            marker = "D"
        ax.errorbar(
            row["exp_beta"],
            y[idx],
            xerr=[[row["exp_beta"] - row["ci_low"]], [row["ci_high"] - row["exp_beta"]]],
            fmt=marker,
            color=color,
            ecolor=color,
            elinewidth=1.2,
            capsize=3,
            markersize=5.2,
        )
    ax.axvline(1.0, color="#888888", linestyle="--", linewidth=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{row['label']} (n={row['n']})" for row in rows], fontsize=8.5)
    ax.set_xlabel("exp(beta) per 1% higher CV-BIS")
    ax.set_xlim(0.988, 1.026)
    ax.grid(axis="x", alpha=0.16)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title("Aggregate public forest plot", loc="left", fontweight="bold")
    for idx, row in enumerate(rows):
        ax.text(1.0265, y[idx], f"{row['exp_beta']:.3f} [{row['ci_low']:.3f}, {row['ci_high']:.3f}], {p_label(row['p'])}",
                ha="left", va="center", fontsize=7.7, color="#333333")
    fig.subplots_adjust(left=0.38, right=0.72, top=0.88, bottom=0.16)
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(outdir / f"figure2_forest_public.{ext}", dpi=600, facecolor="white", edgecolor="none")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--estimates", type=Path, default=Path("aggregate_estimates.json"))
    parser.add_argument("--outdir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    data = load_estimates(args.estimates)
    draw_flow(data, args.outdir)
    draw_forest(data, args.outdir)
    print(f"Saved public figures to: {args.outdir}")


if __name__ == "__main__":
    main()
