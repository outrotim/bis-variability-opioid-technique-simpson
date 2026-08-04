#!/usr/bin/env python3
"""Generate aggregate-only public figures from aggregate_estimates.json."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch


COLORS = {
    "ink": "#2F3944",
    "line": "#7A838D",
    "muted": "#66717D",
    "green": "#009E73",
    "tiva": "#0072B2",
    "volatile": "#D55E00",
    "light": "#F4F6F8",
}


def load_estimates(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def save_all(fig, outdir: Path, stem: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for extension in ["png", "pdf", "svg"]:
        fig.savefig(
            outdir / f"{stem}.{extension}",
            dpi=600,
            facecolor="white",
            bbox_inches="tight",
        )


def flow_box(ax, x, y, width, height, title, value, edge):
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.006,rounding_size=0.009",
        facecolor="white",
        edgecolor=edge,
        linewidth=1.0,
        transform=ax.transAxes,
    )
    ax.add_patch(patch)
    ax.text(x + 0.025, y + 0.61 * height, title, transform=ax.transAxes,
            fontsize=8.0, fontweight="bold", color=edge, va="center")
    ax.text(x + 0.025, y + 0.27 * height, f"n = {value:,}", transform=ax.transAxes,
            fontsize=7.2, color=COLORS["ink"], va="center")


def draw_flow(data: dict, outdir: Path) -> None:
    flow = data["cohort_flow"]
    fig = plt.figure(figsize=(170 / 25.4, 112 / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")

    x, width, height = 0.12, 0.56, 0.085
    y_values = [0.86, 0.71, 0.56, 0.41, 0.26]
    keys = [
        "source_cohort",
        "initially_eligible",
        "bis_metrics_processed",
        "primary_analysis",
        "mac_available_subset",
    ]
    for y, key in zip(y_values, keys):
        edge = COLORS["green"] if key in {"primary_analysis", "mac_available_subset"} else COLORS["ink"]
        flow_box(ax, x, y, width, height, flow[key]["label"], flow[key]["n"], edge)
    center = x + width / 2
    for upper, lower in zip(y_values[:-1], y_values[1:]):
        ax.annotate(
            "",
            xy=(center, lower + height),
            xytext=(center, upper),
            xycoords=ax.transAxes,
            arrowprops=dict(arrowstyle="-|>", color=COLORS["line"], lw=0.9),
        )

    exclusions = flow["exclusions"]
    notes = [
        (y_values[1], "Excluded before eligibility", exclusions["before_eligibility"]),
        (y_values[2], "Insufficient BIS data", exclusions["insufficient_bis_data"]),
        (y_values[3], "No positive opioid / missing ASA / signal quality", 196 + 96 + 5),
        (y_values[4], "Missing TWA-MAC", exclusions["missing_twa_mac"]),
    ]
    for y, label, number in notes:
        ax.plot([0.70, 0.73], [y + height / 2] * 2, transform=ax.transAxes,
                color=COLORS["line"], lw=0.8)
        ax.text(0.745, y + height / 2, f"{label}, n = {number:,}", transform=ax.transAxes,
                fontsize=7.0, color=COLORS["muted"], va="center")

    branch_y = 0.08
    flow_box(ax, 0.12, branch_y, 0.25, height, "TIVA", flow["tiva"]["n"], COLORS["tiva"])
    flow_box(ax, 0.43, branch_y, 0.32, height, "Volatile-supplemented", flow["volatile_supplemented"]["n"], COLORS["volatile"])
    ax.plot([center, center], [y_values[-1], 0.20], transform=ax.transAxes, color=COLORS["line"], lw=0.9)
    ax.plot([0.245, 0.59], [0.20, 0.20], transform=ax.transAxes, color=COLORS["line"], lw=0.9)
    for branch_center in [0.245, 0.59]:
        ax.annotate("", xy=(branch_center, branch_y + height), xytext=(branch_center, 0.20),
                    xycoords=ax.transAxes,
                    arrowprops=dict(arrowstyle="-|>", color=COLORS["line"], lw=0.9))

    save_all(fig, outdir, "figure1_flow_public")
    plt.close(fig)


def p_label(value: float) -> str:
    return "P<0.001" if value < 0.001 else f"P={value:.3f}"


def draw_forest(data: dict, outdir: Path) -> None:
    rows = data["core_models"]
    fig, axes = plt.subplots(1, 2, figsize=(170 / 25.4, 102 / 25.4), sharex=True)
    groups = [("OME", "A"), ("Remifentanil only", "B")]
    for ax, (outcome, panel) in zip(axes, groups):
        selected = [row for row in rows if row["outcome"] == outcome]
        y = np.arange(len(selected))[::-1]
        for position, row in zip(y, selected):
            group = row["group"]
            color = COLORS["ink"]
            marker = "D"
            if group == "technique_adjusted":
                color, marker = COLORS["green"], "s"
            elif group == "tiva":
                color, marker = COLORS["tiva"], "o"
            elif group == "volatile":
                color, marker = COLORS["volatile"], "^"
            elif group == "pooled_subset":
                color, marker = COLORS["line"], "D"
            ax.errorbar(
                row["ratio"], position,
                xerr=[[row["ratio"] - row["ci_low"]], [row["ci_high"] - row["ratio"]]],
                fmt=marker, color=color, ecolor=color, markersize=5, capsize=2.5, linewidth=1.1,
                markerfacecolor="white" if group == "pooled_subset" else color,
            )
            ax.text(1.031, position, f"{row['ratio']:.3f} ({row['ci_low']:.3f}-{row['ci_high']:.3f}); {p_label(row['p'])}",
                    va="center", fontsize=6.4, color=COLORS["ink"])
        ax.axvline(1, color=COLORS["line"], lw=0.8, linestyle="--")
        ax.set_yticks(y)
        ax.set_yticklabels([f"{row['label'].split(', ', 1)[1]}\n(n={row['n']:,})" for row in selected], fontsize=7)
        ax.set_xlim(0.988, 1.030)
        ax.set_xlabel("Ratio per 1-point higher CV-BIS", fontsize=7.5)
        ax.tick_params(axis="x", labelsize=7)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.set_title(f"{panel}   {outcome} outcome", loc="left", fontsize=9, fontweight="bold")
    fig.subplots_adjust(left=0.22, right=0.74, wspace=0.62, bottom=0.17, top=0.88)
    save_all(fig, outdir, "figure2_forest_public")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--estimates", type=Path, default=Path("aggregate_estimates.json"))
    parser.add_argument("--outdir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    data = load_estimates(args.estimates)
    draw_flow(data, args.outdir)
    draw_forest(data, args.outdir)
    print(f"Saved public figures to: {args.outdir}")


if __name__ == "__main__":
    main()
