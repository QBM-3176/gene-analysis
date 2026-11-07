#!/usr/bin/env python3
"""RenoRisk: 2-D view of the ML input workbooks (Data/ML_input/all.xlsx, Results/ML_training/all_bn.xlsx).

Projects the five panel-gene expression values to two dimensions with PCA and
t-SNE (features standardised over all samples) and colours each sample by
Disease_Label: blue = healthy control, red = diseased. One figure per workbook:

- all.xlsx: uncorrected across studies, and two cohorts are per-gene centred,
  so between-study differences dominate both projections.
- all_bn.xlsx: Stage23 removeBatchEffect values. The correction
  protected Disease in its design, so any class separation here is partly
  shaped by the labels.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

# 1. Configuration.
cfg = {
    "project_dir": Path("~/Striki/Quant/quant.2/MetaAnalysis").expanduser(),
    "genes": ["C4BPA", "COL13A1", "HABP2", "CFH", "CFI"],
    "random_state": 42,
}
cfg["out_dir"] = cfg["project_dir"] / "Results" / "ML_training"
# input workbook -> (output figure, title note)
cfg["inputs"] = {
    cfg["project_dir"] / "Data" / "ML_input" / "all.xlsx": ("ML_data_2d.png", "not batch-corrected"),
    cfg["out_dir"] / "all_bn.xlsx": ("ML_data_bn_2d.png",
                                  "batch-corrected (removeBatchEffect)"),
}
groups = [
    (0, "Healthy control", "#1f77b4"),
    (1, "Diseased", "#d62728"),
]

# 2. One figure per workbook.
for in_file, (out_png, note) in cfg["inputs"].items():
    data = pd.read_excel(in_file, sheet_name="training_data")
    assert data[cfg["genes"]].notna().all().all(), "missing expression values"
    assert set(data["Disease_Label"]) <= {0, 1}, "unexpected Disease_Label values"

    X = StandardScaler().fit_transform(data[cfg["genes"]].to_numpy())
    y = data["Disease_Label"].to_numpy()

    # 2b. Two 2-D embeddings.
    pca = PCA(n_components=2, random_state=cfg["random_state"])
    pca_xy = pca.fit_transform(X)
    tsne_xy = TSNE(n_components=2, perplexity=30, init="pca",
                   random_state=cfg["random_state"]).fit_transform(X)

    var = pca.explained_variance_ratio_ * 100
    panels = [
        (pca_xy, "PCA", f"PC1 ({var[0]:.1f}%)", f"PC2 ({var[1]:.1f}%)"),
        (tsne_xy, "t-SNE", "t-SNE 1", "t-SNE 2"),
    ]

    # 2c. Plot: controls drawn last so the minority class stays visible.
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), constrained_layout=True)
    for ax, (xy, title, xlab, ylab) in zip(axes, panels):
        for label, name, colour in sorted(groups, key=lambda g: -g[0]):
            mask = y == label
            ax.scatter(xy[mask, 0], xy[mask, 1], s=14, alpha=0.7, c=colour,
                       edgecolors="none", label=f"{name} (n = {mask.sum()})")
        ax.set_title(title)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.legend(frameon=False, loc="best")
    fig.suptitle(
        f"RenoRisk ML input: {len(data)} samples, 5-gene panel "
        f"({', '.join(cfg['genes'])}), standardised, {note}")
    fig.savefig(cfg["out_dir"] / out_png, dpi=150)
    plt.close(fig)
    print(f"wrote {cfg['out_dir'] / out_png}")
