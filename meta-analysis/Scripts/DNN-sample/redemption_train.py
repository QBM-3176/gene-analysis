"""Train the deep neural network on all.xlsx to predict Disease_Label.

Features: Compartment (one-hot), Sex (M = 1, F = 0), C4BPA, COL13A1, HABP2, CFH, CFI.
Validation: LeaveOneGroupOut by StudyID (plan section 6) - each study is held out once
and predicted by a network trained on the others. Outputs go to Results/ML_training/.
"""
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import LeaveOneGroupOut

sys.path.insert(0, str(Path(__file__).resolve().parent))
from deep_neural_network import DeepNeuralNetworkRunner, TARGET

cfg = {
    "project_dir": Path("~/Striki/Quant/quant.2/MetaAnalysis").expanduser(),
    "group_col": "StudyID",
    "epochs": 50,
}
cfg["in_file"] = cfg["project_dir"] / "Data" / "ML_input" / "all.xlsx"
cfg["out_dir"] = cfg["project_dir"] / "Results" / "ML_training"

data = pd.read_excel(cfg["in_file"], sheet_name="training_data")
print(f"Loaded {len(data):,} samples from {data[cfg['group_col']].nunique()} studies")

predictions, per_study = [], []
for train_idx, test_idx in LeaveOneGroupOut().split(data, groups=data[cfg["group_col"]]):
    train, test = data.iloc[train_idx], data.iloc[test_idx]
    study = test[cfg["group_col"]].iloc[0]
    print(f"Held out {study}: {len(train):,} training, {len(test):,} test samples")

    runner = DeepNeuralNetworkRunner(train, test)
    runner.setup()
    runner.train(epochs=cfg["epochs"], verbose=False)
    prob = runner.inference(test)

    y = test[TARGET].to_numpy()
    # AUC is undefined when the held-out study has only one class
    auc = roc_auc_score(y, prob) if len(np.unique(y)) == 2 else np.nan
    per_study.append({"HeldOutStudy": study, "n_control": int((y == 0).sum()),
                      "n_case": int((y == 1).sum()), "AUC": auc,
                      "Accuracy": float(((prob > 0.5) == y).mean())})
    predictions.append(test[["GSM", cfg["group_col"], "Compartment", TARGET]].assign(Probability=prob))

predictions = pd.concat(predictions, ignore_index=True)
per_study = pd.DataFrame(per_study)
summary = pd.DataFrame([{
    "n_studies": len(per_study), "n_studies_with_AUC": int(per_study["AUC"].notna().sum()),
    "mean_study_AUC": per_study["AUC"].mean(),
    "pooled_AUC": roc_auc_score(predictions[TARGET], predictions["Probability"]),
    "pooled_accuracy": float(((predictions["Probability"] > 0.5) == predictions[TARGET]).mean()),
}])

for name, df in {"dnn_lodo_predictions": predictions, "dnn_lodo_per_study": per_study,
                 "dnn_lodo_summary": summary}.items():
    df.to_excel(cfg["out_dir"] / f"{name}.xlsx", index=False)
(cfg["out_dir"] / "dnn_lodo_sessionInfo.txt").write_text(
    f"Python {platform.python_version()}\ntorch {torch.__version__}\n"
    f"pandas {pd.__version__}\nnumpy {np.__version__}\nseed 42\nepochs {cfg['epochs']}\n")

print(f"Saved predictions, per-study metrics")
pd.set_option("display.width", 200)
print(per_study.round(3).to_string(index=False))

print("\nSummary:")
print(summary.round(3).to_string(index=False))
