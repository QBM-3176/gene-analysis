import numpy as np
import pandas as pd
from tqdm.auto import tqdm
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from torch.optim.lr_scheduler import CosineAnnealingLR

GENES = ["C4BPA", "COL13A1", "HABP2", "CFH", "CFI"]
COMPARTMENTS = ["Glomerular", "Tubulointerstitial", "Whole biopsy"]
SEX_CODES = {"M": 1.0, "F": 0.0}  # missing sex -> 0.5 (unknown)
TARGET = "Disease_Label"


def encode_features(df):
    """Numeric feature matrix: one-hot Compartment, Sex as 0/1, the five panel genes."""
    compartment = pd.get_dummies(
        pd.Categorical(df["Compartment"], categories=COMPARTMENTS), prefix="Compartment", dtype=float)
    sex = df["Sex"].map(SEX_CODES).fillna(0.5).to_numpy(float)[:, None]
    return np.hstack([compartment.to_numpy(), sex, df[GENES].to_numpy(float)])


class ResidualBlock(nn.Module):
    def __init__(self, hidden_size, dropout_prob):
        super(ResidualBlock, self).__init__()
        self.block = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.ReLU(),
            nn.Dropout(dropout_prob),
            nn.Linear(hidden_size, hidden_size),
            nn.LayerNorm(hidden_size),
        )
        self.relu = nn.ReLU()

    def forward(self, x):
        residual = x
        out = self.block(x)
        out += residual  # Skip connection
        return self.relu(out)


class DeepNeuralNetwork(nn.Module):
    # Sized for ~850 samples x 9 features; the original 10 x 4096 had ~150M parameters
    def __init__(self, input_size, num_layers=4, hidden_size=64, dropout_prob=0.2):
        super(DeepNeuralNetwork, self).__init__()

        # First layer
        self.input_layer = nn.Sequential(
            nn.Linear(input_size, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.ReLU(),
            nn.Dropout(dropout_prob),
        )

        # Residual blocks
        self.residual_blocks = nn.ModuleList()
        for i in range(num_layers - 2):
            self.residual_blocks.append(ResidualBlock(hidden_size, dropout_prob))

        # Output layer: one logit for Disease_Label
        self.output_layer = nn.Linear(hidden_size, 1)

    def forward(self, x):
        x = self.input_layer(x)

        for block in self.residual_blocks:
            x = block(x)

        return self.output_layer(x)


class DeepNeuralNetworkRunner:
    """Binary classifier of Disease_Label from rows of all.xlsx (pandas DataFrames)."""

    def __init__(self, train, val):
        self.train_data = train
        self.val_data = val
        self.x_mean = None
        self.x_std = None
        self.model = None
        self.device = None
        self.loss_function = None
        self.optimizer = None
        self.train_dataset = None
        self.train_loader = None

        np.random.seed(42)
        torch.manual_seed(42)
        torch.cuda.manual_seed(42)

    def transform(self, df):
        """Encode and standardise with the training-set mean/SD (no peeking at val/test)."""
        return torch.FloatTensor((encode_features(df) - self.x_mean) / self.x_std)

    def setup(self):
        X_train_np = encode_features(self.train_data)
        self.x_mean = X_train_np.mean(axis=0)
        self.x_std = X_train_np.std(axis=0)
        self.x_std[self.x_std == 0] = 1.0  # e.g. a compartment absent from training

        self.X_train = self.transform(self.train_data)
        self.y_train = torch.FloatTensor(self.train_data[TARGET].to_numpy(float)).unsqueeze(1)
        self.X_val = self.transform(self.val_data)
        self.y_val = torch.FloatTensor(self.val_data[TARGET].to_numpy(float)).unsqueeze(1)

        self.model = DeepNeuralNetwork(self.X_train.shape[1])
        total_params = sum(p.numel() for p in self.model.parameters() if p.requires_grad)
        print(f"Deep Neural Network created with {total_params:,} parameters")

        if torch.cuda.is_available():
            self.device = torch.device("cuda")
        elif torch.backends.mps.is_available():
            self.device = torch.device("mps")
        else:
            self.device = torch.device("cpu")

        print(f"Using {self.device}")

        self.model.to(self.device)
        # Weight the positives so the controls (~1 in 5 samples) are not ignored
        n_pos = self.y_train.sum()
        pos_weight = (len(self.y_train) - n_pos) / n_pos.clamp(min=1)
        self.loss_function = nn.BCEWithLogitsLoss(pos_weight=pos_weight.reshape(1).to(self.device))
        self.optimizer = optim.AdamW(self.model.parameters(), lr=0.001, weight_decay=0.01)

        self.train_dataset = TensorDataset(self.X_train, self.y_train)
        self.train_loader = DataLoader(self.train_dataset, batch_size=64, shuffle=True)

    def train(self, epochs=50, verbose=True):
        scheduler = CosineAnnealingLR(self.optimizer, T_max=epochs, eta_min=0)
        for epoch in range(1, epochs + 1):
            self.model.train()
            train_losses = []

            for batch_X, batch_y in tqdm(self.train_loader, disable=not verbose, leave=False):
                batch_X = batch_X.to(self.device)
                batch_y = batch_y.to(self.device)

                self.optimizer.zero_grad()
                outputs = self.model(batch_X)
                loss = self.loss_function(outputs, batch_y)
                loss.backward()

                # Gradient clipping
                torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)

                self.optimizer.step()
                train_losses.append(loss.item())

            # Validation
            self.model.eval()
            with torch.no_grad():
                val_outputs = self.model(self.X_val.to(self.device))
                val_loss = self.loss_function(val_outputs, self.y_val.to(self.device))
                val_pred = (val_outputs > 0).float()
                accuracy = (val_pred == self.y_val.to(self.device)).float().mean()

            if verbose:
                avg_train_loss = np.mean(train_losses)
                print(f"Epoch [{epoch}/{epochs}]")
                print(f"Train Loss: {avg_train_loss:.4f}, Val Loss: {val_loss.item():.4f}")
                print(f"Val accuracy: {accuracy.item():.3f}")
                print(f"Learning rate: {scheduler.get_last_lr()[0]:.6f}")

            scheduler.step()

    def save(self, path):
        torch.save(self.model.state_dict(), path)

    def load(self, path, device="mps"):
        self.model.load_state_dict(torch.load(path, map_location=device))
        self.model.to(self.device)

    def inference(self, df):
        """Predicted probability of Disease_Label = 1 for each row of df."""
        self.model.eval()
        with torch.no_grad():
            logits = self.model(self.transform(df).to(self.device))
            return torch.sigmoid(logits).squeeze(1).cpu().numpy()
