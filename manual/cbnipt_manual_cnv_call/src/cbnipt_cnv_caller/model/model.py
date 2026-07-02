import torch
import torch.nn as nn
import torch.nn.functional as F


class ChromosomeAwareAutoEncoder(nn.Module):
    """
    Multi-feature bin-level convolutional autoencoder.

    Input:
        x:    [B, C, L]
              B = chromosome 개수 또는 sample-chromosome batch
              C = feature 개수
              L = bin length

        mask: [B, 1, L]
              valid bin = 1
              padding / invalid bin = 0

    Output:
        recon: [B, C, L]
        latent: [B, hidden_dim, L/2]
        attn: [B, 1, L/2]
    """

    def __init__(self, in_channels=4, hidden_dim=64):
        super().__init__()

        self.encoder = nn.Sequential(
            nn.Conv1d(in_channels, 32, kernel_size=7, padding=3),
            nn.GroupNorm(4, 32),
            nn.ReLU(),

            nn.Conv1d(32, 64, kernel_size=5, padding=2),
            nn.GroupNorm(8, 64),
            nn.ReLU(),

            nn.MaxPool1d(2),
        )

        self.bottleneck = nn.Sequential(
            nn.Conv1d(64, hidden_dim, kernel_size=3, padding=1),
            nn.GroupNorm(8, hidden_dim),
            nn.ReLU(),
        )

        self.attn_score = nn.Sequential(
            nn.Conv1d(hidden_dim, 32, kernel_size=1),
            nn.ReLU(),
            nn.Conv1d(32, 1, kernel_size=1),
            nn.Sigmoid()
        )

        self.decoder = nn.Sequential(
            nn.ConvTranspose1d(hidden_dim, 64, kernel_size=2, stride=2),
            nn.GroupNorm(8, 64),
            nn.ReLU(),

            nn.Conv1d(64, 32, kernel_size=5, padding=2),
            nn.GroupNorm(4, 32),
            nn.ReLU(),

            nn.Conv1d(32, in_channels, kernel_size=1)
        )

    def forward(self, x, mask=None):
        """
        x: [B, C, L]
        mask: [B, 1, L]
        """

        z = self.encoder(x)          # [B, hidden, L/2]
        z = self.bottleneck(z)       # [B, hidden, L/2]

        attn = self.attn_score(z)    # [B, 1, L/2]

        # attention은 latent denoising 용도로만 사용
        z_attn = z * attn

        recon = self.decoder(z_attn) # [B, C, L or L±1]

        # ConvTranspose 후 길이 차이 보정
        if recon.shape[-1] != x.shape[-1]:
            recon = F.interpolate(
                recon,
                size=x.shape[-1],
                mode="linear",
                align_corners=False
            )

        if mask is not None:
            recon = recon * mask

        return recon, z, attn