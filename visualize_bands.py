#!/usr/bin/env python3
"""Visualize the 12-band S-field color scheme."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# 12-band color scheme from the MATLAB code
band_cmap = np.array([
    [1.00, 0.00, 0.00],   # 1: bright red       S<1.0
    [1.00, 0.45, 0.00],   # 2: orange-red       1.0<=S<1.3
    [1.00, 0.20, 0.65],   # 3: hot pink         1.3<=S<1.5
    [0.20, 1.00, 0.30],   # 4: bright green     1.5<=S<1.8
    [0.13, 0.57, 0.69],   # 5: dark teal        1.8<=S<=2.0
    [0.99, 0.91, 0.14],   # 6: yellow           2.0<S<=2.2
    [0.00, 0.00, 0.70],   # 7: dark blue        2.2<S<=2.5
    [0.00, 0.45, 0.00],   # 8: dark green       2.5<S<=2.7
    [0.00, 0.00, 0.00],   # 9: black            2.7<S<=3.0
    [0.40, 0.00, 0.55],   # 10: dark purple     3.0<S<=3.5
    [0.80, 0.93, 0.93],   # 11: eggshell blue   3.5<S<=4.0
    [0.72, 0.53, 0.04],   # 12: deep goldenrod  S>4.0
])

labels = [
    'S < 1.0\n(red)',
    '1.0 ≤ S < 1.3\n(orange-red)',
    '1.3 ≤ S < 1.5\n(hot pink)',
    '1.5 ≤ S < 1.8\n(green)',
    '1.8 ≤ S ≤ 2.0\n(teal)',
    '2.0 < S ≤ 2.2\n(yellow)',
    '2.2 < S ≤ 2.5\n(dark blue)',
    '2.5 < S ≤ 2.7\n(dark green)',
    '2.7 < S ≤ 3.0\n(black)',
    '3.0 < S ≤ 3.5\n(purple)',
    '3.5 < S ≤ 4.0\n(eggshell blue)',
    'S > 4.0\n(goldenrod)',
]

# Create figure with color swatches
fig, axes = plt.subplots(3, 4, figsize=(14, 9))
fig.suptitle('12-Band S-field Color Scheme', fontsize=16, fontweight='bold')

for idx, (ax, color, label) in enumerate(zip(axes.flat, band_cmap, labels)):
    ax.add_patch(mpatches.Rectangle((0, 0), 1, 1, facecolor=color, edgecolor='black', linewidth=2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f'Band {idx+1}', fontweight='bold')
    ax.text(0.5, -0.35, label, ha='center', va='top', transform=ax.transAxes, fontsize=10)

plt.tight_layout()
plt.savefig('/home/user/Kick-drift-LSVD/band_swatches.png', dpi=150, bbox_inches='tight')
print("Saved color swatches to band_swatches.png")

# Create a second figure showing the bands in context (simulated S-field)
fig2, axes2 = plt.subplots(1, 2, figsize=(14, 6))

# Left: Raw band image
band_image = np.zeros((200, 400, 3))
boundaries = [1.0, 1.3, 1.5, 1.8, 2.0, 2.2, 2.5, 2.7, 3.0, 3.5, 4.0]
s_values = np.linspace(0.5, 4.5, 400)

for i in range(200):
    band = np.digitize(s_values, boundaries) + 1
    band = np.minimum(band, 12)  # Cap at 12
    for j in range(400):
        band_image[i, j, :] = band_cmap[band[j] - 1]

axes2[0].imshow(band_image)
axes2[0].set_title('Simulated S-field (gradient)', fontweight='bold')
axes2[0].set_xlabel('S value (left to right: 0.5 → 4.5)')
axes2[0].set_ylabel('Spatial position')
axes2[0].set_xticks(np.linspace(0, 400, 10))
axes2[0].set_xticklabels([f'{v:.1f}' for v in np.linspace(0.5, 4.5, 10)])

# Right: Simulated high-stretch zones (bands 7-12)
stretch_zones = np.zeros((200, 400, 3))
s_high = np.linspace(2.2, 4.5, 400)
for i in range(200):
    zone = np.digitize(s_high, boundaries[6:]) + 7  # Start from band 7
    zone = np.minimum(zone, 12)
    for j in range(400):
        stretch_zones[i, j, :] = band_cmap[zone[j] - 1]

axes2[1].imshow(stretch_zones)
axes2[1].set_title('High-stretch zones (Bands 7-12)', fontweight='bold')
axes2[1].set_xlabel('S value (left to right: 2.2 → 4.5)')
axes2[1].set_ylabel('Spatial position')
axes2[1].set_xticks(np.linspace(0, 400, 5))
axes2[1].set_xticklabels([f'{v:.1f}' for v in np.linspace(2.2, 4.5, 5)])

plt.tight_layout()
plt.savefig('/home/user/Kick-drift-LSVD/band_simulation.png', dpi=150, bbox_inches='tight')
print("Saved band simulation to band_simulation.png")

print("\n=== Color Analysis ===")
print("\nPotentially problematic band pairs (low contrast):")
print("1. Bands 8 (dark green [0.00 0.45 0.00]) vs 9 (black [0.00 0.00 0.00])")
print("   - Very dark, may be hard to distinguish at small scales")
print("\n2. Bands 11 (brown [0.50 0.35 0.10]) vs 12 (goldenrod [0.72 0.53 0.04])")
print("   - Both in warm/brown family, but goldenrod is brighter")
print("\n3. Band 7 (dark blue [0.00 0.00 0.70]) next to Band 8 (dark green)")
print("   - Both dark, may lose clarity if regions are thin")
print("\nWell-contrasted pairs (Bands 1-6):")
print("- Bright red, orange-red, hot pink, bright green, dark teal, yellow")
print("- These have excellent color separation and should be clearly visible")

plt.show()
