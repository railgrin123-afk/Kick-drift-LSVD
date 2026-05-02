#!/usr/bin/env python3
"""Visualize proposed dot colors for each S-band, shown against their band backgrounds."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Band background colors (current 12-band scheme)
band_cmap = np.array([
    [1.00, 0.00, 0.00],   # 1: bright red
    [1.00, 0.45, 0.00],   # 2: orange-red
    [1.00, 0.20, 0.65],   # 3: hot pink
    [0.20, 1.00, 0.30],   # 4: bright green
    [0.13, 0.57, 0.69],   # 5: dark teal
    [0.99, 0.91, 0.14],   # 6: yellow
    [0.00, 0.00, 0.70],   # 7: dark blue
    [0.00, 0.45, 0.00],   # 8: dark green
    [0.00, 0.00, 0.00],   # 9: black
    [0.40, 0.00, 0.55],   # 10: dark purple
    [0.80, 0.93, 0.93],   # 11: eggshell blue
    [0.72, 0.53, 0.04],   # 12: deep goldenrod
])

# Proposed dot colors -- each chosen to CONTRAST with its band background
dot_cmap = np.array([
    [0.00, 1.00, 1.00],   # 1: cyan         on red bg
    [0.00, 0.00, 0.60],   # 2: navy         on orange-red bg
    [0.40, 1.00, 0.00],   # 3: lime         on hot pink bg
    [1.00, 0.00, 0.80],   # 4: magenta      on green bg
    [1.00, 0.50, 0.00],   # 5: bright orange on teal bg
    [0.50, 0.00, 0.70],   # 6: deep purple  on yellow bg
    [1.00, 1.00, 0.20],   # 7: neon yellow  on dark blue bg
    [1.00, 0.40, 0.80],   # 8: hot pink     on dark green bg
    [1.00, 1.00, 1.00],   # 9: white        on black bg
    [0.70, 1.00, 0.30],   # 10: yellow-green on purple bg
    [0.80, 0.00, 0.00],   # 11: deep red    on eggshell blue bg
    [0.10, 0.30, 1.00],   # 12: electric blue on goldenrod bg
])

s_ranges = [
    'S < 1.0', '1.0 ≤ S < 1.3', '1.3 ≤ S < 1.5',
    '1.5 ≤ S < 1.8', '1.8 ≤ S ≤ 2.0', '2.0 < S ≤ 2.2',
    '2.2 < S ≤ 2.5', '2.5 < S ≤ 2.7', '2.7 < S ≤ 3.0',
    '3.0 < S ≤ 3.5', '3.5 < S ≤ 4.0', 'S > 4.0',
]

dot_names = [
    'cyan', 'navy', 'lime', 'magenta',
    'bright orange', 'deep purple', 'neon yellow', 'hot pink',
    'white', 'yellow-green', 'deep red', 'electric blue',
]

# Create figure: each row shows the band bg + dot color overlaid
fig, axes = plt.subplots(3, 4, figsize=(16, 11))
fig.suptitle('Proposed Dot Colors for 12-Band S-field\n(dot shown on top of its band background)',
             fontsize=15, fontweight='bold')

for idx, (ax, bgcolor, dotcolor, srange, dname) in enumerate(
        zip(axes.flat, band_cmap, dot_cmap, s_ranges, dot_names)):
    # paint background = band color
    ax.add_patch(mpatches.Rectangle((0, 0), 1, 1, facecolor=bgcolor, edgecolor='black', linewidth=1.5))
    # scatter sample dots in the cell
    rng = np.random.default_rng(42 + idx)
    n_dots = 25
    xs = rng.uniform(0.08, 0.92, n_dots)
    ys = rng.uniform(0.08, 0.92, n_dots)
    ax.scatter(xs, ys, s=80, c=[dotcolor], edgecolors='black', linewidths=0.4, zorder=5)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f'Band {idx+1}:  {srange}', fontweight='bold', fontsize=10)
    ax.text(0.5, -0.12, f'dot: {dname}\nRGB={tuple(round(v,2) for v in dotcolor)}',
            ha='center', va='top', transform=ax.transAxes, fontsize=9)

plt.tight_layout()
plt.savefig('/home/user/Kick-drift-LSVD/dot_color_proposal.png', dpi=150, bbox_inches='tight')
print("Saved dot color proposal to dot_color_proposal.png")

# Second figure: just the 12 dot colors as swatches (so you can see them isolated)
fig2, axes2 = plt.subplots(3, 4, figsize=(14, 9))
fig2.suptitle('Proposed Dot Colors -- Isolated Swatches', fontsize=15, fontweight='bold')

for idx, (ax, dotcolor, srange, dname) in enumerate(
        zip(axes2.flat, dot_cmap, s_ranges, dot_names)):
    ax.add_patch(mpatches.Rectangle((0, 0), 1, 1, facecolor=dotcolor, edgecolor='black', linewidth=2))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f'Band {idx+1}: {dname}', fontweight='bold', fontsize=11)
    ax.text(0.5, -0.18, f'{srange}\nRGB={tuple(round(v,2) for v in dotcolor)}',
            ha='center', va='top', transform=ax.transAxes, fontsize=10)

plt.tight_layout()
plt.savefig('/home/user/Kick-drift-LSVD/dot_color_swatches.png', dpi=150, bbox_inches='tight')
print("Saved dot color swatches to dot_color_swatches.png")
