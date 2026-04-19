import os
from pathlib import Path

from ase.units import kJ, mol
from ase.io import read

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import colorcet as cc
import pandas as pd

from utils import create_label, read_orca


plt.rcParams['legend.borderpad'] = 0.5


if __name__ == '__main__':
    data = np.load('FEPs.npz')

    x = data['cv']
    y = {
        'pbe': data['fep_mace_pbe'],
        'rpa': data['fep_mace_rpa'],
    }

    height = float(13) * 11 / 72  # define 'em' unit based on font size 11
    figure = plt.figure(figsize=(4, height))
    ax = figure.gca()
    ax.grid()
    ax.set_ylabel(create_label('free energy', '[kJ/mol]'))
    ax.set_xlabel(create_label('combined CV', '[-]'))

    symbol_size = float(os.environ['SCATTER_SYMBOL_SIZE'])
    line_width = float(os.environ['MARKER_EDGE_WIDTH'])
    legend = []

    labels = {
        'pbe': 'PBE-D3(BJ)',
        'rpa': 'RPA',
    }
    colors = {'pbe': 'RED', 'rpa': 'BLUE'}

    for method in ['pbe', 'rpa']:
        ax.plot(
            x,
            y[method],
            color=os.environ[colors[method]],
            marker='',
            linestyle='-',
            label=labels[method],
        )

    ax.set_ylim([-5, 140])
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    ax.legend()
    figure.savefig('reaction_profile.svg', bbox_inches='tight', transparent=True, format='svg')

    # surface plot
    data = np.load('fes_rpa.npz')
    x = data['cv1s']
    y = data['cv2s']
    z = data['fs']
    z -= np.min(z[~np.isnan(z)])

    levels = 10 * np.arange(15)
    figure = plt.figure(figsize=(2.3, height))
    ax = figure.gca()
    ax.set_ylabel(create_label('CV 2', '[-]'))
    ax.set_xlabel(create_label('CV 1', '[-]'))

    contour = ax.contourf(
        x,
        y,
        z,
        levels=levels,
        cmap=cc.cm.CET_I1,
        extend='both',
    )
    line_start = (5.7, 0.03)
    line_end = (7.3, 0.9)
    ax.plot([line_start[0], line_end[0]], [line_start[1], line_end[1]], 'k-',linewidth=0.5)

    def add_rotated_text(ax, text, pos, angle, offset):
        display_pos = ax.transData.transform(pos)
        offset_display = ax.transAxes.transform((offset, 0)) - ax.transAxes.transform((0, 0))
        display_pos = display_pos + offset_display
        data_pos = ax.transData.inverted().transform(display_pos)
        ax.text(data_pos[0], data_pos[1], text, rotation=angle + 20,
                rotation_mode='anchor', ha='center', va='center')

    dx = line_end[0] - line_start[0]
    dy = line_end[1] - line_start[1]
    angle = np.degrees(np.arctan2(dy, dx))
    midpoint = ((line_start[0] + line_end[0]) / 2, (line_start[1] + line_end[1]) / 2)

    add_rotated_text(ax, 'combined CV', midpoint, angle, -0.05)

    # mid_point = ((line_start[0] + line_end[0]) / 2, (line_start[1] + line_end[1]) / 2)
    # angle = np.degrees(np.arctan2(line_end[1] - line_start[1], line_end[0] - line_start[0]))
    # angle = np.degrees(np.arctan2(line_end[1] - line_start[1], line_end[0] - line_start[0]))
    # offset = 0.1  # Adjust this value to move text further from or closer to the line
    # dx = offset * np.sin(np.radians(angle))
    # dy = offset * np.cos(np.radians(angle))
    # ax.annotate('combined CV', xy=mid_point, xytext=(dx, dy), 
    #             textcoords='offset points', ha='center', va='center',
    #             rotation=angle, rotation_mode='anchor')

    contour_lines = ax.contour(x, y, z, levels=levels, colors='white', linewidths=0.1)

    cbar_ax = figure.add_axes([0.17, 0.75, 0.45, 0.05])  # [left, bottom, width, height]

    cbar = figure.colorbar(contour, cax=cbar_ax, orientation='horizontal')
    cbar.set_label('free energy\n[kJ/mol]', fontsize=9, labelpad=0.7)  # Reduced labelpad
    cbar.ax.tick_params(pad=0.2)
    cbar.ax.set_xticks([0, 40, 80, 120])
    cbar.ax.xaxis.set_ticks_position('top')

    figure.savefig('reaction_surface.svg', bbox_inches='tight', transparent=True, format='svg')
