import os
import re
from pathlib import Path

import fire
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

names = [
    "Han Chinese (Bejing, China)", "Japanese (Tokyo, Japan)", "Southern Han Chinese",
    "Chinese Dai (Xishuangbanna, China)", "Kinh (Ho Chi Minh, Vietnam)", "NW-Europeans (Utah)", "Toscani (Italia)", "Finnish (Finland)",
    "British (England and Scotland)", "Iberian (Spain)", "Yoruba (Ibadan, Nigeria)", "Luhya (Webuye, Kenya)", "Gambian (Western Divisions, Gambia)",
    "Mende (Sierra Leone)", "Esan (Nigeria)", "African Americans (SW USA)", "African Caribbeans (Barbados)", "Mexican (Los Angeles, USA)",
    "Puerto Ricans (Puerto Rico)", "Colombians (Medellin, Colombia)", "Peruvians (Lima, Peru)", "Gujarati Indian (Houston, TX)",
    "Punjabi (Lahore, Pakistan)", "Bengali (Bangladesh)", "Sri Lankan Tamil (UK)", "Indian Telugu (UK)"
]

groups = '''
CHB
JPT
CHS
CDX
KHV
CEU
TSI
FIN
GBR
IBS
YRI
LWK
GWD
MSL
ESN
ASW
ACB
MXL
PUR
CLM
PEL
GIH
PJL
BEB
STU
ITU
'''.strip().splitlines()

def plot_generator_from_ancestry(eigenvec_path, sample_id, output_dir):
    # sample = os.environ['SAMPLE']
    # id = os.environ['ID']
    # dir = os.environ['DIR']

    df = pd.read_csv(eigenvec_path,header=None,sep=' ')
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True,exists_ok=True)

    # Define the color palette
    my_palette = sns.color_palette("Paired", n_colors=26)

    colors = sns.color_palette("Paired", n_colors=26)
    group2color = dict(zip(groups, colors))
    group2color['sample'] = 'black'
    short2full = dict(zip(groups, names))
    short2full['sample'] = sample_id
    full2color = dict(zip(names, colors))
    full2color[sample_id] = 'black'

    df['group'] = df[0].apply(lambda x: x.split('_')[0])
    df['group'] = df['group'].mask(df[0] == sample_id, 'sample')
    df['color'] = df['group'].map(group2color)
    df['full_name'] = df['group'].map(short2full)
    df['x'] = df[2] + 0.5 * df[3] - 0.5 * df[4]
    df['y'] = np.sin(np.pi / 3) * df[3] + np.sin(np.pi / 3) * df[4]

    sample = df[df['group'] == 'sample'].iloc[0]

    loc_x = sample['x']
    loc_y = sample['y']
    top_x = loc_x + ((df['x'].max() - df['x'].min()) / 30)
    bottom_x = loc_x - ((df['x'].max() - df['x'].min()) / 30)
    top_y = loc_y + ((df['y'].max() - df['y'].min()) / 30)
    bottom_y = loc_y - ((df['y'].max() - df['y'].min()) / 30)

    df_zoom = df[(df['x'] < top_x) & (df['x'] > bottom_x) & (df['y'] < top_y) & (df['y'] > bottom_y)]
    print(df_zoom.shape)
    df_zoom.head(2)

    fig = plt.figure(figsize=(12, 10))
    ax_main = fig.add_axes([0, 0, 1, 1])
    df_without_sample = df[df['group'] != 'sample']
    df_sample = df[df['group'] == 'sample']
    title = f'Ancestry {sample_id}'

    # main
    for name in names:
        df_sub = df[df['full_name'] == name]
        color = full2color.get(name, 'black')
        ax_main.scatter(df_sub['x'], df_sub['y'], s=20, marker='o', label=name, color=color)
    ax_main.scatter(df_sample['x'], df_sample['y'], s=120, marker='*', label=sample_id, color='black')
    ax_main.set_xticks([])
    ax_main.set_yticks([])
    ax_main.legend(loc='lower left', ncol=2)
    ax_main.patch.set_edgecolor('black')
    ax_main.patch.set_linewidth(1)
    ax_main.set_title(title, fontsize=20)

    # rect
    rect = patches.Rectangle((bottom_x, bottom_y), top_x - bottom_x, top_y - bottom_y, linewidth=3, edgecolor='r',
                             facecolor='none')
    ax_main.add_patch(rect)

    # zoom
    ax_zoom = fig.add_axes([0, 0.7, 0.3, 0.3])
    for name in names:
        df_sub = df_zoom[df_zoom['full_name'] == name]
        color = full2color.get(name, 'black')
        ax_zoom.scatter(df_sub['x'], df_sub['y'], s=40, marker='o', label=name, color=color)
    ax_zoom.scatter(df_sample['x'], df_sample['y'], s=240, marker='*', label=sample_id, color='black')
    ax_zoom.set_xticks([])
    ax_zoom.set_yticks([])

    fig.savefig(f"{output_dir}/AncestryPlot_py.pdf", format='pdf', dpi=150, bbox_inches='tight')

if __name__ == '__main__':
    fire.Fire(plot_generator_from_ancestry)