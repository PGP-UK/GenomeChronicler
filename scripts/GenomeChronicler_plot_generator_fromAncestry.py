import os
import re

import fire
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def _str2color(name, my_palette):
    all_colors = []
    temp = name
    if re.search("^CHB_[HG|NA]", temp):
        all_colors.append(my_palette[0])
    elif re.search("^JPT_[HG|NA]", temp):
        all_colors.append(my_palette[1])
    elif re.search("^CHS_[HG|NA]", temp):
        all_colors.append(my_palette[2])
    elif re.search("^CDX_[HG|NA]", temp):
        all_colors.append(my_palette[3])
    elif re.search("^KHV_[HG|NA]", temp):
        all_colors.append(my_palette[4])
    elif re.search("^CEU_[HG|NA]", temp):
        all_colors.append(my_palette[5])
    elif re.search("^TSI_[HG|NA]", temp):
        all_colors.append(my_palette[6])
    elif re.search("^FIN_[HG|NA]", temp):
        all_colors.append(my_palette[7])
    elif re.search("^GBR_[HG|NA]", temp):
        all_colors.append(my_palette[8])
    elif re.search("^IBS_[HG|NA]", temp):
        all_colors.append(my_palette[9])
    elif re.search("^YRI_[HG|NA]", temp):
        all_colors.append(my_palette[10])
    elif re.search("^LWK_[HG|NA]", temp):
        all_colors.append(my_palette[11])
    elif re.search("^GWD_[HG|NA]", temp):
        all_colors.append(my_palette[12])
    elif re.search("^MSL_[HG|NA]", temp):
        all_colors.append(my_palette[13])
    elif re.search("^ESN_[HG|NA]", temp):
        all_colors.append(my_palette[14])
    elif re.search("^ASW_[HG|NA]", temp):
        all_colors.append(my_palette[15])
    elif re.search("^ACB_[HG|NA]", temp):
        all_colors.append(my_palette[16])
    elif re.search("^MXL_[HG|NA]", temp):
        all_colors.append(my_palette[17])
    elif re.search("^PUR_[HG|NA]", temp):
        all_colors.append(my_palette[18])
    elif re.search("^CLM_[HG|NA]", temp):
        all_colors.append(my_palette[19])
    elif re.search("^PEL_[HG|NA]", temp):
        all_colors.append(my_palette[20])
    elif re.search("^GIH_[HG|NA]", temp):
        all_colors.append(my_palette[21])
    elif re.search("^PJL_[HG|NA]", temp):
        all_colors.append(my_palette[22])
    elif re.search("^BEB_[HG|NA]", temp):
        all_colors.append(my_palette[23])
    elif re.search("^STU_[HG|NA]", temp):
        all_colors.append(my_palette[24])
    elif re.search("^ITU_[HG|NA]", temp):
        all_colors.append(my_palette[25])
    else:
        all_colors.append("black")
    return all_colors[0]


def plot_generator_from_ancestry(sample, id, dir):
    # sample = os.environ['SAMPLE']
    # id = os.environ['ID']
    # dir = os.environ['DIR']

    data = pd.read_csv(f"{dir}/results/results_{sample}/temp/{sample}_1kGP_pruned_pca_20.eigenvec",header=None,sep=' ')

    # Define the color palette
    my_palette = sns.color_palette("Paired", n_colors=26)

    # Create a list to store the colors for each sample
    all_colors = []

    for i in range(len(data)):
        temp = data.iloc[i, 1]
        color = _str2color(temp, my_palette)
        all_colors.append(color)

    pgp = np.where(data[1].str.contains(id))[0][0]

    loc_x = data.iloc[pgp, 2] + 0.5 * data.iloc[pgp, 3] - 0.5 * data.iloc[pgp, 4]
    loc_y = np.sin(np.pi / 3) * data.iloc[pgp, 3] + np.sin(np.pi / 3) * data.iloc[pgp, 4]

    max_x = data[2] + 0.5 * data[3] - 0.5 * data[4]
    min_x = data[2] + 0.5 * data[3] - 0.5 * data[4]
    max_y = np.sin(np.pi / 3) * data[3] + np.sin(np.pi / 3) * data[4]
    min_y = np.sin(np.pi / 3) * data[3] + np.sin(np.pi / 3) * data[4]

    top_x = loc_x + ((max_x.max() - min_x.min()) / 30)
    bottom_x = loc_x - ((max_x.max() - min_x.min()) / 30)
    top_y = loc_y + ((max_y.max() - min_y.min()) / 30)
    bottom_y = loc_y - ((max_y.max() - min_y.min()) / 30)

    # data2<-which(data[, 3] + 0.5 * data[, 4] - 0.5 * data[, 5] < top_x & data[, 3] + 0.5 * data[, 4] - 0.5 * data[, 5] > bottom_x & sin(pi / 3) * data[,4] + sin(pi / 3) * data[, 5] < top_y & sin(pi / 3) * data[,4] + sin(pi / 3) * data[, 5] > bottom_y)
    # data2 = np.where((data[2] + 0.5 * data[3] - 0.5 * data[4] < top_x) & (data[2] + 0.5 * data[3] - 0.5 * data[4] > bottom_x) & (np.sin(np.pi / 3) * data[3] + np.sin(np.pi / 3) * data[4] < top_y)
    data2 = np.where((data[2] + 0.5 * data[3] - 0.5 * data[4] < top_x) & (data[2] + 0.5 * data[3] - 0.5 * data[4] > bottom_x) & (np.sin(np.pi / 3) * data[3] + np.sin(np.pi / 3) * data[4] < top_y) & (np.sin(np.pi / 3) * data[3] + np.sin(np.pi / 3) * data[4] > bottom_y))[0]
    zoom_data = data.loc[data2]

    # pgp2 = np.where([re.search(id, x) is not None for x in zoom_data[1]])[0]
    # zoom_data = np.vstack((zoom_data.loc[np.setdiff1d(np.arange(zoom_data.shape[0]), pgp2)], zoom_data.loc[pgp2], zoom_data.loc[pgp2]))
    pgp2 = zoom_data[1].str.contains(id)
    # print(zoom_data.drop(index=pgp2).shape)
    # print(pgp2.shape)
    # exit()
    zoom_data = pd.concat([zoom_data[~pgp2], zoom_data[pgp2], zoom_data[pgp2]])
    zoom_data = zoom_data.reset_index(drop=True)

    zoom_colors = []
    for i in range(len(zoom_data)):
        temp = zoom_data.iloc[i, 0]
        color = _str2color(temp, my_palette)
        zoom_colors.append(color)

    # data is a pandas dataframe, all.colors is a list of colors, pgp is a list of indices,
    # bottom_x, bottom_y, top_x, and top_y are values defining the rectangular bounds
    data_sub = data.iloc[pgp]
    plt.scatter(data[2] + 0.5 * data[3] - 0.5 * data[4], np.sin(np.pi / 3) * data[3] + np.sin(np.pi / 3) * data[4], c=all_colors, marker='o', s=30, edgecolors='none', alpha=0.5)
    plt.scatter(data_sub[2] + 0.5 * data_sub[3] - 0.5 * data_sub[4], np.sin(np.pi / 3) * data_sub[3] + np.sin(np.pi / 3) * data_sub[4], c='black', marker='*', s=100, edgecolors='none', alpha=1)
    plt.scatter(data_sub[2] + 0.5 * data_sub[3] - 0.5 * data_sub[4], np.sin(np.pi / 3) * data_sub[3] + np.sin(np.pi / 3) * data_sub[4], c='black', marker='o', s=30, edgecolors='none', alpha=1)

    # plt.legend(["Han Chinese (Bejing, China)", "Japanese (Tokyo, Japan)", "Southern Han Chinese",
    #     "Chinese Dai (Xishuangbanna, China)", "Kinh (Ho Chi Minh, Vietnam)", "NW-Europeans (Utah)", "Toscani (Italia)", "Finnish (Finland)",
    #     "British (England and Scotland)", "Iberian (Spain)", "Yoruba (Ibadan, Nigeria)", "Luhya (Webuye, Kenya)", "Gambian (Western Divisions, Gambia)",
    #     "Mende (Sierra Leone)", "Esan (Nigeria)", "African Americans (SW USA)", "African Caribbeans (Barbados)", "Mexican (Los Angeles, USA)",
    #     "Puerto Ricans (Puerto Rico)", "Colombians (Medellin, Colombia)", "Peruvians (Lima, Peru)", "Gujarati Indian (Houston, TX)", 
    #     "Punjabi (Lahore, Pakistan)", "Bengali (Bangladesh)", "Sri Lankan Tamil (UK)", "Indian Telugu (UK)", id], loc='lower left', ncol=2, markerscale=0.5, scatterpoints=1)

    plt.savefig(f"{dir}/results/results_{sample}/temp/AncestryPlot_py1.pdf",format='pdf')

    rect = plt.Rectangle((bottom_x, bottom_y), top_x - bottom_x, top_y - bottom_y, linewidth=2, edgecolor='red', facecolor='none')
    plt.gca().add_patch(rect)

    pchs = [19 if i < len(zoom_data) - 2 else 16 for i in range(len(zoom_data))]
    cexs = [0.5 if i < len(zoom_data) - 2 else 1.2 for i in range(len(zoom_data))]
    cexs = np.array(cexs)

    fig, ax = plt.subplots()
    ax.scatter(zoom_data[2] + 0.5 * zoom_data[3] - 0.5 * zoom_data[4], np.sin(np.pi / 3) * zoom_data[3] + np.sin(np.pi / 3) * zoom_data[4], c=zoom_colors, marker='o', s=30 * cexs, edgecolors='none', alpha=0.5)
    ax.set_xlim(bottom_x, top_x)
    ax.set_ylim(bottom_y, top_y)
    ax.set_xticks([])
    ax.set_yticks([])

    plt.savefig(f"{dir}/results/results_{sample}/temp/AncestryPlot_py2.pdf",format='pdf')

if __name__ == '__main__':
    fire.Fire(plot_generator_from_ancestry)