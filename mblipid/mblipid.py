import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mbread import read_annol
from mbvolc import check_volc_kegg, list_volc_sig


def get_lipids(ind, obj_names, fcthr=2, pthr=0.05, direc='up', met='neg'):

    annol = read_annol(ind, met)
    annol.set_index('Name', inplace=True)

    volc_list = list_volc_sig(obj_names,
                              fcthr=fcthr,
                              pthr=pthr,
                              direc=direc,
                              met=met)

    kegg_list = check_volc_kegg(ind, volc_list, obj_names).index.to_frame()
    kegg_list.index.name = ''

    resl = kegg_list.join(annol, on='Name', how='inner')
    lipid_df = resl[['CATEGORY', 'MAIN_CLASS', 'SUB_CLASS']]
    lipid_df.index.name = 'Name'
    return lipid_df


def get_lipid_groups(obj_names, fcthr=2, pthr=0.05, direc='up', met='neg'):

    lipid_groups = []

    for ind, _ in enumerate(obj_names, start=1):

        lipid_df = get_lipids(ind,
                              obj_names,
                              fcthr=fcthr,
                              pthr=pthr,
                              direc=direc,
                              met=met)
        lipid_group = lipid_df.groupby(by='MAIN_CLASS').count()[['SUB_CLASS']]
        lipid_groups.append(lipid_group)

    return lipid_groups


def join_lipids(obj_names, titles, fcthr=2, pthr=0.05, direc='up', met='neg'):

    lipid_groups = []

    for ind, _ in enumerate(obj_names, start=1):

        lipid_df = get_lipids(ind,
                              obj_names,
                              fcthr=fcthr,
                              pthr=pthr,
                              direc=direc,
                              met=met)
        lipid_group = lipid_df.groupby(by='MAIN_CLASS').count()[['SUB_CLASS']]
        lipid_groups.append(lipid_group)

    lipid_total = pd.concat(lipid_groups, axis=1, join='outer')
    lipid_total.columns = titles
    lipid_total.fillna(0, inplace=True)
    lipid_total.sort_index(inplace=True)
    if direc == 'up':
        lipid_total = lipid_total.assign(direction=["up"] *
                                         len(lipid_total.index))
    else:
        lipid_total = lipid_total.assign(direction=["down"] *
                                         len(lipid_total.index))
    return lipid_total


def plot_lipids_count(ax, df, titles, subtitle, baralpha=0.5):

    barWidth = 0.2
    lips = df.index.tolist()

    r1 = np.arange(len(lips))
    rs = [r1]

    for _ in range(len(titles)):
        ri = [x + barWidth for x in rs[-1]]
        rs.append(ri)

    for ind, title in enumerate(titles):
        bars = df[title].values

        ax.barh(rs[ind],
                bars,
                height=barWidth,
                alpha=baralpha,
                label=titles[ind])

    ax.set_yticks([x + barWidth / 2 - 0.005 for x in rs[1]])
    ax.set_title(subtitle)
