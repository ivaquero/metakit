import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mbread import read_anno, read_volcano


def add_volc_distance(volc_df, fcthr=2, pthr=0.05):

    volc_df['log2(FC)'] = np.log2(volc_df['FC'])
    volc_df['-log10(p)'] = -np.log10(volc_df['raw.pval'])

    dist1 = ((volc_df['log2(FC)'] - fcthr)**2 +
             (volc_df['-log10(p)'] - pthr)**2)**0.5
    dist2 = ((volc_df['log2(FC)'] + fcthr)**2 +
             (volc_df['-log10(p)'] - pthr)**2)**0.5

    volc_df['distance'] = np.where(volc_df['FC'] >= fcthr, dist1, dist2)
    volc_df['significant'] = np.where(volc_df['raw.pval'] <= pthr, "True",
                                      "False")

    return volc_df.sort_values(by='distance', ascending=False)


def count_vol(obj_names, fcthr=2, pthr=0.05, met='neg'):

    nms = pd.DataFrame()

    for ind, obj in enumerate(obj_names, start=1):
        volc_df = read_volcano(ind, met)
        vol = add_volc_distance(volc_df, fcthr=fcthr,
                                pthr=pthr).query('raw.pval < pthr')
        volc_up = vol.query('FC>4').shape[0]
        volc_down = vol.query('FC<0.25').shape[0]
        nms.index = ['up regulated', 'down regulated', 'total']
        nms[obj] = [volc_up, volc_down, volc_up + volc_down]

    return nms


def list_volc_sig(obj_names, fcthr=2, pthr=0.05, direc='up', met='neg'):

    checks = []

    for ind, _ in enumerate(obj_names, start=1):
        volc_df = read_volcano(ind, met)
        vf = add_volc_distance(volc_df, fcthr=fcthr, pthr=pthr)
        vf = vf[vf['raw.pval'] < pthr]

        if direc == 'up':
            vf = vf[vf['FC'] > fcthr**2]
        elif direc == 'down':
            vf = vf[vf['FC'] < 1 / fcthr**2]

        vf = vf.reset_index()['Name']
        checks.append(vf)

    volc_list = pd.concat(checks, axis=1)
    volc_list.columns = obj_names

    return volc_list


def list_volc_sig2(obj_names, fcthr=2, pthr=0.05, met='neg'):

    volc_list = list_volc_sig(obj_names,
                              fcthr=fcthr,
                              pthr=pthr,
                              direc='up',
                              met=met)
    volc_list_dn = list_volc_sig(obj_names,
                                 fcthr=fcthr,
                                 pthr=pthr,
                                 direc='down',
                                 met=met)
    return [volc_list, volc_list_dn]


def check_volc_kegg(ind, volc_list, obj_names, reduced=True, met='neg'):

    anno_df = read_anno(ind, met)
    ch = volc_list[[obj_names[ind - 1]]]
    ch.columns.values[0] = 'Name'
    ch.set_index('Name', inplace=True)

    res = ch.join(anno_df, how='inner')
    res = res[res['Kegg_map'].str.contains("nan") == False].replace('map\d+',
                                                                    '',
                                                                    regex=True)

    if reduced:
        res = res['Kegg_map'].str.replace(
            'Metabolic pathways;',
            '').str.replace('Metabolic pathways', '').str.replace(
                'Microbial metabolism',
                '').str.replace('in diverse environments', '').str.replace(
                    'Biosynthesis of secondary metabolites',
                    '').str.replace('     ;', '').str.replace('   ;',
                                                              '').str.strip()

    return pd.DataFrame(res)


def check_volc_kegg2(ind,
                     obj_names,
                     fcthr=2,
                     pthr=0.05,
                     reduced=True,
                     met='neg'):

    volc_list, volc_list_dn = list_volc_sig2(obj_names,
                                             fcthr=fcthr,
                                             pthr=pthr,
                                             met=met)
    kegg_list = check_volc_kegg(ind,
                                volc_list,
                                obj_names,
                                reduced=reduced,
                                met=met)
    kegg_list_dn = check_volc_kegg(ind,
                                   volc_list_dn,
                                   obj_names,
                                   reduced=reduced,
                                   met=met)
    return [kegg_list, kegg_list_dn]


def list_volc_kegg(volc_list, obj_names, reduced=True, met='neg'):

    res_names = pd.DataFrame()

    for ind, _ in enumerate(obj_names, start=1):
        anno_df = read_anno(ind, met)
        resup_s = check_volc_kegg(ind, volc_list, obj_names, reduced,
                                  met).reset_index()['Name']
        res_names = pd.concat([res_names, resup_s], axis=1)

    res_names.columns = obj_names
    return res_names


def tabularize_volc_kegg(volc_list, obj_names, reduced=True, met='neg'):

    res_table = pd.DataFrame()

    for ind, obj in enumerate(obj_names, start=1):
        nnn = pd.DataFrame({
            'Name': [obj_names[ind - 1]],
            'Kegg_map': '-' * 10
        })
        res_table = res_table.append(nnn)
        ress = check_volc_kegg(ind, volc_list, obj_names, reduced,
                               met).reset_index()
        res_table = res_table.append(ress)

    return res_table.dropna()


def sec_vol(volc_df, fcthr=2, pthr=0.05):

    volc_df['-log10(p)'] = -np.log10(volc_df['raw.pval'])

    v1 = volc_df[(volc_df['log2(FC)'] > fcthr) & (volc_df['raw.pval'] < pthr)]
    v2 = volc_df[(volc_df['log2(FC)'] < -fcthr) & (volc_df['raw.pval'] < pthr)]
    v3 = volc_df[(volc_df['log2(FC)'] <= fcthr)
                 & (volc_df['log2(FC)'] >= -fcthr)]
    v4 = volc_df[(volc_df['raw.pval'] > pthr)]

    return [v1, v2, v3, v4]


def plot_volcano(ax, v, title, fcthr=2, pthr=0.05, n_labeled=0):

    seced_vol = sec_vol(v, fcthr=fcthr, pthr=pthr)
    colors = ['red', 'green', 'grey', 'grey']

    for vv, color in zip(seced_vol, colors):
        ax.scatter(vv['log2(FC)'], vv['-log10(p)'], c=color, s=20)

    if n_labeled > 0:
        for vc in seced_vol[0:2]:
            vc = vc.copy()

            for vi in range(n_labeled):
                ax.text(vc.iloc[vi, :]['log2(FC)'],
                        vc.iloc[vi, :]['-log10(p)'] + 0.25,
                        vc.index[vi],
                        bbox=dict(
                            boxstyle="round",
                            ec=(1., 0.5, 0.5),
                            fc=(1., 0.9, 0.5),
                        ))

    ax.vlines([-2, 2], -1, 5, colors='black', linestyles='-.', linewidth=.5)
    ax.hlines([0.05], -20, 20, colors='black', linestyles='-.', linewidth=.5)
    ax.set_title(title, fontsize='large')
    return ax
