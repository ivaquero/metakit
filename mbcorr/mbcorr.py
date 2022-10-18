import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from mbread import read_upload_df
from mbvolc import list_volc_sig


def get_corr(ind, obj_names, fcthr=2, pthr=0.05, direc='up', met='neg'):

    volc_list = list_volc_sig(obj_names,
                              fcthr=fcthr,
                              pthr=pthr,
                              direc=direc,
                              met=met)
    obj = obj_names[ind - 1]
    obj_list = volc_list[obj].dropna().tolist()
    obj_df = read_upload_df(ind=ind, obj=obj[0])
    return obj_df[obj_list].T


def get_diff_corr(ind,
                  obj_names,
                  n=10,
                  n_sample=6,
                  fcthr=2,
                  pthr=0.05,
                  direc='up',
                  met='neg'):

    corr = get_corr(ind,
                    obj_names,
                    fcthr=fcthr,
                    pthr=pthr,
                    direc=direc,
                    met=met)
    obj_name = obj_names[ind - 1][0]

    corr1 = corr.T.iloc[:n_sample, :].values
    corr2 = corr.T.iloc[n_sample:, :].values
    Δcorr = pd.DataFrame(corr1 - corr2, columns=corr.T.columns)
    Δcorr.index = [f"Δ{obj_name}{i}" for i in range(1, n_sample + 1)]

    Δcorr = Δcorr.T
    Δcorr['total'] = Δcorr.sum(axis=1)
    if direc == 'up':
        Δcorr.sort_values('total', ascending=False, inplace=True)
    else:
        Δcorr.sort_values('total', ascending=True, inplace=True)
    Δcorr.drop('total', axis=1, inplace=True)
    return Δcorr.head(n)


def get_diff_corr_all(ind, obj_names, fcthr, pthr, met, n=10, n_sample=6):

    diff1 = get_diff_corr(
        ind,
        obj_names,
        n,
        n_sample=n_sample,
        fcthr=fcthr,
        pthr=pthr,
        met=met,
        direc='up',
    )

    diff2 = get_diff_corr(ind,
                          obj_names,
                          n,
                          n_sample=n_sample,
                          fcthr=fcthr,
                          pthr=pthr,
                          met=met,
                          direc='down')
    return diff1.append(diff2)


def get_diff_chem(diff_corr, n_reduced=40):

    name_len = [len(i) for i in diff_corr.index]
    to_reduce = diff_corr.index[np.array(name_len) > n_reduced]
    diff_corr_t = diff_corr.T.drop(to_reduce, axis=1)

    corr_chem = diff_corr_t.corr()
    corr_chem['total'] = corr_chem.sum(axis=1)
    corr_chem.sort_values('total', ascending=False, inplace=True)
    corr_chem.drop('total', axis=1, inplace=True)
    corr_chem_t = corr_chem.T
    corr_chem_t['total'] = corr_chem_t.sum(axis=1)
    corr_chem_t.sort_values('total', inplace=True)
    corr_chem_t.drop('total', axis=1, inplace=True)
    return corr_chem_t


def plot_diff_chem(ax,
                   corr,
                   cmap='RdBu_r',
                   annot=False,
                   corrange=(-1, 1),
                   fmt='.1f',
                   fontsize='xx-large',
                   cbar=False,
                   mask=True,
                   xangle=45,
                   yangle=0):

    if mask:
        mask = np.zeros_like(corr)
        mask[np.triu_indices_from(mask)] = True
        sns.heatmap(corr,
                    mask=mask,
                    square=True,
                    cmap=cmap,
                    cbar=cbar,
                    annot=annot,
                    fmt=fmt,
                    vmin=corrange[0],
                    vmax=corrange[1],
                    ax=ax)
    else:
        sns.heatmap(corr,
                    cmap=cmap,
                    cbar=cbar,
                    annot=annot,
                    fmt=fmt,
                    vmin=corrange[0],
                    vmax=corrange[1],
                    ax=ax)

    plt.setp(ax.get_xticklabels(),
             rotation=xangle,
             fontsize=fontsize,
             horizontalalignment='right')
    plt.setp(ax.get_yticklabels(),
             rotation=yangle,
             fontsize=fontsize,
             horizontalalignment='right')
