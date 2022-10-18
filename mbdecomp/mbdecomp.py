import matplotlib.pyplot as plt
import numpy as np

from mbread import read_oplsda, read_oplsda_model


def sort_oplsda(ind, obj_names, kegg_lists):

    om = read_oplsda(ind)
    kegg_list, kegg_list_dn = kegg_lists

    omup = om.join(kegg_list, on='Name').query('V1 > 1')
    omup.dropna(inplace=True)
    omup['direction'] = ['up' for i in range(len(omup.V1))]

    omdn = om.join(kegg_list_dn, on='Name').query('V1 > 1')
    omdn.dropna(inplace=True)
    omdn['direction'] = ['down' for i in range(len(omdn.V1))]
    omreg = omup.append(omdn)
    omreg = omreg.sort_values('V1', ascending=False).reset_index()
    return omreg


def plot_vip(ax,
             omreg,
             title,
             markerlabels,
             markercolors=('red', 'green'),
             markersize=5,
             line=True,
             linestyle='dotted',
             lineweight=0.7,
             linecolor='grey',
             legendloc='lower right',
             legendsize='x-large',
             ticksize='x-large',
             titlesize='xx-large'):

    om1 = omreg.query('direction=="up"')
    om2 = omreg.query('direction=="down"')

    y1 = omreg.index.max() - om1.index
    y2 = omreg.index.max() - om2.index

    if line:
        for yi in omreg.index:
            ax.axhline(yi,
                       xmin=0,
                       xmax=2,
                       ls=linestyle,
                       lw=lineweight,
                       color=linecolor)

    ax.scatter(om1.V1,
               y1,
               s=markersize,
               color=markercolors[0],
               label=markerlabels[0])
    ax.scatter(om2.V1,
               y2,
               s=markersize,
               color=markercolors[1],
               label=markerlabels[1])
    ax.set_yticks(omreg.index[::-1], )
    ax.set_yticklabels(omreg.Name, rotation=0, fontsize=ticksize)
    plt.setp(ax.get_xticklabels(), fontsize=ticksize)
    ax.legend(loc=legendloc, fontsize=legendsize)
    ax.set_title(title, fontsize=titlesize)

    return ax


def plot_oplsda_model(ax,
                      ind,
                      subtitle,
                      legendloc='lower right',
                      legendsize='x-large',
                      ticksize='x-large',
                      titlesize='xx-large'):

    df = read_oplsda_model(ind)
    barWidth = 0.2

    xs = np.arange(len(df.index))

    for ind, col in enumerate(df.columns):
        bars = ax.bar(xs + ind * barWidth,
                      df[col].values,
                      width=barWidth,
                      alpha=0.5,
                      label=col)
        ax.bar_label(bars)

    ax.set_xticks([x + 0.5 * barWidth for x in xs])
    ax.set_xticklabels(df.index, rotation=0, fontsize=ticksize)
    ax.set_title(subtitle, fontsize=titlesize)
    plt.setp(ax.get_yticklabels(), fontsize=ticksize)
    ax.legend(loc=legendloc, fontsize=legendsize)

    return ax
