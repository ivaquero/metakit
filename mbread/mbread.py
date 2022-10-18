import pandas as pd


def unicodify(df):
    df['Name'] = df['Name'].str.replace('伪', 'α')
    df['Name'] = df['Name'].str.replace('¦Â', 'α')
    df['Name'] = df['Name'].str.replace('¦Á', 'α')
    df['Name'] = df['Name'].str.replace('尾', 'β')
    df['Name'] = df['Name'].str.replace('¦´', 'γ')
    df['Name'] = df['Name'].str.replace('螖', 'δ')
    df['Name'] = df['Name'].str.replace('¦¤', 'δ')
    df['Name'] = df['Name'].str.replace('未', 'Δ')
    df['Name'] = df['Name'].str.replace('¦Ä', 'Δ')
    df['Name'] = df['Name'].str.replace('卤', '±')
    df['Name'] = df['Name'].str.replace('¡À', '±')
    return df


def read_file(ind: int, met='neg', unicode=True):

    path = f'Result-{ind}/1.MetQuant/meta_intensity_{met}.xls'
    df = pd.read_table(path, sep='\t', encoding_errors='ignore', index_col=0)
    if unicode:
        df = unicodify(df)
    return df


def get_rsd(ind: int, met='neg', unicode=True):

    df = read_file(ind, met, unicode=unicode)
    qcs = [f'{met}_QC{i}' for i in range(1, 4)]
    rsd = df[qcs].apply(lambda x: x.std() / x.mean(), axis=1)

    thres = rsd[rsd < 0.3].shape[0] / rsd.shape[0]

    qcsm = df[qcs].mean(axis=1)

    return thres, qcsm


def get_values(ind, obj, met='neg', unicode=True):
    df = read_file(ind, met, unicode=unicode)
    col_1 = [f'{met}_{obj}r{i}' for i in range(1, 7)]
    col_2 = [f'{met}_{obj}{i}' for i in range(1, 7)]
    return df[col_1 + col_2]


def transform_df(ind, obj, met='neg', standard=True):
    _, qcsm = get_rsd(ind, met)
    g = get_values(ind, obj, met)
    transf_df = g.copy()
    for col in g.columns:
        transf_df[col] = g[col] / qcsm
    if standard:
        transf_df = transf_df.apply(lambda x: (x - x.mean()) / x.std())
    return transf_df


def add_label_row(ind: int,
                  obj: str,
                  g_names: list,
                  n: int,
                  met='neg',
                  standard=True,
                  unicode=True):

    df = read_file(ind, met, unicode=unicode)
    transf_df = transform_df(ind, obj, met, standard=True)

    labs = [f'{obj}-{g_names[0]}'] * n + [f'{obj}-{g_names[1]}'] * n
    transf_df.index = df.Name

    labels = pd.DataFrame([labs], columns=transf_df.columns)
    df_new = pd.concat([labels, transf_df])
    df_new.columns = df_new.columns.str.replace(met + '_', '')
    return df_new


def save_df(df_new, ind, obj, met='neg'):
    path = f'Result-{ind}/{obj}-{met}.csv'
    df_new.to_csv(path)


def read_upload_df(ind, obj, met='neg'):
    path = f'Result-{ind}/{obj}-{met}.csv'
    return pd.read_csv(path, index_col=0)


def read_anno(ind: int, met='neg'):
    path = f'Result-{ind}/2.MetAnnotation/KEGG/meta_{met}_kegg_anno.xls'
    df = pd.read_table(path, sep='\t', encoding_errors='ignore', index_col=0)
    return df.set_index('Name')


def read_volcano(ind: int, met='neg'):
    path = f'Result-{ind}/volcano-{met}.csv'
    volc_df = pd.read_csv(path, index_col=0)
    volc_df.index.name = 'Name'
    return volc_df.sort_values(by='FC', ascending=False)


def read_oplsda(ind, met='neg'):
    path = f'Result-{ind}/oplsda_vip-{met}.csv'
    df = pd.read_csv(path, index_col=0)
    df.index.name = 'Name'
    return df.reset_index()


def read_oplsda_model(ind, met='neg'):
    path = f'Result-{ind}/oplsda_model-{met}.csv'
    df = pd.read_csv(path, index_col=0)
    df.index.name = 'Name'
    return df


def read_annol(ind, met='neg'):
    path = f'Result-{ind}/2.MetAnnotation/Lipidmaps/meta_{met}_lipidmaps_anno.xls'
    return pd.read_table(path, sep='\t', encoding_errors='ignore', index_col=0)
