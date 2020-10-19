""" bayesian optimization plot module """
import copy
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
import GPy
from combination import combine

def add_graph(res, n, data, next_step, values, index, columns):
    """ モデルの可視化 """
    plots = data.dropna()
    kargs = {
        'values':values, 'index':[index], 'columns':[columns],
        'aggfunc':np.mean, 'dropna':False
    }

    table = pd.pivot_table(data, **kargs)
    table = table.fillna(table.mean())

    fig = plt.figure(figsize=(5, 4))

    ax = fig.add_subplot(111)
    ax.set_xlabel(columns)
    ax.set_ylabel(index)
    ax.set_title(values)

    contour = ax.contourf(
        *np.meshgrid(*[x.values.astype(np.float32) for x in [table.columns, table.index]]), table.values,
        cmap='nipy_spectral', levels=100, alpha=0.9
    )
    fig.colorbar(contour)
    contour.set_clim(vmin=data[values].min(), vmax=data[values].max())

    plt.scatter(plots[columns], plots[index], s=5, c='green')
    plt.scatter(next_step[columns], next_step[index], s=5, c='blue')

    fig.savefig(f'{res}/graphs/{values}/{values}_{str(n)}.png')
    plt.close()


def scaling(_df, _item):
    """ scaling """
    scaler = MinMaxScaler(feature_range=(0, 1)).fit(np.array([_item['min'], _item['max']]).astype(float).reshape(-1, 1))  # scaling
    _df[_item['name']] = scaler.transform(_df[_item['name']].values.reshape(-1, 1))                     # scalingした探索点
    return _df[_item['name']]

def search(res, params_csv, steps_csv):
    """ search next point """
    params_df = pd.read_csv(params_csv, header=0)
    steps_df = pd.read_csv(steps_csv, header=0)

    ### 探索項目の取得 ###
    params = list(params_df.T.to_dict().values())
    keys = [d.get('name') for d in params]
    dim = len(keys)
    spaces = {}
    for item in params:
        item['div'] = int((item['max']-item['min'])/item['step'])+1     # 探索分割数
        item['space'] = 1/(item['div']-1)                               # 探索間隔
        spaces[item['name']] = item['space']

    steps_list = list(steps_df.T.to_dict().values())

    init_num = 4
    for index in range(init_num, len(steps_list)):
        print(index, end=', ')

        ### 探索記録の取得 ###
        df = pd.DataFrame(steps_list[:index])
        n = len(df) # データサイズ

        ### 全探索点の列挙 ###
        combi = combine(params)
        for com in combi:
            for i, row in df.iterrows():
                if all([math.isclose(com[k], row[k], abs_tol=spaces[k]/1e3) for k in keys]):
                    com['score'] = row['score']

        rest = combi
        for i, row in df.iterrows():        # 探索済みの組み合わせを削除
            rest = [item for item in rest if not all([math.isclose(item[k], row[k], abs_tol=spaces[k]/1e3) for k in keys])]


        ### 探索項目の整形とスケーリング ###
        s_df = copy.copy(df)
        for item in params:
            item['df'] = scaling(df, item)  # scalingした探索記録
        x_train = np.stack([d.get('df') for d in params], axis=1)
        y_train = df["score"].values

        combi_scaled = copy.deepcopy(combi)
        for com in combi_scaled:
            for k in keys:
                scaler = MinMaxScaler(feature_range=(0, 1)).fit(np.array([params_df[params_df['name']==k]['min'].values[0], params_df[params_df['name']==k]['max'].values[0]]).astype(float).reshape(-1, 1))  # scaling
                com[k] = scaler.transform([[com[k]]])[0][0]                     # scalingした探索点

        ### 探索記録をガウス過程回帰してモデルを最適化 ###
        kern = GPy.kern.RBF(dim, ARD=True)
        model = GPy.models.GPRegression(X=x_train.reshape(-1, dim), Y=y_train.reshape(-1, 1), kernel=kern, normalizer=True)
        model.optimize(max_iters=1e5)

        ### 探索点xに対する予測出力yから回帰モデルをもとに獲得関数acqを最大化し、次に探索すべき点を探す ###
        for i, com in enumerate(combi):
            x_pred = np.array([np.array([combi_scaled[i][k] for k in keys])])
            com['y_mean'], com['y_var'] = [y[0][0] for y in list(model.predict(x_pred))]
            y_mean, y_var = com['y_mean'], com['y_var']
            com['acq'] = (y_mean + ((np.log(n) / n) ** 0.5 * y_var))

            for j, row in s_df.iterrows():
                if all([math.isclose(com[k], row[k], abs_tol=spaces[k]/1e3) for k in keys]):
                    com['acq'] = 0.0
        next_step = (steps_df.loc[index, [k for k in keys]]).to_dict()

        ### 探索点の可視化 ###
        data = pd.DataFrame(combi)
        _ = [add_graph(res, n, data, next_step, **kargs) for kargs in [
            {'values':'score', 'index':keys[0], 'columns':keys[1]}, # 探索点の計算値 score
            {'values':'y_mean', 'index':keys[0], 'columns':keys[1]},# 予測出力 y_mean
            {'values':'acq', 'index':keys[0], 'columns':keys[1]}    # 獲得関数 acq
        ]]

    #return next_step
