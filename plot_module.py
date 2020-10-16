""" bayesian optimization plot module """
import numpy as np
import pandas as pd
import copy
from sklearn.preprocessing import MinMaxScaler
import GPy
from combination import combine
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.cm import ScalarMappable

def gaussian_2dim(res, combi, _keys, step):
    """ モデルの可視化 """
    data = pd.DataFrame(combi)

    table = pd.pivot_table(data, values='score', index=[_keys[0]], columns=[_keys[1]], aggfunc=np.mean, dropna=False)
    table = table.fillna(table.mean())

    Xgrid = table.columns.values.astype(np.float32)
    Ygrid = table.index.values.astype(np.float32)
    X, Y = np.meshgrid(Xgrid, Ygrid)
    Z = table.values
    
    fig = plt.figure(figsize=(12, 4))
    ax = fig.add_subplot(121)
    ax.set_title("score")
    contour = ax.contourf(X, Y, Z)
    fig.colorbar(contour)
    contour.set_clim(vmin=0, vmax=data['score'].max())

    ax.set_xlabel(table.columns.name)
    ax.set_ylabel(table.index.name)

    plots = data.dropna()
    plt.scatter(plots[table.columns.name], plots[table.index.name], s=5, c='black')

    #plt.show()
    fig.savefig(res+'/graphs/scores/score_'+str(step)+'.png')
    plt.close()

def heatmap(res, combi, keys, step, next_step):
    """ モデルの可視化 """
    data = pd.DataFrame(combi)
    plots = data.dropna()

    table = pd.pivot_table(data, values='y_mean', index=[keys[0]], columns=[keys[1]], aggfunc=np.mean)
    #print(table)

    Xgrid = table.columns.values.astype(np.float32)
    Ygrid = table.index.values.astype(np.float32)
    X, Y = np.meshgrid(Xgrid, Ygrid)
    Z = table.values
    
    fig = plt.figure(figsize=(12, 4))
    ax = fig.add_subplot(121)
    ax.set_title("contour")
    contour = ax.contourf(X, Y, Z)
    fig.colorbar(contour)
    contour.set_clim(vmin=data['y_mean'].min(), vmax=data['y_mean'].max())

    ax.set_xlabel(table.columns.name)
    ax.set_ylabel(table.index.name)

    plt.scatter(plots[table.columns.name], plots[table.index.name], s=5, c='black')
    plt.scatter(next_step[table.columns.name], next_step[table.index.name], s=5, c='red')

    #plt.show()
    fig.savefig(res+'/graphs/heatmaps/heatmap_'+str(step)+'.png')
    plt.close()

    table = pd.pivot_table(data, values='acq', index=[keys[0]], columns=[keys[1]], aggfunc=np.mean)
    #print(table)

    Xgrid = table.columns.values.astype(np.float32)
    Ygrid = table.index.values.astype(np.float32)
    X, Y = np.meshgrid(Xgrid, Ygrid)
    Z = table.values
    
    fig = plt.figure(figsize=(12, 4))
    ax = fig.add_subplot(121)
    ax.set_title("contour")
    contour = ax.contourf(X, Y, Z)
    fig.colorbar(contour)
    contour.set_clim(vmin=np.min(data['acq'].values[np.nonzero(data['acq'].values)]), vmax=data['acq'].max())

    ax.set_xlabel(table.columns.name)
    ax.set_ylabel(table.index.name)

    plt.scatter(plots[table.columns.name], plots[table.index.name], s=5, c='black')
    plt.scatter(next_step[table.columns.name], next_step[table.index.name], s=5, c='red')

    #plt.show()
#    plt.figure()
#    sns.heatmap(table)
    fig.savefig(res+'/graphs/acqs/acq_'+str(step)+'.png')
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

    params_list = params_df.reset_index().T.reset_index().T.values.tolist()
    params = []
    for i in range(len(params_list)-1):
        item = {}
        for j in range(1, len(params_list[0])):
            item[params_list[0][j]] = params_list[i+1][j]
        params.append(item)

    spaces = {}
    for item in params:
        item['div'] = int((item['max']-item['min'])/item['step'])+1     # 探索分割数
        item['space'] = 1/(item['div']-1)                               # 探索間隔
        spaces[item['name']] = item['space']

    steps_list = steps_df.reset_index().T.reset_index().T.values.tolist()

    init_num = 5
    for index in range(init_num, len(steps_list)):
        print(index)

        steps = []
        for i in range(index):
            item = []
            for j in range(1, len(steps_list[0])):
                item.append(steps_list[i+1][j])
            steps.append(item)

        ### 探索項目の取得 ###
        keys = [d.get('name') for d in params]
        dim = len(keys)
        col = np.append(keys, 'score')
        df = pd.DataFrame(np.arange(len(col)).reshape(1, len(col)), columns=col)

        ### 探索記録の取得 ###
        for i, step in enumerate(steps):
            df.loc[i] = step
        n = len(df) # データサイズ

        ### 全探索点の列挙 ###
        combi = combine(params)
        combi = [item for item in combi if not any([float(item[k]) == 0.0 for k in keys])]  # 値が0の探索点を削除
        rest = combi
        for i, row in df.iterrows():                                                        # 探索済みの組み合わせを削除
            rest = [item for item in rest if not all([abs(round(item[k] - row[k], 5)) < float(spaces[k]) for k in keys])]

        for com in combi:
            for i, row in df.iterrows():
                if all([abs(round(com[k] - row[k], 5)) < float(spaces[k]) for k in keys]):
                    com['score'] = row['score']

        gaussian_2dim(res, combi, keys, n)      # 探索点の可視化

        ### 探索項目の整形とスケーリング ###
        s_df = copy.copy(df)
        for item in params:
            item['df'] = scaling(df, item)                                  # scalingした探索記録
        x_train = np.stack([d.get('df') for d in params], axis=1)
        y_train = df["score"].values

        combi_scaled = copy.deepcopy(combi)
        for com in combi_scaled:
            for k in keys:
                scaler = MinMaxScaler(feature_range=(0, 1)).fit(np.array([params_df[params_df['name']==k]['min'].values[0], params_df[params_df['name']==k]['max'].values[0]]).astype(float).reshape(-1, 1))  # scaling
                com[k] = scaler.transform([[com[k]]])[0][0]                     # scalingした探索点
        rest_scaled = copy.deepcopy(combi_scaled)

        ### 探索記録をガウス過程回帰してモデルを最適化 ###
        kern = GPy.kern.RBF(dim, ARD=True)
        model = GPy.models.GPRegression(X=x_train.reshape(-1, dim), Y=y_train.reshape(-1, 1), kernel=kern, normalizer=True)
        model.optimize(max_iters=1e5)

        ### 探索点xに対する予測出力yから回帰モデルをもとに獲得関数acqを最大化し、次に探索すべき点を探す ###
        next_step = {}
        if len(rest) > 0:
            max_acq = {'acq':0.0}
            for k in keys:
                max_acq[k] = rest[0][k]
            for i, com in enumerate(rest):
                x_pred = np.array([np.array([rest_scaled[i][k] for k in keys])])
                y_mean, y_var = model.predict(x_pred)
                com['acq'] = (y_mean + ((np.log(n) / n) ** 0.5 * y_var)).tolist()[0][0] #need for optimization
                if max_acq['acq'] < com['acq']:
                    max_acq['acq'] = com['acq']
                    for k in keys:
                        max_acq[k] = com[k]

            #for item in params:
            #    item['next'] = max_acq[item['name']]
            #    next_step[item['name']] = item['next']
            next_step = (steps_df.loc[index, [k for k in keys]]).to_dict()

            ### heatmap ###
            for i, com in enumerate(combi):
                x_pred = np.array([np.array([combi_scaled[i][k] for k in keys])])
                y_mean, y_var = model.predict(x_pred)
                com['y_mean'] = y_mean[0][0]
                com['y_var'] = y_var[0][0]
                com['acq'] = (y_mean + ((np.log(n) / n) ** 0.5 * y_var)).tolist()[0][0] #need to optimize

                for j, row in s_df.iterrows():
                    if all([abs(round(com[k] - row[k], 5)) < float(spaces[k]) for k in keys]):
                        com['acq'] = 0.0

            heatmap(res, combi, keys, len(x_train), next_step)

    return next_step
