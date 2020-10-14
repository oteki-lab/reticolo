""" bayesian optimization """
import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
import GPy
from combination import combine

def scaling(_df, _item):
    """ scaling """
    scaler = MinMaxScaler(feature_range=(0, 1)).fit(np.array([_item['min'], _item['max']]).astype(float).reshape(-1, 1))  # scaling
    _df[_item['name']] = scaler.transform(_df[_item['name']].values.reshape(-1, 1))                     # scalingした探索点
    return _df[_item['name']]

def search(params_csv, steps_csv):
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
    steps = []
    for i in range(len(steps_list)-1):
        item = []
        for j in range(1, len(steps_list[0])):
            item.append(steps_list[i+1][j])
        steps.append(item)

    ### 探索項目の取得 ###
    keys = [d.get('name') for d in params]
    col = np.append(keys, 'score')
    df = pd.DataFrame(np.arange(len(col)).reshape(1, len(col)), columns=col)

    ### 探索記録の取得 ###
    for i, step in enumerate(steps):
        df.loc[i] = step
    n = len(df) # データサイズ

    ### 全探索点の列挙 ###
    combi = combine(params)
    combi = [item for item in combi if not any([float(item[k]) == 0.0 for k in keys])]  # 値が0の探索点を削除
    for i, row in df.iterrows():                                                        # 探索済みの組み合わせを削除
        combi = [item for item in combi if not all([abs(float(item[k]) - float(row[k])) < spaces[k] for k in keys])]
    """combi = [
        [
            item for item in combi if (not any([float(item[k]) == 0.0 for k in keys]) and (not all([float(item[k]) == float(row[k]) for k in keys])))
        ] for row in df.iterrows()
    ]    # 探索済みの組み合わせを削除"""

    ### 探索項目の整形とスケーリング ###
    for item in params:
        item['df'] = scaling(df, item)                                  # scalingした探索記録
    x_train = np.stack([d.get('df') for d in params], axis=1)
    y_train = df["score"].values

    ### 探索記録をガウス過程回帰してモデルを最適化 ###
    kern = GPy.kern.RBF(2, ARD=True)
    model = GPy.models.GPRegression(X=x_train.reshape(-1, 2), Y=y_train.reshape(-1, 1), kernel=kern, normalizer=True)
    model.optimize()

    ### 探索点xに対する予測出力yから回帰モデルをもとに獲得関数acqを最大化し、次に探索すべき点を探す ###
    for com in combi:
        x_pred = np.array([np.array([com[k] for k in keys]).T])
        y_mean, y_var = model.predict(x_pred)
        com['acq'] = (y_mean + ((np.log(n) / n) ** 0.5 * y_var)).tolist()[0][0]
    max_acq = [com for com in combi if com['acq'] == max(com['acq'] for com in combi)][0]

    next_step = {}
    for item in params:
        item['next'] = max_acq[item['name']]
        next_step[item['name']] = item['next']

    #from plot_module import gaussian_2dim, gaussian_1dim
    #gaussian_2dim(model, keys, params, n, y_train)      # モデルの可視化
    #x = np.linspace(0, 1.0, params[0]['div'])
    #gaussian_1dim(model, keys, params, n, x)            # 獲得関数の可視化

    return next_step
