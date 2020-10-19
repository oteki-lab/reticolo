""" bayesian optimization """
import copy
import math
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

    ### 探索項目の取得 ###
    params = list(params_df.T.to_dict().values())
    keys = [d.get('name') for d in params]
    dim = len(keys)
    spaces = {}
    for item in params:
        item['div'] = int((item['max']-item['min'])/item['step'])+1     # 探索分割数
        item['space'] = 1/(item['div']-1)                               # 探索間隔
        spaces[item['name']] = item['space']

    ### 探索記録の取得 ###
    steps_list = list(steps_df.T.to_dict().values())
    df = pd.DataFrame(steps_list)
    n = len(df) # データサイズ

    ### 全探索点の列挙 ###
    combi = combine(params)
    for i, row in df.iterrows():                                                        # 探索済みの組み合わせを削除
        combi = [item for item in combi if not all([math.isclose(item[k], row[k], abs_tol=spaces[k]/1e3) for k in keys])]

    combi_scaled = copy.deepcopy(combi)
    for com in combi_scaled:
        for k in keys:
            scaler = MinMaxScaler(feature_range=(0, 1)).fit(np.array([params_df[params_df['name']==k]['min'].values[0], params_df[params_df['name']==k]['max'].values[0]]).astype(float).reshape(-1, 1))  # scaling
            com[k] = scaler.transform([[com[k]]])[0][0]                     # scalingした探索点


    ### 探索項目の整形とスケーリング ###
    for item in params:
        item['df'] = scaling(df, item)                                  # scalingした探索記録
    x_train = np.stack([d.get('df') for d in params], axis=1)
    y_train = df["score"].values

    ### 探索記録をガウス過程回帰してモデルを最適化 ###
    kern = GPy.kern.RBF(dim, ARD=True)
    model = GPy.models.GPRegression(X=x_train.reshape(-1, dim), Y=y_train.reshape(-1, 1), kernel=kern, normalizer=None)
    model.optimize(max_iters=1e5)

    ### 探索点xに対する予測出力yから回帰モデルをもとに獲得関数acqを最大化し、次に探索すべき点を探す ###
    max_acq = {'acq':0.0}
    for k in keys:
        max_acq[k] = combi[0][k]
    for i, com in enumerate(combi):
        x_pred = np.array([np.array([combi_scaled[i][k] for k in keys])])
        com['y_mean'], com['y_var'] = [y[0][0] for y in list(model.predict(x_pred))]
        y_mean, y_var = model.predict(x_pred)
        com['acq'] = (y_mean + ((np.log(n) / n) ** 0.5 * y_var)).tolist()[0][0] #need for optimization
        if max_acq['acq'] < com['acq']:
            max_acq['acq'] = com['acq']
            print(x_pred, y_mean[0][0], ((np.log(n) / n) ** 0.5 * y_var)[0][0])
            for k in keys:
                max_acq[k] = com[k]
    #max_acq = [com for com in combi if com['acq'] == max(com['acq'] for com in combi)][0]

    next_step = {}
    for item in params:
        item['next'] = max_acq[item['name']]
        next_step[item['name']] = item['next']

    return next_step
