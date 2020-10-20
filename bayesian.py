""" bayesian optimization """
import math
import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
import GPy
from combination import combine

def acq(mean, var, N):
    """ Acquisition function """
    return mean + ((np.log(N) / N) ** 0.5 * var)

def scaling(com, pdf, k):
    """ scaling function """
    scaler = MinMaxScaler(feature_range=(0, 1)).fit(np.array([pdf[pdf['name'] == k]['min'].values[0], pdf[pdf['name'] == k]['max'].values[0]]).astype(float).reshape(-1, 1))
    return scaler.transform([[com]])[0][0]

def search(res, params_csv, steps_csv):
    """ predict next search point """
    ### 探索項目の取得 ###
    params_df = pd.read_csv(params_csv, header=0)
    params = list(params_df.T.to_dict().values())
    keys = [d.get('name') for d in params]
    dim = len(keys)
    spaces = {}
    for item in params:
        item['div'] = int((item['max']-item['min'])/item['step'])+1     # 探索分割数
        item['space'] = 1/(item['div']-1)                               # 探索間隔
        spaces[item['name']] = item['space']

    ### 探索記録の取得 ###
    steps_df = pd.read_csv(steps_csv, header=0)
    list_df = pd.DataFrame(list(steps_df.T.to_dict().values()))
    step = len(list_df)                                                 # データサイズ
    for item in params:                                                 # 探索記録のスケーリング
        scaler = MinMaxScaler(feature_range=(0, 1)).fit(np.array([item['min'], item['max']]).astype(float).reshape(-1, 1))
        item['list'] = scaler.transform(list_df[item['name']].values.reshape(-1, 1))

    ### 探索点の設定 ###
    combi = combine(params)                                             # 全探索点の列挙
    for com in combi:                                                   # 探索済点に探索記録のスコアをセット
        for _, row in list_df.iterrows():
            if all([math.isclose(com[k], row[k], abs_tol=spaces[k]/1e3) for k in keys]):
                com['score'] = row['score']

    ### 探索記録をガウス過程回帰してモデルを最適化 ###
    kargs = {
        'X': np.stack([d.get('list') for d in params], axis=1).reshape(-1, dim),
        'Y': list_df["score"].values.reshape(-1, 1),
        'kernel': GPy.kern.RBF(dim, ARD=True),
        'normalizer': None
    }
    model = GPy.models.GPRegression(**kargs)
    model.optimize(max_iters=1e5)

    ### 探索点xに対する予測出力yから回帰モデルをもとに獲得関数acqを最大化し、次に探索すべき点を探す ###
    for _, com in enumerate(combi):
        x_pred = np.array([np.array([scaling(com[k], params_df, k) for k in keys])])
        com['y_mean'], com['y_var'] = [y[0][0] for y in list(model.predict(x_pred))]
        com['acq'] = acq(com['y_mean'], com['y_var'], step)

        for _, row in list_df.iterrows():
            if all([math.isclose(com[k], row[k], abs_tol=spaces[k]/1e3) for k in keys]):
                com['acq'] = 0.0

    data = pd.DataFrame(combi)                                                    # all data
    data.to_pickle(f'{res}/bo_data/bo_data_{step}.pkl')                           # save data
    return {p['name']: data.loc[data['acq'].idxmax()][p['name']] for p in params} # return next step
