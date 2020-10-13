""" bayesian optimization """
import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
import GPy

def scaling(_df, _item):
    """ scaling """
    scaler = MinMaxScaler().fit(np.array([_item['min'], _item['max']]).astype(float).reshape(-1, 1))  # scaling
    _df[_item['name']] = scaler.transform(_df[_item['name']].values.reshape(-1, 1))                     # scalingした探索点
    return _df[_item['name']]

def search(params, combi, steps):
    """ search next point """
    ### 探索項目の取得 ###
    keys = [d.get('name') for d in params]
    col = np.append(keys, 'score')
    df = pd.DataFrame(np.arange(len(col)).reshape(1, len(col)), columns=col)

    ### 探索記録の取得 ###
    for i, step in enumerate(steps):
        df.loc[i] = step
    n = len(df) # データサイズ

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

    for item in params:
        item['next'] = max_acq[item['name']]
        print(item['name'], ":", item['next'], "(", item['next']*item['max'], "ml)")
    from plot_module import gaussian_2dim, gaussian_1dim
    gaussian_2dim(model, keys, params, n, y_train)      # モデルの可視化
    x = np.linspace(0, 1.0, params[0]['div'])
    gaussian_1dim(model, keys, params, n, x)            # 獲得関数の可視化

    return max_acq
