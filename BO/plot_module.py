""" bayesian optimization plot module """
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
import numpy as np

def gaussian_2dim(_model, _keys, _params, _n, _y_train):
    """ モデルの可視化 """
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)

    _model.plot_mean(ax=ax, cmap="jet")  # カーネル最適化後の予測

    ### データ点の描画(defaultでは見辛い) ###
    for xi, yi in zip(_params[0]['df'].values, _params[1]['df'].values):
        ax.plot(xi, yi, marker=".", color="k", markersize=10)
    ax.set_ylabel("normalized_lemon", fontsize=18)
    ax.set_xlabel("normalized_coke", fontsize=20)
    ax.tick_params(labelsize=20)

    ### color bar追加 ###
    #axpos = ax.get_position()
    cbar_ax = fig.add_axes([1, 0.15, 0.02, 0.8])
    norm = colors.Normalize(vmin=_y_train.min(), vmax=_y_train.max())
    mappable = ScalarMappable(cmap='jet', norm=norm)
    mappable.set_array([])
    cbar_ax.tick_params(labelsize=10)
    fig.colorbar(mappable, cax=cbar_ax)

    #fig.tight_layout()
    plt.savefig(f"2dim_gaussian_n={_n}.png")

def gaussian_1dim(_model, _keys, _params2, _n, _x):
    """ 獲得関数の可視化 """
    ### 入力を1次元固定して予測 ###
    _x_pred = np.array([_x, np.full(_params2[0]['div'], _params2[1]['next'])]).T
    _y_mean, _y_var = _model.predict(_x_pred)
    _acq = (_y_mean + ((np.log(_n) / _n) ** 0.5 * _y_var)) / 5   # 獲得関数acq

    _model.plot(fixed_inputs=[(1, _params2[1]['next'])], plot_data=False, plot_limits=[-.01, 1.01])
    plt.xlabel("normalized coke", fontsize=14)
    plt.ylabel("taste score [ g ]", fontsize=14)

    plt.plot(np.linspace(-.01, 1.01, 11), _acq, color="g")
    plt.plot(_acq.argmax() * 0.1, _acq.max(), marker=".", color="r", markersize=14)

    plt.legend(["Mean", "Acquisition", "Acq Max", "Confidence"])
    plt.savefig(f"1dim_gaussian_n={_n}_lemon={_params2[1]['next']}.png")
