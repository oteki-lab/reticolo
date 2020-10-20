""" bayesian optimization plot module """
import os
import re
from PIL import Image
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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


def draw_graph(res, keys):
    """ draw graphs """
    steps_df = pd.read_csv(f'{res}/hp_steps.csv', header=0)

    bo_data_path = res+'\\bo_data'
    bo_data = sorted(os.listdir(bo_data_path), key=lambda s: int(re.search(r'\d+', s).group()))
    bo_data = [f'{bo_data_path}/{d}' for d in bo_data]

    for path in bo_data:
        n = int(re.findall(r'bo_data_(\d+)', path)[0])
        data = pd.read_pickle(path)
        next_step = (steps_df.loc[n, [k for k in keys]]).to_dict()
        _ = [add_graph(res, n, data, next_step, **kargs) for kargs in [
            {'values':'score', 'index':keys[0], 'columns':keys[1]}, # 探索点の計算値 score
            {'values':'y_mean', 'index':keys[0], 'columns':keys[1]},# 予測出力 y_mean
            {'values':'acq', 'index':keys[0], 'columns':keys[1]}    # 獲得関数 acq
        ]]

    image_dir = res+"/graphs/"
    os.makedirs(image_dir + '/gifs', exist_ok=True)
    image_path = ["score", "y_mean", "acq"]

    for path in image_path:
        X = []
        img_paths = os.listdir(image_dir + '/' + path)
        img_paths = sorted(img_paths, key=lambda s: int(re.search(r'\d+', s).group()))

        for img_name in img_paths:
            X.append(Image.open(image_dir + '/' + path + '/' + img_name))

        X[0].save(image_dir + '/gifs/' + path + '.gif', save_all=True, append_images=X[1:], optimize=False, duration=1000, loop=0)

if __name__ == "__main__":
    RES = 'Results\\20201020180847496'
    keys_list = ['h2', 'h']
    draw_graph(RES, keys_list)
