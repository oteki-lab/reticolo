""" bayesian optimization plot module """
import os
import re
from PIL import Image
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def add_graph(res, n, data, next_step, values, index, columns):
    """ モデルの可視化 """
    #data = data[data[values]>0.0]
    plots = data.dropna()
    kargs = {
        'values':values, 'index':[index], 'columns':[columns],
        'aggfunc':np.mean, 'dropna':False
    }

    table = pd.pivot_table(data, **kargs)
    #table = table.fillna(table.mean())

    fig = plt.figure(figsize=(5, 4))

    ax = fig.add_subplot(111)
    ax.set_xlabel(columns)
    ax.set_ylabel(index)
    ax.set_title(f'{values} step:{str(n)}')

    contour = ax.contourf(
        *np.meshgrid(*[x.values.astype(np.float32) for x in [table.columns, table.index]]), table.values,
        cmap='jet', levels=200, alpha=0.9
    )
    fig.colorbar(contour)
    contour.set_clim(vmin=data[values].min(), vmax=63)
    contour.set_clim(data[values].min(), 63)

    #plt.scatter(plots[columns], plots[index], s=5, c='none', edgecolors='purple', alpha=0.7, linewidth=0.4)
    plt.scatter(next_step[columns], next_step[index], s=5, c='pink', edgecolors='pink')

    plt.xlim(0.0, 3.0)
    plt.ylim(0.0, 3.0)
    
    fig.savefig(f'{res}/{values}_{str(n)}.png')
    plt.close()


def draw_graph(res, keys):
    """ draw graphs """
    steps_df = pd.read_csv(res, header=0)
    for n, step in enumerate(steps_df.T):
        if n > 35 and 5*n < len(steps_df):
            data = steps_df[0:5*n]
            next_step = steps_df[5*n-1:5*n]
            _ = [add_graph('grid/', n, data, next_step, **kargs) for kargs in [
                {'values':'eff', 'index':keys[0], 'columns':keys[1]}, # 探索点の計算値 score
            ]]

    make_gif(res, keys)

def make_gif(res, keys):
    """ make gif """
    image_dir = 'grid'
    os.makedirs('grid/gifs', exist_ok=True)
    image_path = ["g"]

    for path in image_path:
        X = []
        img_paths = os.listdir(image_dir + '/' + path)
        img_paths = sorted(img_paths, key=lambda s: int(re.search(r'\d+', s).group()))

        for i, img_name in enumerate(img_paths):
            if i < 4001:
                X.append(Image.open(image_dir + '/' + path + '/' + img_name))

        X[0].save(image_dir + '/gifs/' + path + '.gif', save_all=True, append_images=X[1:], optimize=False, duration=200, loop=1)

if __name__ == "__main__":
    steps_path = 'DB_map_max.csv'
    keys_list = ['Eic', 'Eg']
    #draw_graph(steps_path, keys_list)
    make_gif(steps_path, keys_list)