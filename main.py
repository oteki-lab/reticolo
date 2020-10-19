""" bayesian optimization main """
import os
import re
from PIL import Image
from plot_module import search

if __name__ == "__main__":
    res_dir = 'Results\\20202619111027347'

    params = res_dir+'\\hp_table.csv'
    steps = res_dir+'\\hp_steps.csv'

    image_dir = res_dir+"/graphs/"
    image_path = ["score", "y_mean", "acq"]

    res = search(res_dir, params, steps)


    os.makedirs(image_dir + '/gifs', exist_ok=True)
    for path in image_path:
        X = []
        img_paths = os.listdir(image_dir + '/' + path)
        img_paths = sorted(img_paths, key=lambda s: int(re.search(r'\d+', s).group()))

        for img_name in img_paths:
            X.append(Image.open(image_dir + '/' + path + '/' + img_name))

        X[0].save(image_dir + '/gifs/' + path + '.gif', save_all=True, append_images=X[1:], optimize=False, duration=1000, loop=0)
