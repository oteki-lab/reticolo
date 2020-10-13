""" bayesian optimization main """
from high_params import params
from step_list import steps
from combination import combine
from bayesian_optimization import search

if __name__ == "__main__":

    ### 探索項目の整形 ###
    for item in params:
        item['div'] = int((item['max']-item['min'])/item['step'])+1     # 探索分割数
        item['space'] = 1/(item['div']-1)                               # 探索間隔

    combi = combine(params)

    res = search(params, combi, steps)
    print(res)
