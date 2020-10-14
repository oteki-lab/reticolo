""" bayesian optimization main """
from bayesian import search

if __name__ == "__main__":
    params = 'Results\\20202114191036805\\hp_table.csv'
    steps = 'Results\\20202114191036805\\hp_steps.csv'
    res = search(params, steps)
