""" combination """
import itertools
import numpy as np

def combine(_params):
    """get params"""
    x_dict = {}
    for item in _params:
        x_dict[item['name']] = np.linspace(item['min'], item['max'], item['div'])

    ### Tacking all parameters by parsing command line arguments
    # All static parameters in parameters
    in_param = dict((k, v) for k, v in x_dict.items() if not hasattr(v, '__iter__'))

    ### Listing iterable parameters
    # All iterable parameters that vary for each experiment
    ite_vals = dict((k, sorted(v)) for k, v in x_dict.items() if k not in in_param)
    # All keys of iterable parameters
    ite_keys = ite_vals.keys()
    # Generate all combinations of iterable parameters
    ite_comb = list(itertools.product(*[ite_vals[k] for k in ite_keys]))
    # list combinations of iterable parameters for each experiment
    ite_list = [dict() for values in ite_comb]

    ite_list = [dict(zip(ite_keys, tuple(map(lambda x: isinstance(x, float) and round(x, 5) or x, values)))) for values in ite_comb]
    return ite_list
    