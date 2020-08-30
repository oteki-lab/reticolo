""" reticolo wrapper """
#!/usr/bin/env python3

import copy
import itertools
import csv
from importlib import import_module
from importlib.abc import MetaPathFinder
from importlib.util import spec_from_file_location
import sys
import argparse
import pickle
import numpy, scipy.io

def write_csv(csv_file, save_dict):
    """Save parameter list in mat"""
    scipy.io.savemat('input.mat', mdict={'input_data': save_dict})
    print("Output combination list file")

class Listing:
    """Parse command line arguments and Take an input file and outputs a matrix with a column"""
    def __init__(self, input_file):
        # Set up
        for k in ["in_param", "ite_list"]:
            setattr(self, k, None)

        def import_file(name, location):
            """Load input file as module"""
            class Finder(MetaPathFinder):
                """Load module class"""
                @staticmethod
                def find_spec(fullname, *_):
                    """Load module method"""
                    if fullname == name:
                        return spec_from_file_location(name, location)
            finder = Finder()
            sys.meta_path.insert(0, finder)
            try:
                return import_module(name)
            finally:
                sys.meta_path.remove(finder)

        ### Tacking all parameters by parsing command line arguments
        # load input file as module
        mod = import_file('input_file', input_file)
        # All parameters in module
        all_para = dict((k, v) for k, v in mod.__dict__.items() if isinstance(v, (int, float, list, tuple, str, bool, dict)) and not k.startswith("__"))
        # All static parameters in parameters
        in_param = dict((k, v) for k, v in all_para.items() if not hasattr(v, '__iter__'))

        ### Listing iterable parameters
        # All iterable parameters that vary for each experiment
        ite_vals = dict((k, sorted(v)) for k, v in all_para.items() if k not in in_param)
        # All keys of iterable parameters
        ite_keys = ite_vals.keys()
        # Generate all combinations of iterable parameters
        ite_comb = list(itertools.product(*[ite_vals[k] for k in ite_keys]))
        # list combinations of iterable parameters for each experiment
        ite_list = [dict(zip(ite_keys, values)) for values in ite_comb]

        ### Listing combinations of parameters for each experiment
        params_list = []
        for _, ite_para in enumerate(copy.deepcopy(ite_list)):
            ite_para.update(in_param)
            params_list.append(ite_para)
        write_csv("file.pickle", params_list)

def main(arguments=None):
    """Main function"""
    def parsing(arguments):
        """Parsing arguments"""
        parser = argparse.ArgumentParser(description="Calculates IBSC in DB model.")
        parser.add_argument('input_file', nargs="?", default="input.py")
        return parser.parse_args(arguments.split() if arguments else None)

    args = parsing(arguments)
    try:
        Listing(args.input_file)
    except NameError:
        print("Input Files Error")

if __name__ == "__main__":
    main() #py main.py input.py
