#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import inspect

from mule.sdc.qmatrix import getSetup

descr = getSetup.__doc__

parser = argparse.ArgumentParser(
    prog='generateSWEETFileDict_SDC',
    description=descr,
    formatter_class=argparse.RawDescriptionHelpFormatter)

setup = inspect.signature(getSetup)

parser.add_argument(
    f'--fileName', default="paramsSDC.sweet", 
    help="Name of the SWEETFileDict file (default: paramsSDC.sweet)")

for name, val in setup.parameters.items(): 
    parser.add_argument(
        f'--{name}', default=val.default, type=val.annotation,
        help='see documentation above')

args = parser.parse_args()

params = {name: val for name, val in args._get_kwargs()}
params.pop('fileName')

params = getSetup(**params)
params.writeToFile('paramsSDC.sweet')
