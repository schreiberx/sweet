#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import inspect

from mule.sdc.qmatrix import getSetup

# Extract parameter documentation from function docstring
doc = getSetup.__doc__.splitlines()
paramsDoc = {}
for i in range(len(doc)):
    if doc[i].strip().startswith('Parameters'):
        break
i += 2
while i < len(doc):
    if doc[i].strip().startswith('Example'):
        break
    if doc[i].strip().endswith('optional'):
        name = doc[i].strip().split(":")[0].strip()
        paramsDoc[name] = []
        i += 1
        while doc[i].startswith("        ") or doc[i].strip() == "":
            if doc[i].strip() != "":
                paramsDoc[name].append(doc[i])
            i += 1
        paramsDoc[name] = '\n'.join(paramsDoc[name])

# Specific formater with a bit more width ...
class HelpFormatter(argparse.RawTextHelpFormatter):
    def __init__(
        self, prog, indent_increment=2, max_help_position=36, width=None):
        super().__init__(
            prog, indent_increment, max_help_position, width)


# Program parser
parser = argparse.ArgumentParser(
    prog='generateSWEETFileDict_SDC',
    description="Generate SWEETFileDict parameter file for SDC",
    formatter_class=HelpFormatter)

setup = inspect.signature(getSetup)

parser.add_argument(
    f'--showDefault', action='store_true', help="show default configuration")

parser.add_argument(
    f'--fileName', default="params_SDC.sweet", 
    help="name of the SWEETFileDict file (default: params_SDC.sweet)")

parser.add_argument(
    f'--preset', default=None, 
    help="unique ID for preset parameters")

parser.add_argument(
    f'--showPreset', action='store_true', help="show preset configurations")

PRESET_LIST = {
    'P1' : {
        'nNodes': 4,
        'nodeType': 'RADAU-RIGHT',
        'nIter': 3,
        'qDeltaImplicit': 'BEPAR',
        'qDeltaExplicit': 'PIC',
        'initialSweepType': 'QDELTA',
        'useEndUpdate': False,
        'qDeltaInitial': 'BEPAR',
        'diagonal': True,
    },
    'P2' : {
        'nNodes': 3,
        'nodeType': 'RADAU-RIGHT',
        'nIter': 3,
        'qDeltaImplicit': 'OPT-QMQD-0',
        'qDeltaExplicit': 'PIC',
        'initialSweepType': 'QDELTA',
        'useEndUpdate': False,
        'qDeltaInitial': 'BEPAR',
        'diagonal': True,
    },
    'P3' : {
        'nNodes': 3,
        'nodeType': 'RADAU-RIGHT',
        'nIter': 3,
        'qDeltaImplicit': 'BEPAR',
        'qDeltaExplicit': 'PIC',
        'initialSweepType': 'QDELTA',
        'useEndUpdate': True,
        'qDeltaInitial': 'BEPAR',
        'diagonal': True,
    },
}

# Add parser argument for each function parameter
defaults = {}
for name, val in setup.parameters.items(): 
    parser.add_argument(
        f'--{name}', default=None, type=val.annotation,
        help=paramsDoc[name])
    defaults[name] = val.default

args = parser.parse_args()

if args.showDefault:
    print('-'*80)
    print('Preset configurations')
    print('-'*80)
    for name, val in defaults.items(): 
        print(f'--{name} {val}')
    print('-'*80)
    

if args.showPreset:
    print('-'*80)
    print('Preset configurations')
    print('-'*80)
    for key, dico in PRESET_LIST.items():
        print(f'{key} :')
        for name, value in dico.items():
            print(f'--{name} {value}')
        print('-'*80)

if args.showPreset or args.showDefault:
    exit()

params = {name: val for name, val in args._get_kwargs()}
params.pop('fileName')
params.pop('preset')
params.pop('showPreset')
params.pop('showDefault')

if args.preset is not None:
    preset = args.preset
    if preset in PRESET_LIST.keys():
        for name, val in params.items():
            if val is not None:
                PRESET_LIST[preset][name] = val
        params = PRESET_LIST[preset]
    else:
        raise ValueError(f'no {name} preset config found, choose between {list(PRESET_LIST.keys())}')
for name, val in params.items():
    if val is None:
        params[name] = defaults[name]
params = getSetup(**params)
params.writeToFile(args.fileName)
print(f'Written SDC parameter file : {args.fileName}')
