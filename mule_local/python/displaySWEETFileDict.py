#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse

from mule.SWEETFileDict import SWEETFileDict

parser = argparse.ArgumentParser(
    prog='displaySWEETFileDict',
    description="Display the content of a SWEETFileDict binary file",
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument(
    'file', help="The SWEETFileDict to display"
)

args = parser.parse_args()

params = SWEETFileDict(args.file)

print(f'Content of SWEETFileDict : {args.file}')
for key, val in params.dict.items():
    print('-'*80)
    val = str(val)
    if '\n' in val:
        print(f' -- {key}:\n{val}')
    else:
        print(f' -- {key}: {val}')
print('-'*80)