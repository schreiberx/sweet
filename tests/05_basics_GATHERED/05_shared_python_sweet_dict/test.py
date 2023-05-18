#! /usr/bin/env python3

import sys
import os
import numpy as np

from mule.utils import exec_program
exec_program('mule.benchmark.cleanup_all', catch_output=False)

#Dummy email for tests on the linux cluster platform
os.environ['MULE_USER_EMAIL']='dummy_email@linuxcluster.com'

from mule.SWEETFileDict import *

shared_dict = SWEETFileDict()

shared_dict.set("str", "hello world")

shared_dict.set("int64", 123)
shared_dict.set("float64", 123.456)
shared_dict.set("complex128", 123.456+(123.456+100)*1j)

shared_dict.set("array_float64_1d", np.array([1., 2., 3.]))
shared_dict.set("array_float64_2d", np.array([[1., 2., 3.], [4., 5., 6.]]))
shared_dict.set("array_float64_3d", np.array([[[1., 2., 3.], [4., 5., 6.]], [[7., 8., 9.], [10., 11., 12.]]]))

shared_dict.set("array_complex128_1d", np.array([1.+102j, 2.+103j, 3.+104j]))
shared_dict.set("array_complex128_2d", np.array([[1.+102j, 2.+103j, 3.+104j], [4.+105j, 5.+106j, 6.+107j]]))
shared_dict.set("array_complex128_3d", np.array([[[1.+102j, 2.+103j, 3.+104j], [4.+105j, 5.+106j, 6.+107j]], [[7.+108j, 8.+109j, 9.+110j], [10.+111j, 11.+112j, 12.+113j]]]))

print(shared_dict)

shared_dict_filename = "output_shared_dict.sweet"

print("*"*80)
print(f"Writing to file '{shared_dict_filename}'")
shared_dict.writeToFile(shared_dict_filename)
print("*"*80)

print("*"*80)
print(f"Loading from file '{shared_dict_filename}'")
new_dict = SWEETFileDict(shared_dict_filename)
print("*"*80)

print(new_dict)


if shared_dict == new_dict:
    print("Both dictionaries are the same! SUCCESS!")
else:
    raise Exception("Dictionary problem!")




"""
Next, we run the unit test from the SWEET C++ side
"""

from mule.JobMule import *
from mule.InfoError import *


jg = JobGeneration()
jg.compile.mode = "debug"
jg.compile.program = "tests/core_fileDict"
jg.runtime.user_defined_parameters['sweet-file-dict'] = {'id': '', 'value': "../"+shared_dict_filename, 'option': '--sweet-file-dict='}
jg.runtime.user_defined_parameters['unit-test-sweet-file-dict'] = {'id': '', 'value': "1", 'option': '--run-unit-test-sweet-file-dict='}

jg.gen_jobscript_directory()

exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
    sys.exit(exitcode)

print("Benchmarks successfully finished")

exec_program('mule.benchmark.cleanup_all', catch_output=False)

