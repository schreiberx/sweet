#! /usr/bin/env python3

from mule.JobGeneration import *

jg = JobGeneration()

#
# Compile options
#
jg.compile.program = 'tutorial_classDict_ProgramArguments'

# Set some user parameter - just because we can
jg.runtime.user_defined_parameters['grav_param'] = {'id': '', "option": "--gravitation=", "value": 666}

# Generate jobs
jg.gen_jobscript_directory()
