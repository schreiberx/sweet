#! /usr/bin/env python3

import sys
import os

if not 'MULE_SOFTWARE_ROOT' in os.environ:
    print("No SWEET environment variables detected, skipping Travis updates")
    sys.exit(0)

os.chdir(os.environ['MULE_SOFTWARE_ROOT'])


import glob
import re

from itertools import product

travis_file=".travis.yml"

verbosity = 10
if len(sys.argv) > 1:
    verbosity = int(sys.argv[1])

if verbosity >= 10:
    print("Working directory: "+os.path.abspath(os.curdir))
    print("Setting up tests in travis file '"+travis_file+"'")

tests = glob.glob('./tests/??_*/test.sh')
tests += glob.glob('./tests/??_*/test.py')


if verbosity >= 10:
    for test in tests:
        print(" + Found test script '"+test+"'")



if verbosity >= 10:
    print("Writing content to file '"+travis_file+"'")

with open(travis_file, 'w') as f:
    f.write("""#
# Script for Travis CI
#
# See doc/travis_ci.txt for more information
#

language: cpp

dist: trusty

#
# We want to setup the 3rd party libraries and test SWEET
# with different software consellations
#
# This is called a build matrix and generates different build
# environments
#
jobs:
  include:

""")

    jobs_list = []


    if False:
    #if True:
        jobs_list += [
"""
    # Test with G++-9
    - os: linux
      addons:
        apt:
          sources:
            - ubuntu-toolchain-r-test
          packages:
            - texinfo
            - g++-9
            - gfortran-9
      env:
        - MATRIX_EVAL="export CC=gcc-9 && export CXX=g++-9 && export FC=gfortran-9 && export F90=gfortran-9"
"""]
    #if False:
    if True:
        # texinfo is required for compiling 'make'
        jobs_list += [
"""
    # Test with G++-8
    - os: linux
      addons:
        apt:
          sources:
            - ubuntu-toolchain-r-test
          packages:
            - texinfo
            - g++-8
            - gfortran-8
      env:
        - MATRIX_EVAL="export CC=gcc-8 && export CXX=g++-8 && export FC=gfortran-8 && export F90=gfortran-8"
"""]

    if False:
    #if True:
        jobs_list += [
"""
    # Test with G++-7
    - os: linux
      addons:
        apt:
          sources:
            - ubuntu-toolchain-r-test
          packages:
            - texinfo
            - g++-7
            - gfortran-7
      env:
        - MATRIX_EVAL="export CC=gcc-7 && export CXX=g++-7 && export FC=gfortran-7 && export F90=gfortran-7"
"""]


    c = 0
    for (j, test) in product(jobs_list, tests):
        if True:
            # This version allows reutilizing the cache
            f.write(j)
            if True:
                f.write("      script: "+test)
            else:
                f.write("      script:\n")
                f.write("        - cd \""+os.path.dirname(test)+"\"\n")
                f.write("        - ./"+os.path.basename(test)+"\n")

        else:
            j = j.replace('MATRIX_EVAL="', 'MATRIX_EVAL="TESTSCRIPT='+test+' && ')
            f.write("      script: $TESTSCRIPT")
        f.write("\n")
        f.write("\n")
        c += 1


    f.write("""

#
# Install dependencies
#
# See https://docs.travis-ci.com/user/installing-dependencies/
#
before_install:
  # Load matrix environment
  - echo "${MATRIX_EVAL}"
  - eval "${MATRIX_EVAL}"

  # Debug output
  - hostname

  # Load SEET environment variables
  - cd local_software || exit 1
  - source env_vars.sh || exit 1

  # Setup additional SWEET software packages
  - ./setup_local_software.sh || exit 1

  # Go back to SWEET's root directory
  - cd $MULE_SOFTWARE_ROOT



#
# SWEET requires binaries compiled individually for each test
# Skip installation phase (install: true)
#
install: true



#
# Cache installed software
#
# After restoring the cache, the install scripts check for found
# software and avoid recompiling and installing it.
#
# See https://docs.travis-ci.com/user/caching/
#
# The cache is setup amongst others based on the environment variables
cache:
  directories:
    # Cache the install directory for SWEET's 3rd party software
    local_software/local

""")

if verbosity >= 10:
    print("Job combinations: "+str(c))
