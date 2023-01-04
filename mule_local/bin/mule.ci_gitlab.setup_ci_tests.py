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

gitlab_ci_file=".gitlab-ci.yml"

verbosity = 10
if len(sys.argv) > 1:
    verbosity = int(sys.argv[1])


if verbosity >= 10:
    print("Working directory: "+os.path.abspath(os.curdir))
    print("Setting up tests in Gitlab CI file '"+gitlab_ci_file+"'")

tests = glob.glob('./tests/??_*/test.sh')
tests += glob.glob('./tests/??_*/test.py')


if verbosity >= 10:
    for test in tests:
        print(" + Found test script '"+test+"'")



if verbosity >= 10:
    print("Writing content to file '"+gitlab_ci_file+"'")

with open(gitlab_ci_file, 'w') as f:
    f.write("""
#
# Gitlab CI file for SWEET
#

# bionic: 18.04
# focal: 20.04
# jammy: 22.04
image: ubuntu:focal

before_script:
    # Install required packages
    # Avoid interaction by setting 
    - export DEBIAN_FRONTEND=noninteractive
    # Install packages
    - apt-get update -qq && apt-get install -y -qq git make automake g++ cmake libfftw3-dev libopenmpi-dev libeigen3-dev liblapack-dev numactl python3 curl

stages:          # List of stages for jobs, and their order of execution
  - setup-software
  - tests


setup-job:       # This job runs in the build stage, which runs first.
  stage: setup-software
  script:
    # Setup SWEET

    # Get only head by specifying -- depth 1
    - git clone --depth 1 https://github.com/schreiberx/sweet.git
    - cd sweet
    - source ./activate.sh
    - echo "SWEET environment activated"

    # Compile other required software
    - echo "Fetching and compiling software"
    - cd local_software
    - ./setup_local_software.sh
    - cd ../

    # Cleaning up things a little bit
    - echo "Cleaning up..."
    - rm -rf .git
    - rm -rf local_software/local_src
    - echo "Finished"

    # We're ready to use SWEET
  
  # Now we just want to store all the software as an artifact to use it for the next steps
  artifacts:
    # Just get everything for the next stage
    untracked: true
    expire_in: 1 day
  
  #artifacts:
  #  paths:
  #    - sweet/local_software/local/


""")

    c = 0
    for test in tests:

        f.write(f"""

job-test-{c}:   # This job runs in the test stage.
  stage: tests  # It only starts when the job in the build stage completes successfully.
  dependencies:
    - "setup-job"
  script:
    - cd sweet
    - source ./activate.sh
    - {test}
    - echo "TEST {c} SUCCESSFULLY FINISHED"

""")
        c += 1


if verbosity >= 10:
    print("Job combinations: "+str(c))
