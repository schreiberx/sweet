#! /usr/bin/env python3

import sys
import os

if not 'MULE_SOFTWARE_ROOT' in os.environ:
    print("No SWEET environment variables detected, skipping CI updates")
    sys.exit(0)



# Which gitlab CI to use
#gitlab_ci_kind = 'gricad'
gitlab_ci_kind = 'inria'



"""
List of dictionaries describing the test environments
"""
test_environment_ = []


ubuntu_version_ = [18, 20, 22]

# Setup configurations for all different ubuntu versions and compilers
for ubuntu_version in ubuntu_version_:

    # Available Ubuntu compiler versions
    if ubuntu_version == 18:
        gcc_version_ = [6,7,8]
    elif ubuntu_version == 20:
        gcc_version_ = [8,9,10]
    elif ubuntu_version == 22:
        gcc_version_ = [9,10,11,12]

    # Particular image version - named differently on different CI systems
    image_version = None
    if gitlab_ci_kind == "gricad":
        # bionic: 18.04
        # focal: 20.04
        # jammy: 22.04
        if ubuntu_version == 18:
            image_version = "ubuntu:bionic"
        elif ubuntu_version == 20:
            image_version = "ubuntu:focal"
        elif ubuntu_version == 22:
            image_version = "ubuntu:jammy"
        else:
            raise Exception("Unsupported version")

    elif gitlab_ci_kind == "inria":
        if ubuntu_version == 18:
            image_version = "library/ubuntu:18.04"
        elif ubuntu_version == 20:
            image_version = "library/ubuntu:20.04"
        elif ubuntu_version == 22:
            image_version = "library/ubuntu:22.04"
        else:
            raise Exception("Unsupported version")

    # Graphics packages vary in Ubuntu versions
    apt_get_graphics_pkgs = None
    if ubuntu_version == 18:
        apt_get_graphics_pkgs = "pkg-config libfreetype6-dev libglu1-mesa-dev"
    elif ubuntu_version == 20:
        apt_get_graphics_pkgs = "pkg-config libfreetype-dev libopengl-dev libglu1-mesa-dev libxext-dev"
    elif ubuntu_version == 22:
        apt_get_graphics_pkgs = "pkg-config libfreetype-dev libopengl-dev libglu1-mesa-dev libxext-dev"
    else:
        raise Exception("Unsupported version")

    # Special tags for INRIA CI
    tags = ""
    if gitlab_ci_kind == "inria":
        tags = """\
  tags:
    - ci.inria.fr
    - large
"""


    # Standard apt-get packages
    apt_get_packages = "git make automake cmake python3 wget"

    for gcc_version in gcc_version_:

        #
        # Generate individual CI tests
        # only for these particular environments
        #
        gen_ci_tests = False
        if ubuntu_version == 22 and gcc_version == 12:
            gen_ci_tests = True
        #elif ubuntu_version == 20 and gcc_version == 10:
        #    gen_ci_tests = True
        #elif ubuntu_version == 18 and gcc_version == 8:
        #    gen_ci_tests = True
        elif ubuntu_version == 18 and gcc_version == 6:
            gen_ci_tests = True

        test_environment_ += [{
                    'id': f"ubuntu{ubuntu_version}-gcc{gcc_version}",
                    'ubuntu_version': ubuntu_version,
                    'gcc_version': gcc_version,
                    'apt_get': f"apt-get install -y -qq {apt_get_packages} g++-{gcc_version} gcc-{gcc_version} gfortran-{gcc_version} {apt_get_graphics_pkgs}",
                    'image_version': image_version,
                    'tags': tags,
                    'gen_ci_tests': gen_ci_tests,
                }]




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

# Sort tests
tests.sort()



job_counter = 0

if verbosity >= 10:
    for test in tests:
        print(" + Found test script '"+test+"'")

if verbosity >= 10:
    print("Writing content to file '"+gitlab_ci_file+"'")


content = ""

content += """\
#
# Gitlab CI file for SWEET
#

stages:          # List of stages for jobs, and their order of execution
  - setup_local_software
"""

for te in test_environment_:
    content += f"""\
  - stage_tests_{te['id']}
"""

content += f"""\

"""

#
# Setup local software for each environment
#
for te in test_environment_:

    script_header = f"""\
    # Avoid interaction by setting 
    - export DEBIAN_FRONTEND=noninteractive
    # Install packages
    - apt-get update -qq
    - {te['apt_get']}
    - export CC=gcc-{te['gcc_version']}
    - export CXX=g++-{te['gcc_version']}
    - export F90=gfortran-{te['gcc_version']}
    - export FC=gfortran-{te['gcc_version']}
    - export LD=ld
"""

    content += f"""\

setup-local-software-{te['id']}:       # This job runs in the build stage, which runs first.
  stage: setup_local_software
  image: {te['image_version']}
{te['tags']}\
  script:
{script_header}\
    # Setup SWEET

    # Get only head by specifying -- depth 1
    - git clone --depth 1 https://github.com/schreiberx/sweet.git
    - cd sweet
    - source ./activate.sh
    - echo "SWEET environment activated"

    # Compile other required software
    - echo "Fetching and compiling software"
    - cd local_software
    - time ./setup_local_software.sh
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
    name: $CI_JOB_NAME
  
"""

    job_counter += 1


#
# Generate individual test cases
#

#for te in test_environment_:
    if te['gen_ci_tests']:
        for test in tests:

            r = re.match(r".*/([0-9_a-z]*)/test.*", test)
            job_id = r.group(1)

            content += f"""\

job-test-{te['id']}-{job_id}:
  stage: stage_tests_{te['id']}
  image: {te['image_version']}
{te['tags']}\
  # Job dependency
  needs: ["setup-local-software-{te['id']}"]
  # Download files from artifact of this build
  dependencies:
    - "setup-local-software-{te['id']}"
  script:
{script_header}\
    - ls
    - cd sweet
    - ls
    - source ./activate.sh
    - ls
    - {test}
    - echo "TEST {job_id} SUCCESSFULLY FINISHED"

"""
            job_counter += 1


print(content)


with open(gitlab_ci_file, 'w') as f:
    f.write(content)


if verbosity >= 10:
    print(f"Job counter: {job_counter}")

