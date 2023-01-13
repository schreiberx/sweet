#! /bin/bash


# Removing all python cache files (__pycache__, *.pyc, *.pyo)
find . | grep -E "(/__pycache__$|\.pyc$|\.pyo$)" | xargs rm -rf

# Removing all empty directories
find . -empty -type d -delete

