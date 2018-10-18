#! /usr/bin/env python3

import glob
import re

travis_file=".travis.yml"

print("Setting up tests in travis file '"+travis_file+"'")


files = glob.glob('./tests/??_*/test.sh')
files += glob.glob('./tests/??_*/test.py')

tests = ""

for f in files:
	print(" + Found test script "+f)
	tests += '  - '+f+" || exit 1\n"


with open(travis_file, 'r') as f:
	travis_content = f.read()

travis_content = re.sub(r'(###SCRIPT_START###).*(###SCRIPT_END###)', "\\1\nscript:\n"+tests+"\\2", travis_content, flags=re.MULTILINE | re.DOTALL)

print("Writing content to file '"+travis_file+"'")
with open(travis_file, 'w') as f:
	f.write(travis_content)
