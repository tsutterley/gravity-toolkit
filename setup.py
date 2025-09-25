import os
from setuptools import setup, find_packages

# list of all scripts to be included with package
scripts = []
for dir in ['access','dealiasing','geocenter','mapping','scripts','utilities']:
    scripts.extend([os.path.join(dir,f) for f in os.listdir(dir) if f.endswith('.py')])
scripts.append(os.path.join('gravity_toolkit','grace_date.py'))
scripts.append(os.path.join('gravity_toolkit','grace_months_index.py'))

setup(
    name='gravity-toolkit',
    scripts=scripts,
)
