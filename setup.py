import os
from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('requirements.txt') as fh:
    install_requires = fh.read().splitlines()

setup(
    name='read-GRACE-harmonics',
    version='1.0.1.14',
    description='Reads Level-2 spherical harmonic coefficients from the NASA/DLR GRACE and NASA/GFZ GRACE Follow-on missions',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/tsutterley/read-GRACE-harmonics',
    author='Tyler Sutterley',
    author_email='tsutterl@uw.edu',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
    ],
    keywords='GRACE, GRACE-FO, Gravity, satellite geodesy, spherical harmonics',
    packages=find_packages(),
    install_requires=install_requires,
    dependency_links=['https://github.com/tsutterley/read-GRACE-geocenter/tarball/master'],
    scripts=[os.path.join('scripts',f) for f in os.listdir('scripts')],
)
