from setuptools import setup, find_packages
setup(
    name='read-GRACE-harmonics',
    version='1.0.1.11',
    description='Reads Level-2 spherical harmonic coefficients from the NASA/DLR GRACE and NASA/GFZ GRACE Follow-on missions',
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
    install_requires=['numpy','pyyaml','lxml','future','matplotlib','cartopy','netCDF4','h5py'],
    dependency_links=['https://github.com/tsutterley/read-GRACE-geocenter/tarball/master'],
)
