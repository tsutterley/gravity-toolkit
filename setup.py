from setuptools import setup, find_packages
setup(
	name='read-GRACE-harmonics',
	version='1.0.0.3',
	description='Reads Level-2 spherical harmonic coefficients from the NASA/DLR Gravity Recovery and Climate Experiment (GRACE) mission',
	url='https://github.com/tsutterley/read-GRACE-harmonics',
	author='Tyler Sutterley',
	author_email='tyler.c.sutterley@nasa.gov',
	license='MIT',
	classifiers=[
		'Development Status :: 3 - Alpha',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Physics',
		'License :: OSI Approved :: MIT License',
		'Programming Language :: Python :: 2',
		'Programming Language :: Python :: 2.7',
	],
	keywords='GRACE Gravity Spherical Harmonics',
	packages=find_packages(),
	install_requires=['numpy','pyyaml'],
)
