#!/usr/bin/env python
u"""
read_GRACE_harmonics.py
Written by Tyler Sutterley (05/2018)

Reads GRACE datafile and extracts spherical harmonic data and drift rates (RL04)
Adds drift rates to clm and slm for release 4 harmonics
Correct GSM data for drift in pole tide following Wahr et al. (2015)
Extracts date of GRACE datafile from file name and calculates mean of range

INPUTS:
	input_file: GRACE Level-2 spherical harmonic datafile
	LMAX: Maximum degree of spherical harmonics (degree of truncation)

OPTIONS:
	MMAX: Maximum order of spherical harmonics (order of truncation)
		default is the maximum spherical harmonic degree
	POLE_TIDE: correct GSM data for pole tide drift following Wahr et al. (2015)

OUTPUTS:
	time: mid-month date in year-decimal
	start: start date of range as Julian day
	end: end date of range as Julian day
	clm: cosine spherical harmonics of input data (LMAX,MMAX)
	slm: sine spherical harmonics of input data (LMAX,MMAX)
	eclm: cosine spherical harmonic uncalibrated standard deviations (LMAX,MMAX)
	eslm: sine spherical harmonic uncalibrated standard deviations (LMAX,MMAX)

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python (http://www.numpy.org)

UPDATE HISTORY:
	Updated 05/2018: updates to file name structure with release 6 and GRACE-FO
	Written 10/2017 for public release
"""
import os
import re
import gzip
import numpy as np

#-- PURPOSE: read Level-2 GRACE spherical harmonic files
def read_GRACE_harmonics(input_file, LMAX, MMAX=None, POLE_TIDE=False):
	#-- compile numerical expression operator for parameters from files
	regex_pattern = ('(.*?)-2_(\d+)-(\d+)_(.*?)_(UTCSR|EIGEN|JPLEM|JPLMSC)_'
		'(.*?)_(\d+)(.*?)(\.gz|\.gfc)?$')
	rx = re.compile(regex_pattern, re.VERBOSE)
	#-- extract parameters from input filename
	PFX,SD,ED,N,PRC,F1,DRL,F2,SFX=rx.findall(os.path.basename(input_file)).pop()

	#-- JPL Mascon solutions
	if PRC in ('JPLMSC'):
		DSET = 'GSM'
		DREL = np.int(DRL)
		FLAG = 'GRCOF2'
	#-- Kusche et al. (2009) DDK filtered solutions 10.1007/s00190-009-0308-3
	elif PFX.startswith('kfilter_DDK'):
		DSET = 'GSM'
		DREL = np.int(DRL)
		FLAG = 'gfc'
	#-- Standard GRACE solutions
	else:
		DSET = PFX
		DREL = np.int(DRL)
		FLAG = 'GRCOF2'

	#-- output python dictionary with GRACE data and date information
	grace_L2_input = {}
	#-- extract GRACE date information from input file name
	start_yr = np.float(SD[:4])
	end_yr = np.float(ED[:4])
	start_day = np.float(SD[4:])
	end_day = np.float(ED[4:])
	#-- calculate mid-month date taking into account if measurements are
	#-- on different years
	if ((start_yr % 4) == 0):#-- Leap Year (% = modulus)
		dpy = 366.0
	else:#-- Standard Year
		dpy = 365.0
	#-- For data that crosses years (end_yr - start_yr should be at most 1)
	end_day = ((end_yr-start_yr)*dpy+end_day) if (start_yr!=end_yr) else end_day
	#-- Calculation of Mid-month value
	mid_day = np.mean([start_day, end_day])
	#-- Calculating the mid-month date in decimal form
	grace_L2_input['time'] = start_yr + mid_day/dpy
	#-- Calculating the Julian dates of the start and end date
	grace_L2_input['start'] = np.float(367.0*start_yr - np.floor(7.0*start_yr/4.0) -
		np.floor(3.0*(np.floor((start_yr - 8.0/7.0)/100.0) + 1.0)/4.0) +\
		np.floor(275.0/9.0) + start_day + 1721028.5)
	grace_L2_input['end'] = np.float(367.0*end_yr - np.floor(7.0*end_yr/4.0) -
		np.floor(3.0*(np.floor((end_yr - 8.0/7.0)/100.0) + 1.0)/4.0) +
		np.floor(275.0/9.0) + end_day + 1721028.5)

	#-- set maximum spherical harmonic order
	MMAX = np.copy(LMAX) if (MMAX is None) else MMAX
	#-- Spherical harmonic coefficient matrices to be filled from data file
	grace_L2_input['clm'] = np.zeros((LMAX+1,MMAX+1))
	grace_L2_input['slm'] = np.zeros((LMAX+1,MMAX+1))
	#-- spherical harmonic uncalibrated standard deviations
	grace_L2_input['eclm'] = np.zeros((LMAX+1,MMAX+1))
	grace_L2_input['eslm'] = np.zeros((LMAX+1,MMAX+1))
	if ((DREL == 4) and (DSET == 'GSM')):
		#-- clm and slm drift rates for RL04
		drift_c = np.zeros((LMAX+1,MMAX+1))
		drift_s = np.zeros((LMAX+1,MMAX+1))

	#-- Opening data file to extract spherical harmonic coefficients
	#-- check if file is compressed (read with gzip if gz)
	if (SFX == '.gz'):
		#-- GRACE file is compressed (gz) file
		with gzip.open(os.path.expanduser(input_file),'r') as f:
			file_contents = f.read().splitlines()
	else:
		#-- GRACE file is standard ascii file
		with open(os.path.expanduser(input_file),'r') as f:
			file_contents = f.read().splitlines()

	#-- for each line in the GRACE file
	for line in file_contents:
		#-- find if line starts with data marker flag (e.g. GRCOF2)
		if bool(re.match(FLAG,line)):
			#-- split the line into individual components
			line_contents = line.split()
			#-- degree and order for the line
			l1 = np.int(line_contents[1])
			m1 = np.int(line_contents[2])
			#-- if degree and order are below the truncation limits
			if ((l1 <= LMAX) and (m1 <= MMAX)):
				grace_L2_input['clm'][l1,m1] = np.float(line_contents[3])
				grace_L2_input['slm'][l1,m1] = np.float(line_contents[4])
				grace_L2_input['eclm'][l1,m1] = np.float(line_contents[5])
				grace_L2_input['eslm'][l1,m1] = np.float(line_contents[6])
		#-- find if line starts with drift rate flag
		elif bool(re.match('GRDOTA',line)):
			#-- split the line into individual components
			line_contents = line.split()
			l1 = np.int(line_contents[1])
			m1 = np.int(line_contents[2])
			#-- Reading Drift rates for low degree harmonics
			drift_c[l1,m1] = np.float(line_contents[3])
			drift_s[l1,m1] = np.float(line_contents[4])

	#-- Adding drift rates to clm and slm for RL04
	#-- if drift rates exist at any time, will add to harmonics
	#-- Will convert the secular rates into a stokes contribution
	#-- Currently removes 2003.3 to get the temporal average close to 0.
	#-- note: += means grace_xlm = grace_xlm + drift_x
	if ((DREL == 4) and (DSET == 'GSM')):
		#-- time since 2003.3
		dt = (grace_L2_input['time']-2003.3)
		grace_L2_input['clm'][:,:] += dt*drift_c[:,:]
		grace_L2_input['slm'][:,:] += dt*drift_s[:,:]

	#-- Correct Pole Tide following Wahr et al. (2015) 10.1002/2015JB011986
	if POLE_TIDE and (DSET == 'GSM'):
		#-- time since 2000.0
		dt = (grace_L2_input['time']-2000.0)
		#-- CSR and JPL Pole Tide Correction
		if PRC in ('UTCSR','JPLEM','JPLMSC'):
			#-- values for IERS mean pole [2010]
			if (grace_L2_input['time'] < 2010.0):
				a = np.array([0.055974,1.8243e-3,1.8413e-4,7.024e-6])
				b = np.array([-0.346346,-1.7896e-3,1.0729e-4,0.908e-6])
			elif (grace_L2_input['time'] >= 2010.0):
				a = np.array([0.023513,7.6141e-3,0.0,0.0])
				b = np.array([-0.358891,0.6287e-3,0.0,0.0])
			#-- calculate m1 and m2 values
			m1 = np.copy(a[0])
			m2 = np.copy(b[0])
			for x in range(1,4):
				m1 += a[x]*dt**x
				m2 += b[x]*dt**x
			#-- pole tide values for CSR and JPL
			#-- CSR and JPL both remove the IERS mean pole from m1 and m2
			#-- before computing their harmonic solutions
			C21_PT = -1.551e-9*(m1 - 0.62e-3*dt) - 0.012e-9*(m2 + 3.48e-3*dt)
			S21_PT = 0.021e-9*(m1 - 0.62e-3*dt) - 1.505e-9*(m2 + 3.48e-3*dt)
			#-- correct GRACE spherical harmonics for pole tide
			#-- note: -= means grace_xlm = grace_xlm - PT
			grace_L2_input['clm'][2,1] -= C21_PT
			grace_L2_input['slm'][2,1] -= S21_PT
		#-- GFZ Pole Tide Correction
		elif ((PRC == 'EIGEN') and (DSET == 'GSM')):
			#-- pole tide values for GFZ
			#-- GFZ removes only a constant pole position
			C21_PT = -1.551e-9*(-0.62e-3*dt) - 0.012e-9*(3.48e-3*dt)
			S21_PT = 0.021e-9*(-0.62e-3*dt) - 1.505e-9*(3.48e-3*dt)
			#-- correct GRACE spherical harmonics for pole tide
			#-- note: -= means grace_xlm = grace_xlm - PT
			grace_L2_input['clm'][2,1] -= C21_PT
			grace_L2_input['slm'][2,1] -= S21_PT

	#-- return the GRACE data, GRACE date (mid-month in decimal), and the
	#-- start and end days as Julian dates
	return grace_L2_input
