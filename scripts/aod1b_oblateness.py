#!/usr/bin/env python
u"""
aod1b_oblateness.py
Written by Tyler Sutterley (07/2020)

Reads GRACE/GRACE-FO level-1b dealiasing data files for a specific product
    atm: atmospheric loading from ECMWF
    ocn: oceanic loading from OMCT
    glo: global atmospheric and oceanic loading
    oba: ocean bottom pressure from OMCT

Creates monthly files of oblateness (C20) variations at 6-hour intervals

NOTE: this reads the GFZ AOD1B files downloaded from PO.DAAC
https://podaac-uat.jpl.nasa.gov/drive/files/allData/grace/L1B/GFZ/AOD1B/RL05/

CALLING SEQUENCE:
    aod1b_oblateness(base_dir, DREL='RL06', DSET='glo', CLOBBER=True)

INPUTS:
    base_dir: working data directory

OPTIONS:
    DREL: GRACE/GRACE-FO data release (RL05 or RL06)
    DSET: GRACE/GRACE-FO dataset (atm, ocn, glo, oba)
    CLOBBER: overwrite existing data
    MODE: Permission mode of directories and files
    VERBOSE: Output information for each output file

COMMAND LINE OPTIONS:
    -D X, --directory=X: Working Data Directory
    -R X, --release=X: GRACE Data Release (RL05 or RL06)
    -C, --clobber: Overwrite existing data
    -M X, --mode=X: Permission mode of directories and files
    -V, --verbose: Output information for each output file

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

UPDATED HISTORY:
    Updated 07/2020: added function docstrings
    Updated 06/2019: using python3 compatible regular expression patterns
    Updated 10/2018: using future division for python3 Compatibility
    Updated 08/2018: using full release string (RL05 instead of 5)
    Updated 06/2018: can read RL06 AOD1B files from PO.DAAC
    Updated 03/2018: can read tar files from GFZ_AOD1b_sync.py
    Updated 04/2017: slight modifications to the regular expression patterns
        to verify that the suffix is the end of a given filename (no .xml files)
    Updated 02/2017: using getopt to set parameters and data directory
        do not extract tar files to temp, extract contents of files to memory
    Updated 05-06/2016: oba=ocean bottom pressure, absolute import of shutil
    Written 05/2016
"""
from __future__ import print_function, division

import sys
import os
import re
import gzip
import getopt
import tarfile
import numpy as np

#-- aod1b data products
product = {}
product['atm'] = 'Atmospheric loading from ECMWF'
product['ocn'] = 'Oceanic loading from OMCT'
product['glo'] = 'Global atmospheric and oceanic loading'
product['oba'] = 'Ocean bottom pressure from OMCT'

#-- program module to read the C20 coefficients of the AOD1b data
def aod1b_oblateness(base_dir, DREL='', DSET='', CLOBBER=False, MODE=0o775,
    VERBOSE=False):
    """
    Creates monthly files of oblateness (C20) variations at 6-hour intervals
    from GRACE/GRACE-FO level-1b dealiasing data files

    Arguments
    ---------
    base_dir: working data directory

    Keyword arguments
    -----------------
    DREL: GRACE/GRACE-FO data release
    DSET: GRACE/GRACE-FO dataset
        atm: atmospheric loading from ECMWF
        ocn: oceanic loading from OMCT/MPIOM
        glo: global atmospheric and oceanic loading
        oba: ocean bottom pressure from OMCT/MPIOM
    CLOBBER: overwrite existing data
    MODE: Permission mode of directories and files
    VERBOSE: Output information for each output file
    """

    #-- compile regular expressions operators for file dates
    #-- will extract the year and month from the tar file (.tar.gz)
    tx = re.compile(r'AOD1B_(\d+)-(\d+)_\d+\.(tar\.gz|tgz)$', re.VERBOSE)
    #-- and the calendar day from the ascii file (.asc or gzipped .asc.gz)
    fx = re.compile(r'AOD1B_\d+-\d+-(\d+)_X_\d+.asc(.gz)?$', re.VERBOSE)
    #-- compile regular expressions operator for the clm/slm headers
    #-- for the specific AOD1b product
    hx = re.compile(r'^DATA.*SET.*{0}'.format(DSET), re.VERBOSE)
    #-- compile regular expression operator to find numerical instances
    #-- will extract the data from the file
    regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?'
    rx = re.compile(regex_pattern, re.VERBOSE)
    #-- output formatting string
    fstr = '{0:4d}-{1:02d}-{2:02d}T{3:02d}:00:00 {4:+16.8E}'

    #-- Maximum spherical harmonic degree (LMAX)
    LMAX = 100
    #-- Calculating the number of cos and sin harmonics up to LMAX
    n_harm = (LMAX**2 + 3*LMAX)//2 + 1

    #-- AOD1B directory and output oblateness directory
    grace_dir = os.path.join(base_dir,'AOD1B',DREL)
    output_dir = os.path.join(grace_dir,'oblateness')
    if not os.access(output_dir, os.F_OK):
        os.mkdir(output_dir, MODE)

    #-- finding all of the tar files in the AOD1b directory
    input_tar_files = [tf for tf in os.listdir(grace_dir) if tx.match(tf)]

    #-- for each tar file
    for i in sorted(input_tar_files):
        #-- extract the year and month from the file
        YY,MM,SFX = tx.findall(i).pop()
        YY,MM = np.array([YY,MM], dtype=np.int)
        #-- output monthly oblateness file
        FILE = 'AOD1B_{0}_{1}_{2:4d}_{3:02d}.txt'.format(DREL,DSET,YY,MM)
        #-- if output file exists: check if input tar file is newer
        TEST = False
        OVERWRITE = ' (clobber)'
        #-- check if output file exists
        if os.access(os.path.join(output_dir,FILE), os.F_OK):
            #-- check last modification time of input and output files
            input_mtime = os.stat(os.path.join(grace_dir,i)).st_mtime
            output_mtime = os.stat(os.path.join(output_dir,FILE)).st_mtime
            #-- if input tar file is newer: overwrite the output file
            if (input_mtime > output_mtime):
                TEST = True
                OVERWRITE = ' (overwrite)'
        else:
            TEST = True
            OVERWRITE = ' (new)'
        #-- As there are so many files.. this will only read the new files
        #-- or will rewrite if CLOBBER is set (if wanting something changed)
        if TEST or CLOBBER:
            #-- if verbose: output information about the oblateness file
            if VERBOSE:
                print('{0}{1}'.format(os.path.join(output_dir,FILE),OVERWRITE))
            #-- open output monthly oblateness file
            f = open(os.path.join(output_dir,FILE), 'w')
            args = ('Oblateness time series',DREL,DSET)
            print('# {0} from {1} AOD1b {2} Product'.format(*args), file=f)
            print('# {0}'.format(product[DSET]), file=f)
            print('# {0:^15}    {1:^15}'.format('ISO-Time','C20'), file=f)

            #-- open the AOD1B monthly tar file
            tar = tarfile.open(name=os.path.join(grace_dir,i), mode='r:gz')

            #-- Iterate over every member within the tar file
            for member in tar.getmembers():
                #-- get calendar day from file
                DD,SFX = fx.findall(member.name).pop()
                DD = np.int(DD)
                #-- open datafile for day
                if (SFX == '.gz'):
                    fid = gzip.GzipFile(fileobj=tar.extractfile(member))
                else:
                    fid = tar.extractfile(member)
                #-- C20 spherical harmonics for day and hours
                C20 = np.zeros((4))
                hours = np.zeros((4),dtype=np.int)

                #-- create counter for hour in dataset
                c = 0
                #-- while loop ends when dataset is read
                while (c < 4):
                    #-- read line
                    file_contents = fid.readline().decode('ISO-8859-1')
                    #-- find file header for data product
                    if bool(hx.search(file_contents)):
                        #-- extract hour from header and convert to float
                        HH, = re.findall(r'(\d+):\d+:\d+',file_contents)
                        hours[c] = np.int(HH)
                        #-- read each line of spherical harmonics
                        for k in range(0,n_harm):
                            file_contents = fid.readline().decode('ISO-8859-1')
                            #-- find numerical instances in the data line
                            line_contents = rx.findall(file_contents)
                            #-- spherical harmonic degree and order
                            l1 = np.int(line_contents[0])
                            m1 = np.int(line_contents[1])
                            if (l1 == 2) and (m1 == 0):
                                C20[c] = np.float(line_contents[2])
                        #-- add 1 to hour counter
                        c += 1
                #-- close the input file for day
                fid.close()
                #-- write to file for each hour
                for h in range(4):
                    print(fstr.format(YY,MM,DD,hours[h],C20[h]),file=f)

            #-- close the tar file
            tar.close()
            #-- close the output file
            f.close()
            #-- set the permissions mode of the output file
            os.chmod(os.path.join(output_dir,FILE), MODE)

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
    print(' -D X, --directory=X\tWorking Data Directory')
    print(' -R X, --release=X\tGRACE Data Release')
    print(' -C, --clobber\t\tOverwrite existing data')
    print(' -M X, --mode=X\t\tPermission mode of directories and files')
    print(' -V, --verbose\t\tOutput information for each output file\n')

#-- Main program that calls aod1b_oblateness()
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','directory=','release=','clobber','mode=','verbose']
    optlist,arglist = getopt.getopt(sys.argv[1:],'hD:R:CM:V',long_options)

    #-- working data directory
    base_dir = os.getcwd()
    #-- Data release
    DREL = 'RL06'
    #-- clobber will overwrite the existing data
    CLOBBER = False
    #-- permissions mode of the local directories and files (number in octal)
    MODE = 0o775
    #-- verbose will output information about each output file
    VERBOSE = False
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("-D","--directory"):
            base_dir = os.path.expanduser(arg)
        elif opt in ("-R","--release"):
            DREL = arg
        elif opt in ("-C","--clobber"):
            CLOBBER = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)
        elif opt in ("-V","--verbose"):
            VERBOSE = True

    #-- check that DSET was entered as system argument
    if not arglist:
        raise Exception('No System Arguments Listed')

    for DSET in arglist:
        #-- verify case
        DSET = DSET.lower() if DSET.isupper() else DSET
        #-- check that DSET was correctly entered
        if DSET not in product.keys():
            for key,val in product.items():
                print('{0}: {1}'.format(key, val))
            raise ValueError('Incorrect Data Product Entered (atm,glo,ocn,oba)')
        #-- run AOD1b oblateness program with parameters
        aod1b_oblateness(base_dir,DREL=DREL,DSET=DSET,CLOBBER=CLOBBER,MODE=MODE,
            VERBOSE=VERBOSE)

#-- run main program
if __name__ == '__main__':
    main()
