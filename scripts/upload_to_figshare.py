#!/usr/bin/env python
u"""
upload_to_figshare.py
Written by Tyler Sutterley (05/2023)

Uploads geocenter files to figshare repository using secure ftp

CALLING SEQUENCE:
    python upload_to_figshare.py --directory <path_to_geocenter_files>

COMMAND LINE OPTIONS:
    -D X, --directory X: working data directory with geocenter files
    -U X, --user X: username for figshare ftp login
    -W X, --password X: password for figshare ftp login
    -r X, --release X: GRACE/GRACE-FO data release
    -t X, --timeout X: Timeout in seconds for blocking operations
    -V, --verbose: increase ftp verbosity level

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 05/2022: use argparse descriptions within documentation
    Written 05/2021
"""
import os
import netrc
import getpass
import pathlib
import argparse
import builtins
import gravity_toolkit as gravtk

# PURPOSE: upload geocenter files to figshare
def upload_to_figshare(grace_dir,DREL,username=None,password=None,
    timeout=None,verbose=False):
    # directory setup
    grace_dir = pathlib.Path(grace_dir).expanduser().absolute()
    # figshare directory for geocenter products
    directory = ('Geocenter Estimates from Time-Variable Gravity '
        'and Ocean Model Outputs')
    # labels for each processing center
    input_flags = {}
    input_flags['CSR'] = ['SLF_iter','SLF_iter_wAOD']
    input_flags['GFZ'] = ['SLF_iter','SLF_iter_wAOD','SLF_iter_wSLR21']
    input_flags['JPL'] = ['SLF_iter','SLF_iter_wAOD']
    # ocean model labels
    model_str = 'OMCT' if DREL in ('RL04','RL05') else 'MPIOM'
    # build list of geocenter files
    geocenter_files = []
    # for each processing center
    for PROC,center_flags in input_flags.items():
        # for each data product flag
        for flag in center_flags:
            f = '{0}_{1}_{2}_{3}.txt'.format(PROC,DREL,model_str,flag)
            geocenter_file = grace_dir.joinpath(f)
            if not geocenter_file.exists():
                raise FileNotFoundError(f'Geocenter file {f} not found')
            else:
                geocenter_files.append(geocenter_file)
    # upload geocenter files to figshare
    gravtk.utilities.to_figshare(geocenter_files,
        username=username,password=password,timeout=timeout,
        directory=directory,verbose=verbose)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Uploads geocenter files to figshare repository
            using secure ftp
            """
    )
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # figshare credentials
    parser.add_argument('--user','-U',
        type=str, default=os.environ.get('FIGSHARE_FTP_USER'),
        help='Username for Figshare ftp login')
    parser.add_argument('--password','-W',
        type=str, default=os.environ.get('FIGSHARE_PASSWORD'),
        help='Password for Figshare ftp login')
    parser.add_argument('--netrc','-N',
        type=pathlib.Path,
        default=pathlib.Path.home().joinpath('.netrc'),
        help='Path to .netrc file for authentication')
    # GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str,
        default='RL06', choices=['RL04','RL05','RL06'],
        help='GRACE/GRACE-FO data release')
    # connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=None,
        help='Timeout in seconds for blocking operations')
    # verbose will output information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # figshare ftp hostname
    HOST = 'ftps.figshare.com'
    # get figshare ftp credentials
    try:
        args.user,_,args.password = netrc.netrc(args.netrc).authenticators(HOST)
    except:
        # check that figshare ftp credentials were entered
        if not args.user:
            prompt = f'Username for {HOST}: '
            args.user = builtins.input(prompt)
        # enter figshare password securely from command-line
        if not args.password:
            prompt = f'Password for {args.user}@{HOST}: '
            args.password = getpass.getpass(prompt)

    # run program with parameters
    upload_to_figshare(args.directory,args.release,
        username=args.user,password=args.password,
        timeout=args.timeout,verbose=args.verbose)

# run main program
if __name__ == '__main__':
    main()
