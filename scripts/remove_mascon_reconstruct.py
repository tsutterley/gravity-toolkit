#!/usr/bin/env python
u"""
remove_mascon_reconstruct.py
Written by Tyler Sutterley (05/2023)
Removes reconstructed mascon files from the index file after running program

COMMAND LINE OPTIONS:
    --help: list the command line options
    --reconstruct-file X: reconstructed mascon time series file
    -V, --verbose: verbose output of processing run

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within documentation
    Updated 04/2022: include utf-8 encoding in reads to be windows compliant
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 11/2021: using python logging for handling verbose output
    Updated 07/2021: switch from parameter files to argparse arguments
    Updated 04/2021: add parser object for removing commented or empty lines
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: flake8 compatible regular expression strings
    Written 12/2019
"""
from __future__ import print_function

import sys
import os
import re
import logging
import pathlib
import argparse
import traceback
import gravity_toolkit as gravtk

# PURPOSE: keep track of threads
def info(args):
    logging.info(pathlib.Path(sys.argv[0]).name)
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: Remove reconstructed spherical harmonic fields to free space
def remove_mascon_reconstruct(RECONSTRUCT_FILE):
    # file parser for reading index files
    # removes commented lines (can comment out files in the index)
    # removes empty lines (if there are extra empty lines)
    parser = re.compile(r'^(?!\#|\%|$)', re.VERBOSE)
    # output directory setup
    RECONSTRUCT_FILE = pathlib.Path(RECONSTRUCT_FILE).expanduser().absolute()
    with RECONSTRUCT_FILE.open(mode='r', encoding='utf8') as f:
        file_list = [l for l in f.read().splitlines() if parser.match(l)]
    # for each valid file in the index (iterate over mascons)
    for fi in file_list:
        # remove reconstructed spherical harmonics file
        fi = pathlib.Path(fi).expanduser().absolute()
        fi.unlink()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Removes reconstructed mascon files from the index file
            after running program
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # mascon reconstruct parameters
    parser.add_argument('--reconstruct-file',
        type=pathlib.Path,
        help='Reconstructed mascon time series file')
    # print information about processing run
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of processing run')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # create logger
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # try to run the analysis with listed parameters
    try:
        info(args)
        # run remove_mascon_reconstruct algorithm with parameters
        remove_mascon_reconstruct(args.reconstruct_file)
    except Exception as exc:
        # if there has been an error exception
        # print the type, value, and stack trace of the
        # current exception being handled
        logging.critical(f'process id {os.getpid():d} failed')
        logging.error(traceback.format_exc())

# run main program
if __name__ == '__main__':
    main()
