#!/usr/bin/env python
u"""
simple_parallel_shell.py
Written by Tyler Sutterley (05/2023)
Runs a shell script in parallel using python multiprocessing

COMMAND LINE OPTIONS:
    -h, --help: Lists the command line options
    -P X, --np X: Run in parallel with X number of processes
    -E X, --environment X: Mapping of environment variables
    -S, --shell: Run specified commands through the shell
    -V, --verbose: Verbose output of processing runs

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 05/2022: use argparse descriptions within documentation
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: added verbose option to print outputs
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: flake8 compatible regular expression strings
    Updated 11/2019: remove commented lines from list of commands
    Written 10/2019
"""
from __future__ import print_function

import sys
import os
import re
import shlex
import logging
import pathlib
import argparse
import traceback
import subprocess
import multiprocessing

# PURPOSE: converts a command line argument into an environment dictionary
def argtoenv(arg):
    output = {}
    for item in arg.split(','):
        key,val = item.split(':')
        output[key] = pathlib.Path(val).expanduser().absolute()
    return output

# PURPOSE: run command and handle error exceptions
def execute_command(comm, shell=False, env=None, verbose=False):
    # create logger
    loglevel = logging.INFO if verbose else logging.CRITICAL
    logging.basicConfig(level=loglevel)
    # try to run the command
    try:
        # run the command
        args = shlex.split(comm.replace('~', pathlib.Path.home()))
        p = subprocess.Popen(args, stdout=subprocess.PIPE, shell=shell, env=env)
        # Wait for process to finish
        p.wait()
        # print outputs from command
        stat, err = p.communicate()
        logging.info(stat.decode('utf-8'))
    except Exception as exc:
        # if there has been an error exception
        # print the type, value, and stack trace of the
        # current exception being handled
        logging.critical(f'process id {os.getpid():d} failed')
        logging.error(traceback.format_exc())

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Runs a shell script in parallel using python
            multiprocessing
            """
    )
    # command line parameters
    parser.add_argument('files',
        type=pathlib.Path, nargs='+',
        help='Shell script files to run')
    # number of processes to run in parallel
    parser.add_argument('--np','-P',
        metavar='PROCESSES', type=int, default=0,
        help='Number of processes to run in parallel')
    # Mapping of environment variables
    parser.add_argument('--environment','-E',
        type=argtoenv,
        help='Mapping of environment variables')
    # Run specified commands through the shell
    parser.add_argument('--shell','-S',
        default=False, action='store_true',
        help='Run specified commands through the shell')
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

    # create list of commands to run from all input shell script files
    command_list = []
    # for each input shell script file
    for fi in args.files:
        # open the shell script file
        input_shell_file = pathlib.Path(fi).expanduser().absolute()
        with input_shell_file.open(mode='r', encoding='utf-8') as f:
            # connect lines that have a line continuation before splitting
            file_contents = f.read().replace('\r\n','\n').replace('\\\n','')
            command_list.extend(file_contents.splitlines())
    # reduce command list to valid and uncommented lines
    command_list = [l for l in command_list if re.match(r'^(?!#)(.+?)',l)]

    # run in parallel with multiprocessing Pool
    pool = multiprocessing.Pool(processes=args.np)
    # for each command
    for comm in command_list:
        kwds = dict(env=args.environment, shell=args.shell,
            verbose=args.verbose)
        pool.apply_async(execute_command, args=(comm,), kwds=kwds)
    # start multiprocessing jobs
    # close the pool
    # prevents more tasks from being submitted to the pool
    pool.close()
    # exit the completed processes
    pool.join()

# run main program
if __name__ == '__main__':
    main()
