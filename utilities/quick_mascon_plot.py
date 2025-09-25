#!/usr/bin/env python
u"""
quick_mascon_plot.py
Written by Tyler Sutterley (06/2024)
Plots a mascon time series file for a particular format

COMMAND LINE OPTIONS:
    -h, --help: Lists the command line options
    -I, --individual: Create individual plots or combine into single
    -O X, --output-file X: Output figure file
    -f X, --figure-format X: Output figure format
    -d X, --figure-dpi X: Output figure resolution in dots per inch
    -H X, --header X: Number of rows of header text to skip
    -m X, --marker X: Plot markers
    -z X, --zorder X: Drawing order of time series
    -t X, --title X: Plot title
    -T X, --time X: Time label for x-axis
    -U X, --units X: Units label for y-axis
    -L X, --legend X: Legend labels for each time series
    -E, --error: Plot mascon errors
    -A, --all: Plot all data without data gap
    -M X, --mode X: Permissions mode of output files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    matplotlib: Python 2D plotting library
        http://matplotlib.org/
        https://github.com/matplotlib/matplotlib

PROGRAM DEPENDENCIES:
    time_series.regress.py: calculates trend coefficients using least-squares

UPDATE HISTORY:
    Updated 06/2024: use wrapper to importlib for optional dependencies
    Updated 08/2023: add option for changing drawing order of time series
    Updated 05/2023: use pathlib to define and operate on paths
        use fit module for getting tidal aliasing terms
    Updated 04/2023: add option to include plot markers on time series plots
    Updated 03/2023: place matplotlib import behind try/except statement
    Updated 01/2023: refactored time series analysis functions
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 05/2022: use argparse descriptions within documentation
    Updated 01/2022: added options for saving files in particular formats
        added option for plotting all data (if not monthly GRACE/GRACE-FO)
    Updated 12/2021: added legend option for labeling each time series
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 10/2020: use argparse to set command line parameters
    Updated 04/2020: added option to skip header text for plotting public data
    Updated 02/2020: made mascon errors optional
    Updated 11/2019: using getopt to set parameters.  use figure tight layout
    Written 10/2019
"""
import sys
import pathlib
import argparse
import numpy as np
import gravity_toolkit as gravtk

# attempt imports
plt = gravtk.utilities.import_dependency('matplotlib.pyplot')

# PURPOSE: read mascon time series file and create plot
def run_plot(i, input_file,
        individual=False,
        header=0,
        marker=None,
        zorder=None,
        title=None,
        time=None,
        units=None,
        legend=None,
        error=False,
        monthly=True
    ):
    """
    Plots a mascon time series file for a particular format

    Arguments
    ---------
    i: iteration of mascon data file
    input_file: mascon time series file

    Keyword arguments
    -----------------
    individual: Create individual plots or combine into single
    header: Number of rows of header text to skip
    marker: Plot marker
    zorder: Drawing order
    title: Region label for title
    time: Time label for x-axis
    units: Units label for y-axis
    legend: Legend labels for each time series
    error: Plot mascon errors
    monthly: Plot monthly data with data gap
    """
    # if creating individual plots
    if individual:
        plt.figure(i+1)
    # read data
    input_file = pathlib.Path(input_file).expanduser().absolute()
    print(input_file.name)
    dinput = np.loadtxt(input_file, skiprows=header)
    # calculate regression
    TERMS = gravtk.time_series.aliasing_terms(dinput[:,1])
    x1 = gravtk.time_series.regress(dinput[:,1],dinput[:,2],ORDER=1,
        CYCLES=[0.5,1.0], TERMS=TERMS)
    x2 = gravtk.time_series.regress(dinput[:,1],dinput[:,2],ORDER=2,
        CYCLES=[0.5,1.0], TERMS=TERMS)
    # print regression coefficients
    args = ('x1',x1['beta'][1],x1['error'][1])
    print('{0}: {1:0.4f} +/- {2:0.4f}'.format(*args))
    args = ('x2',2.0*x2['beta'][2],2.0*x2['error'][2])
    print('{0}: {1:0.4f} +/- {2:0.4f}'.format(*args))
    args = ('AIC',x2['AIC']-x1['AIC'])
    print('{0}: {1:0.4f}'.format(*args))
    # plot all data or monthly data with gap
    if monthly:
        # calculate months and remove GAP from missing
        START_MON,END_MON = (dinput[0,0],dinput[-1,0])
        all_months = np.arange(START_MON,END_MON+1,dtype=np.int64)
        GAP = [187,188,189,190,191,192,193,194,195,196,197]
        MISSING = sorted(set(all_months) - set(dinput[:,0]) - set(GAP))
        months = sorted(set(all_months) - set(MISSING))
        # create a time series with nans for missing months
        tdec = np.full_like(months,np.nan,dtype=np.float64)
        data = np.full_like(months,np.nan,dtype=np.float64)
        if error:
            err = np.full_like(months,np.nan,dtype=np.float64)
        for t,m in enumerate(months):
            valid = np.count_nonzero(dinput[:,0] == m)
            if valid:
                mm, = np.nonzero(dinput[:,0] == m)
                tdec[t] = dinput[mm,1]
                data[t] = dinput[mm,2] - dinput[0,2]
                if error:
                    err[t] = dinput[mm,3]
    else:
        tdec = np.copy(dinput[:,1])
        data = np.copy(dinput[:,2]) - dinput[0,2]
        if error:
            err = np.copy(dinput[:,3])
    # plot all dates
    l, = plt.plot(tdec, data,
        label=legend,
        marker=marker,
        markersize=5,
        zorder=zorder
    )
    # add estimated errors
    if error:
        plt.fill_between(tdec, data-err, y2=data+err,
            color=l.get_color(), alpha=0.25, zorder=zorder)
    # vertical lines for end of the GRACE mission and start of GRACE-FO
    if monthly & ((i == 0) | individual):
        jj, = np.flatnonzero(dinput[:,0] == 186)
        kk, = np.flatnonzero(dinput[:,0] == 198)
        vs = plt.gca().axvspan(dinput[jj,1],dinput[kk,1],
            color='0.5',ls='dashed',alpha=0.15)
        vs._dashes = (4,3)
    # if on the first axes or creating individual plots
    if (i == 0) | individual:
        # add labels
        plt.xlabel(time)
        plt.ylabel(units)
        plt.title(title)
        # use a tight layout to minimize whitespace
        plt.tight_layout()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Plots a mascon time series file for a particular format
            """
    )
    # command line parameters
    parser.add_argument('file',
        type=pathlib.Path, nargs='+',
        help='Mascon data files')
    parser.add_argument('--individual','-I',
        default=False, action='store_true',
        help='Create individual plots or combine into single')
    # output filename, format and dpi
    parser.add_argument('--output-file','-O',
        type=pathlib.Path,
        help='Output figure file')
    parser.add_argument('--figure-format','-f',
        type=str, default='pdf', choices=('pdf','png','jpg','svg'),
        help='Output figure format')
    parser.add_argument('--figure-dpi','-d',
        type=int, default=180,
        help='Output figure resolution in dots per inch (dpi)')
    parser.add_argument('--header','-H',
        type=int, default=0,
        help='Number of rows of header text to skip')
    parser.add_argument('--marker','-m',
        type=str, help='Plot marker')
    parser.add_argument('--zorder','-z',
        type=int, nargs='+',
        help='Drawing order for each time series')
    parser.add_argument('--title','-t',
        type=lambda x: ' '.join(str.split(x,"_")),
        help='Plot title')
    parser.add_argument('--time','-T',
        type=str, default='Time [Yr]',
        help='Time label for x-axis')
    parser.add_argument('--units','-U',
        type=str, default='Mass [Gt]',
        help='Units label for y-axis')
    parser.add_argument('--legend','-L',
        type=str, nargs='+',
        help='Legend labels for each time series')
    parser.add_argument('--error','-E',
        default=False, action='store_true',
        help='Plot mascon errors')
    parser.add_argument('--all','-A',
        default=True, action='store_false',
        help='Plot all data without data gap')
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of output files')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # run plot program for each input file
    for i,f in enumerate(args.file):
        legend = args.legend[i] if args.legend else None
        zorder = args.zorder[i] if args.zorder else None
        run_plot(i, f, individual=args.individual, header=args.header,
            marker=args.marker, title=args.title, time=args.time,
            units=args.units, legend=legend, error=args.error,
            monthly=args.all, zorder=zorder)
    # add legend if applicable
    if args.legend:
        lgd = plt.legend(loc=3,frameon=False)
        lgd.get_frame().set_alpha(1.0)
        for line in lgd.get_lines():
            line.set_linewidth(6)
    # show or save plots
    if args.output_file and (len(plt.get_fignums()) > 1):
        # for each figure number
        for num in plt.get_fignums():
            fig = plt.figure(num)
            output = f'{args.output_file.stem}_{num}{args.output_file.suffix}'
            # save the figure file
            fig.savefig(output,
                metadata={'Title':pathlib.Path(sys.argv[0]).name},
                dpi=args.figure_dpi,
                format=args.figure_format
            )
            # change the permissions mode
            output.chmod(mode=args.mode)
    elif args.output_file:
        # save the figure file
        plt.savefig(args.output_file,
            metadata={'Title':pathlib.Path(sys.argv[0]).name},
            dpi=args.figure_dpi,
            format=args.figure_format
        )
        # change the permissions mode
        args.output_file.chmod(mode=args.mode)
    else:
        # show all plots
        plt.show()
    # close all figure axes and windows
    plt.cla()
    plt.clf()
    plt.close()

# run main program
if __name__ == '__main__':
    main()
