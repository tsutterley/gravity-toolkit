#!/usr/bin/env python
u"""
plot_global_grid_4maps.py
Written by Tyler Sutterley (05/2023)
Creates 4 GMT-like plots in a Plate Carree (Equirectangular) projection

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        http://www.numpy.org
        http://www.scipy.org/NumPy_for_Matlab_Users
    scipy: Scientific Tools for Python
        http://www.scipy.org/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        http://h5py.org
    matplotlib: Python 2D plotting library
        http://matplotlib.org/
        https://github.com/matplotlib/matplotlib
    cartopy: Python package designed for geospatial data processing
        https://scitools.org.uk/cartopy
    pyshp: Python read/write support for ESRI Shapefile format
        https://github.com/GeospatialPython/pyshp

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
        added option to set the input variable names or column order
    Updated 03/2023: switch from parameter files to argparse arguments
        updated inputs to spatial from_ascii function
    Updated 07/2022: place some imports behind try/except statements
    Updated 05/2022: use argparse descriptions within documentation
        use mask, shift grid and interpolation functions from tools
    Updated 12/2021: added color palette table (cpt) file reader from tools
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 10/2020: using spatial utilities for reading
        using argparse to set command-line parameters
    Updated 10/2019: changing Y/N flags to True/False
    Updated 09/2019: added parameter for specifying if netCDF4 or HDF5
    Updated 04/2019: set cap style of cartopy geoaxes outline patch
    Updated 03/2019: replacing matplotlib basemap with cartopy
    Updated 02/2019: saving metadata to output figure
    Updated 12/2018: added parameter CBEXTEND for colorbar extension triangles
    Updated 09/2018: use a custom mask_oceans function instead of default maskoceans
        added REMOVE parameter to subtract a spatial field from each input map
    Updated 05/2018: using __future__ print function for python3 compatibility
    Updated 11/2017: can plot a contour of the global average with MEAN
    Written 08/2017
"""
from __future__ import print_function

import sys
import os
import copy
import logging
import pathlib
import argparse
import warnings
import traceback
import numpy as np
import scipy.ndimage
import scipy.interpolate
import gravity_toolkit as gravtk

# attempt imports
try:
    import cartopy.crs as ccrs
except ModuleNotFoundError:
    warnings.warn("cartopy not available", ImportWarning)
try:
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    import matplotlib.ticker as ticker
    import matplotlib.offsetbox as offsetbox
    matplotlib.rcParams['axes.linewidth'] = 2.0
    matplotlib.rcParams['font.family'] = 'sans-serif'
    matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
    matplotlib.rcParams['mathtext.default'] = 'regular'
except ModuleNotFoundError:
    warnings.warn("matplotlib not available", ImportWarning)
try:
    import shapefile
except ModuleNotFoundError:
    warnings.warn("shapefile not available", ImportWarning)

# cartopy transform for Equirectangular Projection
try:
    projection = ccrs.PlateCarree()
except (NameError,ValueError) as exc:
    pass

# PURPOSE: keep track of threads
def info(args):
    logging.info(pathlib.Path(sys.argv[0]).name)
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE plot coastlines and islands (GSHHS with G250 Greenland)
def plot_coastline(ax, base_dir, LINEWIDTH=0.5):
    # read the coastline shape file
    coastline_dir = base_dir.joinpath('masks','G250')
    coastline_shape_files = []
    coastline_shape_files.append('GSHHS_i_L1_no_greenland.shp')
    coastline_shape_files.append('greenland_coastline_islands.shp')
    for fi,S in zip(coastline_shape_files,[1000,200]):
        coast_shapefile = coastline_dir.joinpath(fi)
        logging.debug(str(coast_shapefile))
        shape_input = shapefile.Reader(str(coast_shapefile))
        shape_entities = shape_input.shapes()
        # for each entity within the shapefile
        for c,ent in enumerate(shape_entities[:S]):
            # extract coordinates and plot
            lon,lat = np.transpose(ent.points)
            ax.plot(lon, lat, c='k', lw=LINEWIDTH, transform=projection)

# PURPOSE: plot Antarctic grounded ice delineation
def plot_grounded_ice(ax, base_dir, LINEWIDTH=0.5):
    grounded_ice_file = ['masks','IceBoundaries_Antarctica_v02',
        'ant_ice_sheet_islands_v2.shp']
    grounded_ice_shapefile = base_dir.joinpath(*grounded_ice_file)
    logging.debug(str(grounded_ice_shapefile))
    shape_input = shapefile.Reader(str(grounded_ice_shapefile))
    shape_entities = shape_input.shapes()
    shape_attributes = shape_input.records()
    i = [i for i,e in enumerate(shape_entities) if (np.ndim(e.points) > 1)]
    # cartopy transform for NSIDC polar stereographic south
    projection = ccrs.Stereographic(central_longitude=0.0,
        central_latitude=-90.0,true_scale_latitude=-71.0)
    for indice in i:
        # extract Polar-Stereographic coordinates for record
        pts = np.array(shape_entities[indice].points)
        ax.plot(pts[:,0], pts[:,1], c='k', lw=LINEWIDTH, transform=projection)

# plot grid program
def plot_grid(base_dir, FILENAMES,
    DATAFORM=None,
    VARIABLES=[],
    MASK=None,
    INTERPOLATION=None,
    DDEG=None,
    INTERVAL=None,
    SCALE_FACTOR=1.0,
    REMOVE_FILE=None,
    COLOR_MAP=None,
    CPT_FILE=None,
    PLOT_RANGE=None,
    BOUNDARY=None,
    ALPHA=1.0,
    CONTOURS=False,
    CONTOUR_RANGE=None,
    MEAN_CONTOUR=False,
    TITLES=None,
    LABELS=None,
    CBEXTEND=None,
    CBTITLE=None,
    CBUNITS=None,
    CBFORMAT=None,
    DRAW_GRID_LINES=False,
    GRID=None,
    FIGURE_FILE=None,
    FIGURE_FORMAT=None,
    FIGURE_DPI=None,
    MODE=0o775):

    # extend list if a single format was entered for all files
    if len(DATAFORM) < len(FILENAMES):
        DATAFORM = DATAFORM*len(FILENAMES)

    # read CPT or use color map
    if CPT_FILE is not None:
        # cpt file
        cmap = gravtk.tools.from_cpt(CPT_FILE)
    else:
        # colormap
        cmap = copy.copy(cm.get_cmap(COLOR_MAP))
    # grey color map for bad values
    cmap.set_bad('lightgray',1.0)

    # set transparency ALPHA
    if BOUNDARY is None:
        # contours
        levels = np.arange(PLOT_RANGE[0], PLOT_RANGE[1]+PLOT_RANGE[2],
            PLOT_RANGE[2])
        norm = colors.Normalize(vmin=PLOT_RANGE[0], vmax=PLOT_RANGE[1])
    else:
        # boundary between contours
        levels = np.array(BOUNDARY, dtype=np.float64)
        norm = colors.BoundaryNorm(BOUNDARY, ncolors=256)

    # convert degree spacing and interval parameters
    # Grid spacing
    dlon,dlat = (DDEG[0],DDEG[0]) if (len(DDEG) == 1) else (DDEG[0],DDEG[1])
    # Grid dimensions
    if (INTERVAL == 1):# (0:360, 90:-90)
        nlon = np.int64((360.0/dlon)+1.0)
        nlat = np.int64((180.0/dlat)+1.0)
    elif (INTERVAL == 2):# degree spacing/2
        nlon = np.int64((360.0/dlon))
        nlat = np.int64((180.0/dlat))

    # interpolation method for image background using transform_scalar
    if (INTERPOLATION == 'nearest'):
        order = 0
    elif (INTERPOLATION == 'bilinear'):
        order = 1
    elif (INTERPOLATION == 'cubic'):
        order = 3

    # create masked array if missing values
    if MASK is not None:
        # Read Land-Sea Mask of specified input file
        # 0=Ocean, 1=Land, 2=Lake, 3=Small Island, 4=Ice Shelf
        # Open the land-sea NetCDF file for reading
        landsea = gravtk.spatial().from_netCDF4(MASK,
            date=False, varname='LSMASK')
        # create land function
        nth,nphi = landsea.shape
        mask = np.zeros((nth,nphi),dtype=bool)
        # combine land and island levels for land function
        indx,indy = np.nonzero((landsea.data >= 1) & (landsea.data <= 3))
        mask[indx,indy] = True

    # remove a spatial field from each input map
    if REMOVE_FILE is not None:
        REMOVE = gravtk.spatial().from_netCDF4(REMOVE_FILE,
            date=False).data[:,:]
    else:
        REMOVE = 0.0

    # image extents
    ax = {}
    # setup Plate Carree projection
    fig, ((ax[0],ax[1]),(ax[2],ax[3])) = plt.subplots(num=1, nrows=2, ncols=2,
        figsize=(10.375,7.0), subplot_kw=dict(projection=projection))
    # WGS84 Ellipsoid parameters
    a_axis = 6378137.0# [m] semimajor axis of the ellipsoid
    flat = 1.0/298.257223563# flattening of the ellipsoid
    # (4pi/3)R^3 = (4pi/3)(a^2)b = (4pi/3)(a^3)(1 -f)
    rad_e = a_axis*(1.0 - flat)**(1.0/3.0)

    for i,ax1 in ax.items():

        # input ascii/netCDF4/HDF5 file
        if (DATAFORM[i] == 'ascii'):
            # ascii (.txt)
            dinput = gravtk.spatial().from_ascii(FILENAMES[i], date=False,
                columns=VARIABLES, spacing=[dlon,dlat], nlat=nlat, nlon=nlon)
        elif (DATAFORM[i] == 'netCDF4'):
            # netCDF4 (.nc)
            field_mapping = gravtk.spatial().default_field_mapping(VARIABLES)
            dinput = gravtk.spatial().from_netCDF4(FILENAMES[i], date=False,
                field_mapping=field_mapping)
        elif (DATAFORM[i] == 'HDF5'):
            # HDF5 (.H5)
            field_mapping = gravtk.spatial().default_field_mapping(VARIABLES)
            dinput = gravtk.spatial().from_HDF5(FILENAMES[i], date=False,
                field_mapping=field_mapping)

        # remove offset and scale to units
        if (REMOVE != 0.0) or (SCALE_FACTOR != 1.0):
            dinput = dinput.offset(-REMOVE).scale(SCALE_FACTOR)

        # update masked values
        if MASK is not None:
            dinput.replace_invalid(fill_value=dinput.fill_value, mask=mask)

        # if dlat is negative
        if (np.sign(dlat) == -1):
            dinput = dinput.flip(axis=0)

        # calculate image coordinates
        xmin,xmax,ymin,ymax = ax1.get_extent()
        mx = np.int64((xmax-xmin)/0.5)+1
        my = np.int64((ymax-ymin)/0.5)+1
        X = np.linspace(xmin,xmax,mx)
        Y = np.linspace(ymin,ymax,my)
        gridx,gridy = np.meshgrid(X,Y)
        # create mesh lon/lat
        points = projection.transform_points(projection,
            gridx.flatten(), gridy.flatten())
        lonsin = points[:,0].reshape(my,mx)
        latsin = points[:,1].reshape(my,mx)

        # interpolate to image coordinates
        if (INTERVAL == 1) and (np.max(dinput.lon) > 180):# (0:360, 90:-90)
            shift_data,lon180 = gravtk.tools.shift_grid(180.0,
                dinput.data,dinput.lon)
            shift_mask,lon180 = gravtk.tools.shift_grid(180.0,
                dinput.mask,dinput.lon)
            img = gravtk.tools.interp_grid(shift_data,lon180,
                dinput.lat,lonsin,latsin,order)
            msk = gravtk.tools.interp_grid(shift_mask,lon180,
                dinput.lat,lonsin,latsin,order)
        elif (INTERVAL == 2) and (np.max(dinput.lon) > 180):# DDEG/2
            shift_data,lon180 = gravtk.tools.shift_grid(180.0+dlon,
                dinput.data,dinput.lon)
            shift_mask,lon180 = gravtk.tools.shift_grid(180.0+dlon,
                dinput.mask,dinput.lon)
            img = gravtk.tools.interp_grid(shift_data,
                lon180,dinput.lat,lonsin,latsin,order)
            msk = gravtk.tools.interp_grid(shift_mask,
                lon180,dinput.lat,lonsin,latsin,order)
        else:# -180:180 or modification of there of
            img = gravtk.tools.interp_grid(dinput.data,
                dinput.lon,dinput.lat,lonsin,latsin,order)
            msk = gravtk.tools.interp_grid(dinput.mask,
                dinput.lon,dinput.lat,lonsin,latsin,order)
        # create masked array of image
        img = np.ma.array(img, mask=msk.astype(bool))

        # plot only grounded points
        if MASK is not None:
            img = gravtk.tools.mask_oceans(lonsin,latsin,
                data=img,order=order,iceshelves=False)
        # plot image with transparency using normalization
        im = ax1.imshow(img, interpolation='nearest',
            extent=(xmin,xmax,ymin,ymax),
            cmap=cmap, norm=norm, alpha=ALPHA,
            origin='lower', transform=projection)

        # create mesh lon/lat
        lon, lat = np.meshgrid(dinput.lon,dinput.lat)
        # recalculate data at zoomed coordinates
        data = np.ma.array(scipy.ndimage.zoom(img.data,5,order=1))
        mask = scipy.ndimage.zoom(np.invert(img.mask),5,order=1,output=bool)
        data.mask = np.invert(mask)
        # plot line contours
        if CONTOURS:
            clevs = np.arange(CONTOUR_RANGE[0],
                CONTOUR_RANGE[1] + CONTOUR_RANGE[2],
                CONTOUR_RANGE[2])
            # remove 0 (will plot in red)
            reduce_clevs = clevs[np.nonzero(clevs)]
            # plot contours
            ax1.contour(data,reduce_clevs,colors='0.2',linestyles='solid',
                extent=(xmin,xmax,ymin,ymax),origin='lower',
                transform=projection)
            ax1.contour(data,[0],colors='red',linestyles='solid',linewidths=1.5,
                extent=(xmin,xmax,ymin,ymax),origin='lower',
                transform=projection)

        # plot line contour for global average
        if MEAN_CONTOUR and CONTOURS:
            # calculate areas of each grid cell
            dphi,dth = (dlon*np.pi/180.0,dlat*np.pi/180.0)
            indy,indx = np.nonzero(np.logical_not(dinput.mask))
            area = (rad_e**2)*dth*dphi*np.cos(lat[indy,indx]*np.pi/180.0)
            # calculate average
            ave = np.sum(area*dinput.data[indy,indx])/np.sum(area)
            # plot line contour of global average
            ax1.contour(data,[ave],colors='blue',linestyles='solid',linewidths=1.5,
                extent=(xmin,xmax,ymin,ymax),origin='lower',
                transform=projection)

        # draw coastlines
        plot_coastline(ax1, base_dir)
        # plot grounded ice exterior for Antarctica
        plot_grounded_ice(ax1, base_dir)

        # draw lat/lon grid lines
        if DRAW_GRID_LINES:
            # meridian and parallel grid spacing
            llx,lly = (GRID[0],GRID[0]) if (len(GRID) == 1) else (GRID[0],GRID[1])
            grid_meridians = np.arange(-180, 180 + llx, llx)
            grid_parallels = np.arange(-90, 90 + lly, lly)
            gl = ax1.gridlines(crs=projection, draw_labels=False,
                linewidth=0.1, color='0.25', linestyle='-')
            gl.xlocator = ticker.FixedLocator(grid_meridians)
            gl.ylocator = ticker.FixedLocator(grid_parallels)

        # add title for each subplot
        if TITLES is not None:
            TITLE = ' '.join(TITLES[i].split('_'))
            ax1.set_title(TITLE.replace('-',u'\u2013'), fontsize=18)
            ax1.title.set_y(1.01)
        # Add figure label for each subplot
        if LABELS is not None:
            at = offsetbox.AnchoredText(LABELS[i],
                loc=2, pad=0, borderpad=0.25, frameon=True,
                prop=dict(size=18,weight='bold',color='k'))
            at.patch.set_boxstyle("Square,pad=0.1")
            at.patch.set_edgecolor("white")
            ax1.axes.add_artist(at)

        # axis = equal
        ax1.set_aspect('equal', adjustable='box')
        # no ticks on the x and y axes
        ax1.get_xaxis().set_ticks([])
        ax1.get_yaxis().set_ticks([])
        # stronger linewidth on frame
        ax1.spines['geo'].set_linewidth(2.0)
        ax1.spines['geo'].set_zorder(10)
        ax1.spines['geo'].set_capstyle('projecting')

    # Add colorbar
    # Add an axes at position rect [left, bottom, width, height]
    cbar_ax = fig.add_axes([0.095, 0.105, 0.81, 0.045])
    # extend = add extension triangles to upper and lower bounds
    # options: neither, both, min, max
    cbar = fig.colorbar(im, cax=cbar_ax, extend=CBEXTEND,
        extendfrac=0.0375, drawedges=False, orientation='horizontal')
    # rasterized colorbar to remove lines
    cbar.solids.set_rasterized(True)
    # Add label to the colorbar
    cbar.ax.set_title(CBTITLE, fontsize=18, rotation=0, y=-1.65, va='top')
    cbar.ax.set_xlabel(CBUNITS, fontsize=18, rotation=0, va='center')
    cbar.ax.xaxis.set_label_coords(1.075, 0.5)
    # Set the tick levels for the colorbar
    cbar.set_ticks(levels)
    cbar.set_ticklabels([CBFORMAT.format(ct) for ct in levels])
    # ticks lines all the way across
    cbar.ax.tick_params(which='both', width=1, length=22, labelsize=18,
        direction='in')

    # adjust subplots within figure
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.16, top=0.97,
        hspace=0.05, wspace=0.05)
    # create output directory if non-existent
    FIGURE_FILE.parent.mkdir(mode=MODE, parents=True, exist_ok=True)
    # save to file
    logging.info(str(FIGURE_FILE))
    plt.savefig(FIGURE_FILE,
        metadata={'Title':pathlib.Path(sys.argv[0]).name},
        dpi=FIGURE_DPI, format=FIGURE_FORMAT)
    plt.clf()
    # change the permissions mode
    FIGURE_FILE.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description=u"""Creates 4 GMT-like plots on a global Plate Carr\u00E9e
            (Equirectangular) projection
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('infile', nargs=4,
        type=pathlib.Path,
        help='Input grid files')
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # Input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, nargs='+',
        default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input data format')
    # variable names (for ascii names of columns)
    parser.add_argument('--variables','-v',
        type=str, nargs='+', default=['lon','lat','z'],
        help='Variable names of data in input file')
    # land-sea mask
    lsmask = gravtk.utilities.get_data_path(['data','landsea_hd.nc'])
    parser.add_argument('--mask',
        type=pathlib.Path, default=lsmask,
        help='Land-sea mask')
    # output grid parameters
    parser.add_argument('--spacing',
        type=float, nargs='+', default=[0.5,0.5], metavar=('dlon','dlat'),
        help='Spatial resolution of input data')
    parser.add_argument('--interval',
        type=int, default=2, choices=[1,2],
        help=('Input grid interval (1: global, 2: centered global)'))
    # Interpolation method
    parser.add_argument('--interpolation','-I',
        type=str, default='bilinear', choices=['nearest','bilinear','cubic'],
        help='Interpolation method')
    # scale factor
    parser.add_argument('--scale-factor','-s',
        type=float, default=1.0,
        help='Multiplicative scale factor for converting to plot units')
    # plot range
    parser.add_argument('--plot-range','-R',
        type=float, nargs=3, metavar=('MIN','MAX','STEP'),
        help='Plot range and step size for normalization')
    parser.add_argument('--boundary','-B',
        type=float, nargs='+',
        help='Plot boundary for normalization')
    # color palette table or named color map
    try:
        cmap_set = set(cm.datad.keys()) | set(cm.cmaps_listed.keys())
    except (ValueError, NameError) as exc:
        cmap_set = []
    parser.add_argument('--colormap','-m',
        metavar='COLORMAP', type=str, default='viridis',
        choices=sorted(cmap_set),
        help='Named Matplotlib colormap')
    parser.add_argument('--cpt-file','-c',
        type=pathlib.Path,
        help='Input Color Palette Table (.cpt) file')
    # color map alpha
    parser.add_argument('--alpha','-a',
        type=float, default=1.0,
        help='Named Matplotlib colormap')
    # plot contour parameters
    parser.add_argument('--plot-contours',
        default=False, action='store_true',
        help='Plot contours')
    parser.add_argument('--contour-range',
        type=float, nargs=3, metavar=('MIN','MAX','STEP'),
        help='Contour range and step size')
    parser.add_argument('--mean-contour',
        default=False, action='store_true',
        help='Plot contours for mean of dataset')
    # title and label
    parser.add_argument('--plot-title', nargs=4,
        type=str, help='Plot title')
    parser.add_argument('--plot-label', nargs=4,
        type=str, help='Plot label')
    # colorbar parameters
    parser.add_argument('--cbextend',
        type=str, default='both',
        choices=['neither', 'both', 'min', 'max'],
        help='Add extension triangles to colorbar')
    parser.add_argument('--cbtitle',
        type=str, default='',
        help='Title label for colorbar')
    parser.add_argument('--cbunits',
        type=str, default='',
        help='Units label for colorbar')
    parser.add_argument('--cbformat',
        type=str, default='{0:3.0f}',
        help='Tick format for colorbar')
    # additional parameters
    parser.add_argument('--draw-grid-lines',
        default=False, action='store_true',
        help='Add map grid lines')
    parser.add_argument('--grid-lines',
        type=float, nargs='+', default=(15,15),
        help='Input grid spacing for meridians and parallels')
    # output file, format and dpi
    parser.add_argument('--figure-file','-O',
        type=pathlib.Path,
        help='Output figure file')
    parser.add_argument('--figure-format','-f',
        type=str, default='png', choices=('pdf','png','jpg','svg'),
        help='Output figure format')
    parser.add_argument('--figure-dpi','-d',
        type=int, default=180,
        help='Output figure resolution in dots per inch (dpi)')
    # print information about each input and output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of run')
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

    # create logger
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # try to run the analysis with listed parameters
    try:
        info(args)
        # run plot program with parameters
        plot_grid(args.directory, args.infile,
            DATAFORM=args.format,
            VARIABLES=args.variables,
            DDEG=args.spacing,
            INTERVAL=args.interval,
            INTERPOLATION=args.interpolation,
            SCALE_FACTOR=args.scale_factor,
            CPT_FILE=args.cpt_file,
            COLOR_MAP=args.colormap,
            PLOT_RANGE=args.plot_range,
            BOUNDARY=args.boundary,
            ALPHA=args.alpha,
            CONTOURS=args.plot_contours,
            CONTOUR_RANGE=args.contour_range,
            MEAN_CONTOUR=args.mean_contour,
            TITLES=args.plot_title,
            LABELS=args.plot_label,
            CBEXTEND=args.cbextend,
            CBTITLE=args.cbtitle,
            CBUNITS=args.cbunits,
            CBFORMAT=args.cbformat,
            DRAW_GRID_LINES=args.draw_grid_lines,
            GRID=args.grid_lines,
            FIGURE_FILE=args.figure_file,
            FIGURE_FORMAT=args.figure_format,
            FIGURE_DPI=args.figure_dpi,
            MODE=args.mode)
    except Exception as exc:
        # if there has been an error exception
        # print the type, value, and stack trace of the
        # current exception being handled
        logging.critical(f'process id {os.getpid():d} failed')
        logging.error(traceback.format_exc())

# run main program
if __name__ == '__main__':
    main()
