#!/usr/bin/env python
u"""
plot_AIS_grid_maps.py
Written by Tyler Sutterley (05/2023)
Creates GMT-like plots for the Antarctic Ice Sheet
on a polar stereographic south (3031) projection

Grounded ice and islands from Eric Rignot grounded ice image
Glacial drainage basins updated from Rignot et al. (2012) and IMBIE-2 (2016)
2003-2004 MODIS mosaic of Antarctica (MOA) from NSIDC (Haran et al., 2013)
    https://nsidc.org/data/moa/

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
    gdal: Pythonic interface to the Geospatial Data Abstraction Library (GDAL)
        https://pypi.python.org/pypi/GDAL/

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
    Updated 05/2018: using __future__ print function for python3 compatibility
    Updated 10/2017: generalizations for different grids. new BOUNDARY parameter
        added IMBIE-2 sub-basins
    Updated 09/2017: use grounded shapefiles. don't plot ice sheet with basins
    Updated 04/2017: version 2 of the Antarctic boundaries file (NSIDC)
    Updated 02/2017: direction="in" for matplotlib2.0 color bar ticks
    Updated 10/2016: added IMBIE-2 basin delineations
    Updated 06/2016: reducing MODIS mosaic when importing with GDAL
    Updated 05/2015: minor update to read and plot of drainage basins
        updates to parallel new plot_AIS_grid_movie.py code
    Written 12/2014
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
    import osgeo.gdal
except ModuleNotFoundError:
    warnings.warn("GDAL not available", ImportWarning)
try:
    import shapefile
except ModuleNotFoundError:
    warnings.warn("shapefile not available", ImportWarning)

# Antarctic 2012 basins
# region directory, filename, title and data type
region_dir = ['masks','Rignot_ANT']
region_title = ['AAp','ApB','BC','CCp','CpD','DDp','DpE','EEp','EpFp',
    'FpG','GH','HHp','HpI','IIpp','IppJ','JJpp','JppK','KKp','KpA']
# regional filenames
region_filename = 'basin_{0}_index.ascii'
# regional datatypes
region_dtype = {'names':('lat','lon'),'formats':('f','f')}

# IMBIE-2 Drainage basins
IMBIE_basin_file = ['masks','ANT_Basins_IMBIE2_v1.6','ANT_Basins_IMBIE2_v1.6.shp']
# basin titles within shapefile to extract
IMBIE_title = ('A-Ap','Ap-B','B-C','C-Cp','Cp-D','D-Dp','Dp-E','E-Ep','Ep-F',
    'F-Fp','F-G','G-H','H-Hp','Hp-I','I-Ipp','Ipp-J','J-Jpp','Jpp-K','K-A')

# background image mosaics
# MODIS mosaic of Antarctica
image_file = ['MOA','moa750_2004_hp1_v1.1.tif']

# Coastlines for antarctica (islands)
coast_file = ['masks','IceBoundaries_Antarctica_v02',
    'ant_ice_sheet_islands_v2.shp']

# Antarctica (AIS)
# x and y limit (modified from Bamber 1km DEM)
xlimits = np.array([-3100000,3100000])
ylimits = np.array([-2600000,2600000])
# cartopy transform for polar stereographic south
try:
    projection = ccrs.Stereographic(central_longitude=0.0,
        central_latitude=-90.0,true_scale_latitude=-71.0)
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

# PURPOSE: plot Rignot 2012 drainage basin polylines
def plot_rignot_basins(ax, base_dir):
    region_directory = base_dir.joinpath(*region_dir)
    # for each region
    for reg in region_title:
        # read the regional polylines
        region_file = region_directory.joinpath(region_filename.format(reg))
        region_ll = np.loadtxt(region_file, dtype=region_dtype)
        # converting region lat/lon into plot coordinates
        points = projection.transform_points(ccrs.PlateCarree(),
            region_ll['lon'], region_ll['lat'])
        ax.plot(points[:,0], points[:,1], color='k', transform=projection)

# PURPOSE: plot Antarctic drainage basins from IMBIE2 (Mouginot)
def plot_IMBIE2_basins(ax, base_dir):
    # read drainage basin polylines from shapefile (using splat operator)
    basin_shapefile = base_dir.joinpath(*IMBIE_basin_file)
    logging.debug(str(basin_shapefile))
    logging.debug(str(basin_shapefile))
    shape_input = shapefile.Reader(str(basin_shapefile))
    shape_entities = shape_input.shapes()
    shape_attributes = shape_input.records()
    # find record index for region by iterating through shape attributes
    # no islands or large regions
    i=[i for i,a in enumerate(shape_attributes) if a[1] in IMBIE_title]
    # for each valid shape entity
    for indice in i:
        # extract Polar-Stereographic coordinates for record
        points = np.array(shape_entities[indice].points)
        # IMBIE-2 basins can have multiple parts
        parts = shape_entities[indice].parts
        parts.append(len(points))
        for p1,p2 in zip(parts[:-1],parts[1:]):
            ax.plot(points[p1:p2,0], points[p1:p2,1], c='k',
                transform=projection)

# PURPOSE: plot Antarctic drainage sub-basins from IMBIE-2 (Mouginot)
def plot_IMBIE2_subbasins(ax, base_dir):
    # read drainage basin polylines from shapefile (using splat operator)
    IMBIE_subbasin_file = ['Basins_20Oct2016_v1.7','Basins_v1.7.shp']
    basin_shapefile = base_dir.joinpath('masks',*IMBIE_subbasin_file)
    logging.debug(str(basin_shapefile))
    logging.debug(str(basin_shapefile))
    shape_input = shapefile.Reader(str(basin_shapefile))
    shape_entities = shape_input.shapes()
    shape_attributes = shape_input.records()
    # iterate through shape entities and attributes
    indices = [i for i,a in enumerate(shape_attributes) if (a[1] != 'Islands')]
    for i in indices:
        # extract Polar-Stereographic coordinates for record
        points = np.array(shape_entities[i].points)
        # IMBIE-2 basins can have multiple parts
        parts = shape_entities[i].parts
        parts.append(len(points))
        for p1,p2 in zip(parts[:-1],parts[1:]):
            ax.plot(points[p1:p2,0], points[p1:p2,1], c='k',
                transform=projection)

# PURPOSE: plot Antarctic grounded ice delineation
def plot_grounded_ice(ax, base_dir, START=1):
    grounded_ice_shapefile = base_dir.joinpath(*coast_file)
    logging.debug(str(grounded_ice_shapefile))
    shape_input = shapefile.Reader(str(grounded_ice_shapefile))
    shape_entities = shape_input.shapes()
    shape_attributes = shape_input.records()
    i = [i for i,e in enumerate(shape_entities) if (np.ndim(e.points) > 1)]
    for indice in i[START:]:
        # extract Polar-Stereographic coordinates for record
        points = np.array(shape_entities[indice].points)
        ax.plot(points[:,0], points[:,1], c='k',
            transform=projection)

# PURPOSE: plot MODIS mosaic of Antarctica as background image
def plot_image_mosaic(ax, base_dir, MASKED=True):
    # read MODIS mosaic of Antarctica
    image_geotiff_file = base_dir.joinpath(*image_file)
    logging.debug(str(image_geotiff_file))
    ds = osgeo.gdal.Open(str(image_geotiff_file))
    # get dimensions
    xsize = ds.RasterXSize
    ysize = ds.RasterYSize
    # get geotiff info
    info_geotiff = ds.GetGeoTransform()
    # calculate image extents
    xmin = info_geotiff[0]
    ymax = info_geotiff[3]
    xmax = xmin + (xsize-1)*info_geotiff[1]
    ymin = ymax + (ysize-1)*info_geotiff[5]
    # read as grayscale image
    mosaic = np.ma.array(ds.ReadAsArray())
    # mask image mosaic
    if MASKED:
        # mask invalid values
        mosaic.fill_value = 0
        # create mask array for bad values
        mosaic.mask = (mosaic.data == mosaic.fill_value)
    # image extents
    extents=(xmin,xmax,ymin,ymax)
    # dataset range
    vmin, vmax = (0, 16386)
    # create color map with transparent bad points
    image_cmap = copy.copy(cm.gist_gray)
    image_cmap.set_bad(alpha=0.0)
    # nearest to not interpolate image
    im = ax.imshow(mosaic, interpolation='nearest', extent=extents,
        cmap=image_cmap, vmin=vmin, vmax=vmax, origin='upper',
        transform=projection)
    im.set_rasterized(True)
    # close the dataset
    ds = None

# PURPOSE: add a plot scale
def add_plot_scale(ax,X,Y,dx,dy,masked,fc1='w',fc2='k'):
    if masked:
        x1,x2,y1,y2 = [X-0.1*dx,X+1.2*dx,Y-2.5*dy,Y+3.2*dy]
        ax.fill([x1,x2,x2,x1,x1], [y1,y1,y2,y2,y1], fc1, zorder=4)
    for i,c in enumerate([fc1,fc2,fc1,fc2]):
        x1,x2,y1,y2 = [X+0.25*i*dx,X+0.25*(i+1)*dx,Y,Y+dy]
        ax.fill([x1,x2,x2,x1,x1], [y1,y1,y2,y2,y1], c, zorder=5)
    ax.plot([X,X+dx,X+dx,X,X], [Y,Y,Y+dy,Y+dy,Y], fc2, zorder=6)
    for i in range(3):
        ax.plot([X+0.5*i*dx,X+0.5*i*dx], [Y,Y-0.5*dy], fc2, zorder=6)
        ax.text(X+0.5*i*dx, Y-0.9*dy, '{0:0.0f}'.format(0.5*i*dx/1e3),
            ha='center', va='top', fontsize=12, color=fc2, zorder=6)
    ax.text(X+0.5*dx, Y+1.3*dy, 'km', ha='center', va='bottom',
        fontsize=12, color=fc2, zorder=6)

# plot grid program
def plot_grid(base_dir, FILENAME,
    DATAFORM=None,
    VARIABLES=[],
    MASK=None,
    INTERPOLATION=None,
    DDEG=None,
    INTERVAL=None,
    SCALE_FACTOR=1.0,
    COLOR_MAP=None,
    CPT_FILE=None,
    PLOT_RANGE=None,
    BOUNDARY=None,
    ALPHA=1.0,
    CONTOURS=False,
    CONTOUR_RANGE=None,
    MEAN_CONTOUR=False,
    TITLE=None,
    LABEL=None,
    CBEXTEND=None,
    CBTITLE=None,
    CBUNITS=None,
    CBFORMAT=None,
    BASEMAP=False,
    BASIN_TYPE=None,
    DRAW_GRID_LINES=False,
    GRID=None,
    DRAW_SCALE=False,
    FIGURE_FILE=None,
    FIGURE_FORMAT=None,
    FIGURE_DPI=None,
    MODE=0o775):

    # read CPT or use color map
    if CPT_FILE is not None:
        # cpt file
        cmap = gravtk.tools.from_cpt(CPT_FILE)
    else:
        # colormap
        cmap = copy.copy(cm.get_cmap(COLOR_MAP))

    # if using MODIS mosaic of Antarctica as basemap
    if BASEMAP:
        # transparent color map for bad values
        cmap.set_bad(alpha=0.0)
    else:
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

    # input ascii/netCDF4/HDF5 file
    if (DATAFORM == 'ascii'):
        # ascii (.txt)
        dinput = gravtk.spatial().from_ascii(FILENAME, date=False,
            columns=VARIABLES, spacing=[dlon,dlat], nlat=nlat, nlon=nlon)
    elif (DATAFORM == 'netCDF4'):
        # netCDF4 (.nc)
        field_mapping = gravtk.spatial().default_field_mapping(VARIABLES)
        dinput = gravtk.spatial().from_netCDF4(FILENAME, date=False,
            field_mapping=field_mapping)
    elif (DATAFORM == 'HDF5'):
        # HDF5 (.H5)
        field_mapping = gravtk.spatial().default_field_mapping(VARIABLES)
        dinput = gravtk.spatial().from_HDF5(FILENAME, date=False,
            field_mapping=field_mapping)

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
        # update mask
        dinput.replace_invalid(fill_value=dinput.fill_value, mask=mask)

    # scale input dataset
    if (SCALE_FACTOR != 1.0):
        dinput = dinput.scale(SCALE_FACTOR)

    # if dlat is negative
    if (np.sign(dlat) == -1):
        dinput = dinput.flip(axis=0)

    # setup stereographic map
    fig, ax1 = plt.subplots(num=1, nrows=1, ncols=1, figsize=(10,7.5),
        subplot_kw=dict(projection=projection))
    # WGS84 Ellipsoid parameters
    a_axis = 6378137.0# [m] semimajor axis of the ellipsoid
    flat = 1.0/298.257223563# flattening of the ellipsoid
    # (4pi/3)R^3 = (4pi/3)(a^2)b = (4pi/3)(a^3)(1 -f)
    rad_e = a_axis*(1.0 - flat)**(1.0/3.0)

    # plot image of MODIS mosaic of Antarctica as base layer
    if BASEMAP:
        # plot MODIS mosaic of Antarctica
        plot_image_mosaic(ax1, base_dir)

    # calculate image coordinates
    mx = np.int64((xlimits[1]-xlimits[0])/1000.)+1
    my = np.int64((ylimits[1]-ylimits[0])/1000.)+1
    X = np.linspace(xlimits[0],xlimits[1],mx)
    Y = np.linspace(ylimits[0],ylimits[1],my)
    gridx,gridy = np.meshgrid(X,Y)
    # create mesh lon/lat
    points = ccrs.PlateCarree().transform_points(projection,
        gridx.flatten(), gridy.flatten())
    lonsin = points[:,0].reshape(my,mx)
    latsin = points[:,1].reshape(my,mx)

    # interpolate to image coordinates
    if (INTERVAL == 1) and (np.max(dinput.lon) > 180):# (0:360, 90:-90)
        shift_data,lon180 = gravtk.tools.shift_grid(180.0,dinput.data,dinput.lon)
        shift_mask,lon180 = gravtk.tools.shift_grid(180.0,dinput.mask,dinput.lon)
        img = gravtk.tools.interp_grid(shift_data,lon180,dinput.lat,lonsin,latsin,order)
        msk = gravtk.tools.interp_grid(shift_mask,lon180,dinput.lat,lonsin,latsin,order)
    elif (INTERVAL == 2) and (np.max(dinput.lon) > 180):# DDEG/2
        shift_data,lon180 = gravtk.tools.shift_grid(180.0+dlon,dinput.data,dinput.lon)
        shift_mask,lon180 = gravtk.tools.shift_grid(180.0+dlon,dinput.mask,dinput.lon)
        img = gravtk.tools.interp_grid(shift_data,lon180,dinput.lat,lonsin,latsin,order)
        msk = gravtk.tools.interp_grid(shift_mask,lon180,dinput.lat,lonsin,latsin,order)
    else:# -180:180 or modification of there of
        img = gravtk.tools.interp_grid(dinput.data,dinput.lon,dinput.lat,lonsin,latsin,order)
        msk = gravtk.tools.interp_grid(dinput.mask,dinput.lon,dinput.lat,lonsin,latsin,order)
    # create masked array of image
    img = np.ma.array(img, mask=msk.astype(bool))

    # plot only grounded points
    if MASK is not None:
        img = gravtk.tools.mask_oceans(lonsin, latsin,
            data=img, order=order, iceshelves=False)
    # plot image with transparency using normalization
    im = ax1.imshow(img, interpolation='nearest', cmap=cmap,
        extent=(xlimits[0],xlimits[1],ylimits[0],ylimits[1]),
        norm=norm, alpha=ALPHA, origin='lower', transform=projection)

    # create mesh lon/lat
    lon, lat = np.meshgrid(dinput.lon,dinput.lat)
    data = dinput.to_masked_array()
    # plot line contours
    if CONTOURS:
        clevs = np.arange(CONTOUR_RANGE[0],
            CONTOUR_RANGE[1] + CONTOUR_RANGE[2],
            CONTOUR_RANGE[2])
        # remove 0 (will plot in red)
        reduce_clevs = clevs[np.nonzero(clevs)]
        # plot contours
        ax1.contour(lon,lat,data,reduce_clevs,colors='0.2',linestyles='solid',
            transform=ccrs.PlateCarree())
        ax1.contour(lon,lat,data,[0],colors='red',linestyles='solid',
            linewidths=1.5,transform=ccrs.PlateCarree())

    # plot line contour for global average
    if MEAN_CONTOUR and CONTOURS:
        # calculate areas of each grid cell
        dphi,dth = (dlon*np.pi/180.0,dlat*np.pi/180.0)
        indy,indx = np.nonzero(np.logical_not(data.mask))
        area = (rad_e**2)*dth*dphi*np.cos(lat[indy,indx]*np.pi/180.0)
        # calculate average
        ave = np.sum(area*data[indy,indx])/np.sum(area)
        # plot line contour of global average
        ax1.contour(lon,lat,data,[ave],colors='blue',linestyles='solid',
            linewidths=1.5,transform=ccrs.PlateCarree())

    # add basins based on BASIN_TYPE (Rignot 2012, IMBIE-2, IMBIE-2 subbasins)
    if (BASIN_TYPE == 'Rignot'):
        plot_rignot_basins(ax1, base_dir)
        start_indice = 1
    elif (BASIN_TYPE == 'IMBIE-2'):
        plot_IMBIE2_basins(ax1, base_dir)
        start_indice = 1
    elif (BASIN_TYPE == 'IMBIE-2_subbasin'):
        plot_IMBIE2_subbasins(ax1, base_dir)
        start_indice = 1
    else:
        start_indice = 0
    # add glaciers/ice caps/islands
    plot_grounded_ice(ax1, base_dir, START=start_indice)

    # draw lat/lon grid lines
    if DRAW_GRID_LINES:
        # meridian and parallel grid spacing
        llx,lly = (GRID[0],GRID[0]) if (len(GRID) == 1) else (GRID[0],GRID[1])
        grid_meridians = np.arange(-180, 180 + llx, llx)
        grid_parallels = np.arange(-90, 90 + lly, lly)
        gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
            linewidth=0.1, color='0.25', linestyle='-')
        gl.xlocator = ticker.FixedLocator(grid_meridians)
        gl.ylocator = ticker.FixedLocator(grid_parallels)

    # Add colorbar and adjust size
    # pad = distance from main plot axis
    # extend = add extension triangles to upper and lower bounds
    # options: neither, both, min, max
    # shrink = percent size of colorbar
    # aspect = lengthXwidth aspect of colorbar
    cbar = plt.colorbar(im, ax=ax1, pad=0.025, extend=CBEXTEND,
        extendfrac=0.0375, shrink=0.925, aspect=20, drawedges=False)
    # rasterized colorbar to remove lines
    cbar.solids.set_rasterized(True)
    # Add label to the colorbar
    cbar.ax.set_ylabel(CBTITLE, labelpad=10, fontsize=24)
    cbar.ax.set_title(CBUNITS, fontsize=24, va='bottom')
    # Set the tick levels for the colorbar
    cbar.set_ticks(levels)
    cbar.set_ticklabels([CBFORMAT.format(ct) for ct in levels])
    # ticks lines all the way across
    cbar.ax.tick_params(which='both', width=1, length=23, labelsize=24,
        direction='in')

    # x and y limits, axis = equal
    ax1.set_xlim(xlimits)
    ax1.set_ylim(ylimits)
    ax1.set_aspect('equal', adjustable='box')
    # no ticks on the x and y axes
    ax1.get_xaxis().set_ticks([])
    ax1.get_yaxis().set_ticks([])

    # add main title
    if TITLE is not None:
        ax1.set_title(TITLE.replace('-',u'\u2013'), fontsize=24)
    # Add figure label
    if LABEL is not None:
        if BASEMAP:
            at = offsetbox.AnchoredText(LABEL,
                loc=2, pad=0, frameon=True,
                prop=dict(size=24,weight='bold'))
            at.patch.set_boxstyle("Square,pad=0.25")
            at.patch.set_edgecolor("white")
        else:
            at = offsetbox.AnchoredText(LABEL,
                loc=2, pad=0, frameon=False,
                prop=dict(size=24,weight='bold'))
        ax1.axes.add_artist(at)

    # draw map scale to corners
    if DRAW_SCALE:
        add_plot_scale(ax1,-295e4,-236e4,1000e3,85e3,False)
    # stronger linewidth on frame
    ax1.spines['geo'].set_linewidth(2.0)
    ax1.spines['geo'].set_zorder(10)
    ax1.spines['geo'].set_capstyle('projecting')

    # adjust subplot within figure
    fig.subplots_adjust(left=0.02,right=0.98,bottom=0.01,top=0.97)
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
        description="""Creates GMT-like plots of the Antarctic ice sheet
            on a polar stereographic south (EPSG 3031) projection
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('infile',
        type=pathlib.Path,
        help='Input grid file')
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # Input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
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
    parser.add_argument('--plot-title',
        type=str, help='Plot title')
    parser.add_argument('--plot-label',
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
    parser.add_argument('--basemap',
        default=False, action='store_true',
        help='Add background basemap image')
    parser.add_argument('--basin-type',
        type=str, default='',
        help='Add delineations for glacier drainage basins')
    parser.add_argument('--draw-grid-lines',
        default=False, action='store_true',
        help='Add map grid lines')
    parser.add_argument('--grid-lines',
        type=float, nargs='+', default=(15,15),
        help='Input grid spacing for meridians and parallels')
    parser.add_argument('--draw-scale',
        default=False, action='store_true',
        help='Add map scale bar')
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
            TITLE=args.plot_title,
            LABEL=args.plot_label,
            CBEXTEND=args.cbextend,
            CBTITLE=args.cbtitle,
            CBUNITS=args.cbunits,
            CBFORMAT=args.cbformat,
            BASEMAP=args.basemap,
            BASIN_TYPE=args.basin_type,
            DRAW_GRID_LINES=args.draw_grid_lines,
            GRID=args.grid_lines,
            DRAW_SCALE=args.draw_scale,
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
