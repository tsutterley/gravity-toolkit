from gravity_toolkit.gauss_weights import gauss_weights
from gravity_toolkit.gen_stokes import gen_stokes
from gravity_toolkit.harmonics import harmonics
from gravity_toolkit.harmonic_summation import harmonic_summation
from gravity_toolkit.plm_holmes import plm_holmes
from gravity_toolkit.read_love_numbers import read_love_numbers
from gravity_toolkit.spatial import spatial
from gravity_toolkit.units import units
from gravity_toolkit.utilities import get_data_path

import numpy as np
import scipy.signal as sg
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.animation as animation
import cartopy.crs as ccrs
from IPython.display import HTML


def create_grid(Ylms, lmax=None, rad=0, destripe=False, unit='cmwe', dlon=0.5, dlat=0.5, bounds=None):
    """
    Function to convert a harmonic object to grid format

    Parameters
    ----------
    Ylms : harmonics object to convert to grid format
    lmax : maximum degree of spherical harmonics used
    rad : radius of the gaussian filter. If set to 0, no gaussian filter is apply
    destripe : boolean to apply or not the destripe method of harmonics
    unit : unit of the grid in ['cmwe', 'cmweEl', 'geoid', 'cmwe_ne', 'microGal']
    dlon : output longitude spacing
    dlat : output latitude spacing
    bounds : list with [lon_max, lon_min, lat_max, lat_min]

    Returns
    -------
    grid : spatial object with the grid converted from the original harmonics object
    """
    # Output spatial data
    grid = spatial()
    grid.time = np.copy(Ylms.time)
    grid.month = np.copy(Ylms.month)

    # Output Degree Interval
    if bounds is None:
        grid.lon = np.arange(-180 + dlon / 2.0, 180 + dlon / 2.0, dlon)
        grid.lat = np.arange(90.0 - dlat / 2.0, -90.0 - dlat / 2.0, -dlat)
    else:
        grid.lon = np.arange(-bounds[1] + dlon / 2.0, bounds[0] + dlon / 2.0, dlon)
        grid.lat = np.arange(bounds[2] - dlat / 2.0, -bounds[3] - dlat / 2.0, -dlat)

    nlon = len(grid.lon)
    nlat = len(grid.lat)

    # update spacing and dimensions
    grid.update_spacing()
    grid.update_extents()
    grid.update_dimensions()

    # Computing plms for converting to spatial domain
    theta = (90.0 - grid.lat) * np.pi / 180.0
    if lmax is None:
        PLM, dPLM = plm_holmes(Ylms.lmax, np.cos(theta))
    else:
        PLM, dPLM = plm_holmes(lmax, np.cos(theta))

    # read load love numbers file
    love_numbers_file = get_data_path(['data', 'love_numbers'])
    # LMAX of load love numbers from Han and Wahr (1995) is 696.
    # from Wahr (2007) linearly interpolating kl worksand ll Love Numbers
    hl, kl, ll = read_love_numbers(love_numbers_file, REFERENCE='CF')

    if unit == 'cmwe':
        dfactor = units(lmax=Ylms.lmax).harmonic(hl, kl, ll).cmwe
    elif unit == 'cmweEl':
        dfactor = units(lmax=Ylms.lmax).harmonic(hl, kl, ll).cmweEL
    elif unit == 'geoid':
        dfactor = units(lmax=Ylms.lmax).harmonic(hl, kl, ll).mmGH
    elif unit == 'cmwe_ne':
        dfactor = units(lmax=Ylms.lmax).harmonic(hl, kl, ll).cmwe_ne
    elif unit == 'microGal':
        dfactor = units(lmax=Ylms.lmax).harmonic(hl, kl, ll).microGal
    else:
        raise ValueError("Unit not accepted, should be either 'cmwe' pr 'cmweEl' or 'cmwe_ne' or 'geoid' or 'microGal'")

    # converting harmonics to truncated, smoothed coefficients in units
    # combining harmonics to calculate output spatial fields
    # output spatial grid
    if not (type(Ylms.month) in [list, np.array]) and len(Ylms.month) == 1:
        grid.data = np.zeros((nlat, nlon))

        if destripe:
            tmp = Ylms.destripe()
        else:
            tmp = Ylms

        if rad != 0:
            wt = 2.0 * np.pi * gauss_weights(rad, lmax)
            tmp.convolve(dfactor * wt)
        else:
            tmp.convolve(dfactor)
        # convert spherical harmonics to output spatial grid
        if lmax is None:
            grid.data[:, :] = harmonic_summation(tmp.clm, tmp.slm,
                                                 grid.lon, grid.lat, LMAX=Ylms.lmax, MMAX=Ylms.mmax, PLM=PLM).T
        else:
            grid.data[:, :] = harmonic_summation(tmp.clm, tmp.slm,
                                                 grid.lon, grid.lat, LMAX=lmax, MMAX=lmax, PLM=PLM).T

    else:
        grid.data = np.zeros((nlat, nlon, len(Ylms.month)))
        for i, grace_month in enumerate(Ylms.month):
            # GRACE/GRACE-FO harmonics for time t
            # convert to output units
            if destripe:
                tmp = Ylms.index(i).destripe()
            else:
                tmp = Ylms.index(i)

            if rad != 0:
                wt = 2.0 * np.pi * gauss_weights(rad, lmax)
                tmp.convolve(dfactor * wt)
            else:
                tmp.convolve(dfactor * np.ones((Ylms.lmax + 1)))
            # convert spherical harmonics to output spatial grid
            if lmax is None:
                grid.data[:, :, i] = harmonic_summation(tmp.clm, tmp.slm,
                                                        grid.lon, grid.lat, LMAX=Ylms.lmax, MMAX=Ylms.mmax, PLM=PLM).T
            else:
                grid.data[:, :, i] = harmonic_summation(tmp.clm, tmp.slm,
                                                        grid.lon, grid.lat, LMAX=lmax, MMAX=lmax, PLM=PLM).T

    grid.mask = np.zeros(grid.data.shape)
    return grid


def grid_to_hs(grid, lmax, mmax=None, unit='cmwe'):
    """
    Function to convert spatial object (grid) to harmonics object (spherical harmonics)

    Parameters
    ----------
    grid : spatial object to convert to harmonics
    lmax : maximal degree of the harmonics object to create
    mmax : maximal order of the harmonics object to create
    unit : unit of the grid in ['cmwe', 'geoid', 'cmwe_ne', 'microGal']

    Returns
    -------
    harmonics : harmonics object
    """
    # -- load love numbers
    hl, kl, ll = read_love_numbers(get_data_path(['data', 'love_numbers']), REFERENCE='CF')

    # -- set maximum spherical harmonic order
    mmax = np.copy(lmax) if (mmax is None) else mmax

    # -- number of dates in data
    if type(grid.time) in [list, np.array] or len(grid.time) != 1:
        n_time = len(grid.time)
    else:
        n_time = 1
    # -- Spherical harmonic coefficient matrices to be filled from data file
    grace_clm = np.zeros((lmax + 1, mmax + 1, n_time))
    grace_slm = np.zeros((lmax + 1, mmax + 1, n_time))
    # -- output dimensions
    lout = np.arange(lmax + 1)
    mout = np.arange(mmax + 1)

    # -- Test to attribute UNITS number
    if unit == 'cmwe':
        UNITS = 1
    elif unit == 'geoid':
        UNITS = 4
    elif unit == 'cmwe_ne':
        UNITS = 6
    elif unit == 'microGal':
        UNITS = 5
    else:
        raise ValueError("Unit not accepted, should be either 'cmwe' or 'cmwe_ne' or 'geoid' or 'microGal'")

    # -- for each date, conversion to spherical harmonics
    if n_time != 1:
        for i in range(n_time):
            harmo = gen_stokes(grid.data[:, :, i],
                               grid.lon[:], grid.lat[:],
                               LMAX=lmax, MMAX=mmax, UNITS=UNITS, LOVE=(hl, kl, ll))

            grace_clm[:, :, i] = harmo.clm
            grace_slm[:, :, i] = harmo.slm

    else:
        print('mono grid_to_hs')
        harmo = gen_stokes(grid.data[:, :],
                           grid.lon[:], grid.lat[:],
                           LMAX=lmax, MMAX=mmax, UNITS=UNITS, LOVE=(hl, kl, ll))

        grace_clm[:, :, 0] = harmo.clm
        grace_slm[:, :, 0] = harmo.slm

    # -- return the GRACE data, GRACE date (mid-month in decimal), and the
    # -- start and end days as Julian dates
    result_dict = {'clm': grace_clm, 'slm': grace_slm, 'time': grid.time, 'month': grid.month,
            'l': lout, 'm': mout, 'title': '', 'directory': ''}

    return harmonics().from_dict(result_dict)


def diff_grid(grid1, grid2):
    """
    Create a grid resulting from the difference between the two given grids

    Parameters
    ----------
    grid1 : spatial object
    grid2 : spatial object to substract to the first

    Returns
    -------
    grid : spatial object with the difference between both grid
    """
    exclude1 = set(grid1.month) - set(grid2.month)

    # Output spatial data
    grid = spatial()
    grid.month = np.array(list(sorted(set(grid1.month) - exclude1)))
    grid.time = np.array([grid1.time[i] for i in range(len(grid1.time)) if not (grid1.month[i] in exclude1)])

    # Output Degree Interval
    grid.lon = grid1.lon
    grid.lat = grid1.lat

    # update spacing and dimensions
    grid.update_spacing()
    grid.update_extents()
    grid.update_dimensions()

    grid.data = np.zeros((grid.lat.shape[0], grid.lon.shape[0], len(grid.month)))
    cmp = 0
    for i in range(len(grid1.month)):
        for j in range(len(grid2.month)):
            if grid1.month[i] == grid2.month[j]:
                grid.data[:, :, cmp] = grid1.data[:, :, i] - grid2.data[:, :, j]
                cmp += 1

    return grid


def filt_Ylms(ylms, filt='low', filt_param=None):
    """
    Apply a temporal filter on harmonics object

    Parameters
    ----------
    ylms : harmonics object to filter
    filt : choice of the filter in ['low', 'band', 'fft']
    filt_param : cut frequency of the filter. For band filter, a list with (f_max, f_min)

    Returns
    -------
    filtered_ylms : temporally filtered harmonics object

    """
    filtered_ylms = ylms.copy()

    # len of the data
    ndata = filtered_ylms.time.shape[0]
    # compute the mean time delta of the object
    dt = float(np.mean((filtered_ylms.time[1:] - filtered_ylms.time[:-1])))

    if filt_param is not None and type(filt_param) != list:
        filt_param = [filt_param]

    if filt == 'low':
        if filt_param is None:
            b, a = sg.butter(10, 0.5, analog=False, fs=1 / dt)
        else:
            b, a = sg.butter(10, filt_param[0], analog=False, fs=1 / dt)

        for i in range(filtered_ylms.clm.shape[0]):
            for j in range(filtered_ylms.clm.shape[1]):
                filtered_ylms.clm[i, j] = sg.filtfilt(b, a, filtered_ylms.clm[i, j])
                filtered_ylms.slm[i, j] = sg.filtfilt(b, a, filtered_ylms.slm[i, j])

    elif filt == 'band':
        if filt_param is None:
            b, a = sg.butter(6, 0.3, analog=False, fs=1 / dt)
            b2, a2 = sg.butter(6, 0.04, btype='highpass', analog=False, fs=1 / dt)
        else:
            b, a = sg.butter(6, filt_param[0], analog=False, fs=1 / dt)
            b2, a2 = sg.butter(6, filt_param[1], btype='highpass', analog=False, fs=1 / dt)

        for i in range(filtered_ylms.clm.shape[0]):
                for j in range(filtered_ylms.clm.shape[1]):
                    filtered_ylms.clm[i, j] = sg.filtfilt(b, a, filtered_ylms.clm[i, j])
                    filtered_ylms.clm[i, j] = sg.filtfilt(b2, a2, filtered_ylms.clm[i, j])
                    filtered_ylms.slm[i, j] = sg.filtfilt(b, a, filtered_ylms.slm[i, j])
                    filtered_ylms.slm[i, j] = sg.filtfilt(b2, a2, filtered_ylms.slm[i, j])

    elif filt == 'fft':
        # zero pad
        n2 = 0
        while ndata > 2 ** n2:
            n2 += 1
        n2 += 1

        fc = np.fft.fft(filtered_ylms.clm, n=2 ** n2, axis=2)
        fs = np.fft.fft(filtered_ylms.slm, n=2 ** n2, axis=2)
        freq = np.fft.fftfreq(2 ** n2, d=dt)
        if filt_param is None:
            to_zero = np.logical_or(freq > 0.5, freq < -0.5)
        else:
            to_zero = np.logical_or(freq > filt_param[0], freq < filt_param[0])
        fc[:, :, to_zero] = 0
        fs[:, :, to_zero] = 0
        filtered_ylms.clm = np.real(np.fft.ifft(fc, axis=2))[:, :, :ndata]
        filtered_ylms.slm = np.real(np.fft.ifft(fs, axis=2))[:, :, :ndata]

    return filtered_ylms


def filt_grid(grid, f_cut=0.5):
    """
    Temporally filter a grid with a truncation in fft at 2 years

    Parameters
    ----------
    grid : spatial object to filter
    f_cut : cutting frequency

    Returns
    -------
    filtered_grid : spatial object filtered

    """
    filtered_grid = grid.copy()
    time = grid.time
    ndata = grid.time.shape[0]

    # zero pad
    n2 = 0
    while ndata > 2 ** n2:
        n2 += 1
    n2 += 1

    # compute the mean time delta of the object
    dt = float(np.mean((time[1:] - time[:-1])))

    f = np.fft.fft(grid.data, n=2 ** n2, axis=2)
    freq = np.fft.fftfreq(2 ** n2, d=dt)

    to_zero = np.logical_or(freq > f_cut, freq < -f_cut)
    f[:, :, to_zero] = 0
    filtered_grid.data = np.real(np.fft.ifft(f, axis=2))[:, :, :ndata]

    return filtered_grid


def save_gif(grid, path, unit='cmwe', bound=None, mask=None, color='viridis'):
    """
    Create a gif of the spatial object

    Parameters
    ----------
    grid : spatial object to convert to gif
    path : path of the future gif (mandatory to end in .gif)
    unit : unit of the grid in ['cmwe', 'mmwe', 'geoid', 'cmwe_ne', 'microGal', 'secacc']
    bound : list with minimal value and maximal value of the colorbar. Default value is None
    mask : np.array corresponding to the mask
    color : matplotlib cmap color of the gif (Recommended: viridis, plasma, RdBu_r)
    """
    matplotlib.rcParams['animation.embed_limit'] = 2**128

    if mask is None:
        data_to_set = grid.data
    else:
        data_to_set = grid.data*mask

    fig, ax1 = plt.subplots(num=1, nrows=1, ncols=1, figsize=(10.375,6.625),
        subplot_kw=dict(projection=ccrs.PlateCarree()))

    # levels and normalization for plot range
    print(np.min(data_to_set), np.max(data_to_set))
    if bound is None:
        vmin, vmax = int(np.min(data_to_set)), int(np.ceil(np.max(data_to_set)))
    else:
        vmin, vmax = bound

    if vmax - vmin >= 3:
        levels = np.arange(vmin, vmax, max(1, int((vmax - vmin)/10)))
    elif vmax - vmin >= 0.3:
        levels = np.arange(vmin, vmax, max(0.1, float('%.1f'%((vmax - vmin)/10))))
    elif vmax - vmin >= 0.03:
        levels = np.arange(vmin, vmax, max(0.1, float('%.2f'%((vmax - vmin)/10))))
    else:
        raise ValueError("The range of data to plot is too small")

    norm = colors.Normalize(vmin=vmin,vmax=vmax)
    cmap = plt.cm.get_cmap(color)
    im = ax1.imshow(np.zeros((np.int(180.0 + 1.0),np.int(360.0 + 1.0))), interpolation='nearest',
        norm=norm, cmap=cmap, transform=ccrs.PlateCarree(),
        extent=grid.extent, origin='upper', animated=True)
    ax1.coastlines('50m')

    # add date label
    time_text = ax1.text(0.025, 0.025, '', transform=fig.transFigure,
        color='k', size=24, ha='left', va='baseline')

    # Add horizontal colorbar and adjust size
    # extend = add extension triangles to upper and lower bounds
    # options: neither, both, min, max
    # pad = distance from main plot axis
    # shrink = percent size of colorbar
    # aspect = lengthXwidth aspect of colorbar
    cbar = plt.colorbar(im, ax=ax1, extend='both', extendfrac=0.0375,
        orientation='horizontal', pad=0.025, shrink=0.85,
        aspect=22, drawedges=False)
    # rasterized colorbar to remove lines
    cbar.solids.set_rasterized(True)
    # Add label to the colorbar
    if unit == "cmwe":
        cbar.ax.set_xlabel('Equivalent Water Thickness', labelpad=10, fontsize=24)
        cbar.ax.set_ylabel('cm', fontsize=24, rotation=0)
    elif unit == "cmwe_ne":
        cbar.ax.set_xlabel('Non elastic Equivalent Water Thickness', labelpad=10, fontsize=24)
        cbar.ax.set_ylabel('cm', fontsize=24, rotation=0)
    elif unit == "mmwe":
        cbar.ax.set_xlabel('Equivalent Water Thickness', labelpad=10, fontsize=24)
        cbar.ax.set_ylabel('mm', fontsize=24, rotation=0)
    elif unit == "geoid":
        cbar.ax.set_xlabel('Geoid Height', labelpad=10, fontsize=24)
        cbar.ax.set_ylabel('mm', fontsize=24, rotation=0)
    elif unit == "microGal":
        cbar.ax.set_xlabel('Acceleration', labelpad=10, fontsize=24)
        cbar.ax.set_ylabel('$\mu Gal$', fontsize=24, rotation=0)
    elif unit == "secacc":
        cbar.ax.set_xlabel('Secular Acceleration', labelpad=10, fontsize=24)
        cbar.ax.set_ylabel('$nT.y^{-2}$', fontsize=24, rotation=0)

    cbar.ax.yaxis.set_label_coords(1.045, 0.1)
    # Set the tick levels for the colorbar
    cbar.set_ticks(levels)
    if vmax - vmin >= 3:
        cbar.set_ticklabels(['{0:d}'.format(ct) for ct in levels])
    elif vmax - vmin >= 0.3:
        cbar.set_ticklabels(['{.1f}'.format(ct) for ct in levels])
    elif vmax - vmin >= 0.03:
        cbar.set_ticklabels(['{.2f}'.format(ct) for ct in levels])
    # ticks lines all the way across
    cbar.ax.tick_params(which='both', width=1, length=26, labelsize=24,
        direction='in')

    # stronger linewidth on frame
    ax1.spines['geo'].set_linewidth(2.0)
    ax1.spines['geo'].set_capstyle('projecting')
    # adjust subplot within figure
    fig.subplots_adjust(left=0.02,right=0.98,bottom=0.05,top=0.98)

    # animate frames
    def animate_frames(i):
        # set image
        im.set_data(data_to_set[:,:,i])
        # add date label
        time_text.set_text('{:.2f}'.format(grid.time[i]))

    # set animation
    anim = animation.FuncAnimation(fig, animate_frames, frames=len(grid.month))
    HTML(anim.to_jshtml())

    anim.save(path, writer='imagemagick', fps=10)
    plt.clf()


def plot_rms_map(grid, path=False, proj=ccrs.PlateCarree(), unit='cmwe', bound=None, mask=None, color='viridis'):
    """
    Create a rms map of the spatial object

    Parameters
    ----------
    grid : spatial object to convert into a rms map
    path : path to save the figure if needed
    proj : projection of the map (Recommended: ccrs.PlateCarree(), ccrs.Mollweide())
    unit : unit of the grid in ['cmwe', 'mmwe', 'geoid', 'cmwe_ne', 'microGal', 'secacc']
    bound : list with minimal value and maximal value of the colorbar. Default value is None
    mask : np.array corresponding to the mask
    color : matplotlib cmap color of the gif (Recommended: viridis, plasma, RdBu_r, OrRd, Blues)
    """
    data_to_set = np.sqrt(np.sum(grid.data ** 2, axis=2) / grid.time.shape[0])

    if mask is not None:
        data_to_set *= mask

    plt.figure()
    matplotlib.rcParams['animation.embed_limit'] = 2 ** 128

    fig, ax1 = plt.subplots(num=1, nrows=1, ncols=1, figsize=(10.375, 6.625),
                            subplot_kw=dict(projection=proj))

    if bound is None:
        vmin, vmax = int(np.min(data_to_set)), int(np.ceil(np.max(data_to_set)))
    else:
        vmin, vmax = bound

    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.cm.get_cmap(color)
    im = ax1.imshow(data_to_set, interpolation='nearest',
                    norm=norm, cmap=cmap, transform=ccrs.PlateCarree(),
                    extent=grid.extent, origin='upper')
    ax1.coastlines('50m')

    # Add horizontal colorbar and adjust size
    # extend = add extension triangles to upper and lower bounds
    # options: neither, both, min, max
    # pad = distance from main plot axis
    # shrink = percent size of colorbar
    # aspect = lengthXwidth aspect of colorbar
    cbar = plt.colorbar(im, ax=ax1, extend='both', extendfrac=0.0375,
                        orientation='horizontal', pad=0.025, shrink=0.85,
                        aspect=22, drawedges=False)
    # rasterized colorbar to remove lines
    cbar.solids.set_rasterized(True)
    # Add label to the colorbar
    if unit == "cmwe":
        cbar.ax.set_xlabel('Equivalent Water Thickness', labelpad=10, fontsize=24)
        cbar.ax.set_ylabel('cm', fontsize=24, rotation=0, labelpad=10)
    elif unit == "cmwe_ne":
        cbar.ax.set_xlabel('Non elastic Equivalent Water Thickness', labelpad=10, fontsize=24)
        cbar.ax.set_ylabel('cm', fontsize=24, rotation=0, labelpad=10)
    elif unit == "mmwe":
        cbar.ax.set_xlabel('Equivalent Water Thickness', labelpad=10, fontsize=24)
        cbar.ax.set_ylabel('mm', fontsize=24, rotation=0, labelpad=10)
    elif unit == "geoid":
        cbar.ax.set_xlabel('Geoid Height', labelpad=10, fontsize=24)
        cbar.ax.set_ylabel('mm', fontsize=24, rotation=0, labelpad=10)
    elif unit == "microGal":
        cbar.ax.set_xlabel('Acceleration', labelpad=10, fontsize=24)
        cbar.ax.set_ylabel('$\mu Gal$', fontsize=24, rotation=0, labelpad=10)
    elif unit == "secacc":
        cbar.ax.set_xlabel('Secular Acceleration', labelpad=10, fontsize=24)
        cbar.ax.set_ylabel('$nT.y^{-2}$', fontsize=24, rotation=0, labelpad=10)

    cbar.ax.yaxis.set_label_coords(1.1, -0.4)
    # Set the tick levels for the colorbar
    # ticks lines all the way across
    cbar.ax.tick_params(which='both', width=1, length=26, labelsize=24,
                        direction='in')

    # stronger linewidth on frame
    ax1.spines['geo'].set_linewidth(2.0)
    ax1.spines['geo'].set_capstyle('projecting')
    # adjust subplot within figure
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.05, top=0.98)

    if path:
        plt.savefig(path, bbox_inches='tight')
    else:
        plt.show()
    plt.close()


def calc_rms_grid(grid, mask=None):
    """
    Compute Root Mean Square (RMS) value of a spatial object

    Parameters
    ----------
    grid : spatial object
    mask : mask to applied before rms computation

    Returns
    -------
    rms : rms of the grid

    """
    if mask is None:
        rms = np.sqrt(np.sum(
            [np.sum(np.cos(lat * np.pi / 180) ** 2 * line ** 2) for lat, line in zip(grid.lat, grid.data)]) / np.sum(
            [np.sum(np.cos(lat * np.pi / 180) ** 2 * line.size) for lat, line in zip(grid.lat, grid.data)]))

    else:
        rms = np.sqrt(np.sum(
            [np.sum(np.cos(lat * np.pi / 180) ** 2 * (line * np.swapaxes(np.tile(line_mask, (line.shape[1], 1)), 0, 1)) ** 2)
                              for lat, line, line_mask in zip(grid.lat, grid.data, mask)]) / np.sum(
            [np.sum(np.cos(lat * np.pi / 180) ** 2 * np.tile(line_mask, (line.shape[1], 1)))
                              for lat, line, line_mask in zip(grid.lat, grid.data, mask)]))
    # attention au cut dans les deux listes
    return rms


def plot_rms_grid(grid, path=False, labels=None, mask=None, unit='cmwe'):
    """
    Create a figure with rms of the grid spatial object in function of time

    Parameters
    ----------
    grid : spatial object or list of spatial object
    path : path to save the figure if needed
    mask : mask to apply on data if needed
    unit : unit of the grid
    """
    if type(grid) != list:
        grid = [grid]

    plot_rms = []
    for g in grid:
        l_rms = []
        for i in range(len(g.time)):
            if mask is None:
                rms = np.sqrt(np.sum([np.sum(np.cos(lat * np.pi / 180) ** 2 * line ** 2) for lat, line in
                                      zip(g.lat, g.data[:, :, i])]) / np.sum(
                    [np.sum(np.cos(lat * np.pi / 180) ** 2 * line.size) for lat, line in
                     zip(g.lat, g.data[:, :, i])]))
            else:
                rms = np.sqrt(np.sum([np.sum(np.cos(lat * np.pi / 180) ** 2 * (line * line_mask) ** 2)
                                      for lat, line, line_mask in zip(g.lat, g.data[:, :, i], mask)])
                              / np.sum([np.sum(np.cos(lat * np.pi / 180) ** 2 * line_mask)
                                        for lat, line_mask in zip(g.lat, mask)]))

            l_rms.append(rms)
        plot_rms.append(l_rms)

    plt.figure()
    if not(type(labels) == list):
        for g, rms in zip(grid, plot_rms):
            plt.plot(g.time, rms)
    else:
        for g, rms, l in zip(grid, plot_rms, labels):
            plt.plot(g.time, rms, label=l)

    plt.xlabel('Time (y)')
    if unit == "cmwe":
        plt.ylabel('cm EWH')
    elif unit == "cmwe_ne":
        plt.ylabel('Non elastic cm EWH')
    elif unit == "mmwe":
        plt.ylabel('mm EWH')
    elif unit == "geoid":
        plt.ylabel('mm Geoid Height')
    elif unit == "microGal":
        plt.ylabel('$\mu Gal$')
    elif unit == "secacc":
        plt.ylabel('$nT.y^{-2}$')
    plt.ylabel('Power (cm EWH)')

    if type(labels) == list:
        plt.legend()

    if path:
        plt.savefig(path, bbox_inches='tight')
    else:
        plt.show()
    plt.close()