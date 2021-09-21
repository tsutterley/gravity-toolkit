#!/usr/bin/env python
u"""
tools.py
Written by Tyler Sutterley (09/2021)
Jupyter notebook, user interface and plotting tools

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    tkinter: Python interface to the Tcl/Tk GUI toolkit
        https://docs.python.org/3/library/tkinter.html
    ipywidgets: interactive HTML widgets for Jupyter notebooks and IPython
        https://ipywidgets.readthedocs.io/en/latest/
    matplotlib: Python 2D plotting library
        http://matplotlib.org/
        https://github.com/matplotlib/matplotlib

PROGRAM DEPENDENCIES:
    grace_find_months.py: finds available months for a GRACE/GRACE-FO dataset
    grace_date.py: reads GRACE index file and calculates dates for each month

UPDATE HISTORY:
    Written 09/2021
"""
import os
import re
import copy
import colorsys
import ipywidgets
import IPython.display
import tkinter.filedialog
import matplotlib.cm as cm
import matplotlib.colors as colors
from gravity_toolkit.grace_find_months import grace_find_months

class widgets:
    def __init__(self, **kwargs):
        """
        Widgets for setting directory and updating local data repository
        """
        # set default keyword arguments
        kwargs.setdefault('directory', os.getcwd())
        kwargs.setdefault('defaults', ['CSR','RL06','GSM',60])
        kwargs.setdefault('style', {})
        # set style
        self.style = copy.copy(kwargs['style'])
        # set the directory with GRACE/GRACE-FO data
        self.directory = ipywidgets.Text(
            value=kwargs['directory'],
            description='Directory:',
            disabled=False,
            style=self.style,
        )
        # button and label for directory selection
        self.directory_button = ipywidgets.Button(
            description="Directory select",
            mustexist="False",
            width="30%",
        )
        # connect directory select button with action
        self.directory_button.on_click(self.select_directory)
        # update local data with PO.DAAC https servers
        self.update = ipywidgets.Checkbox(
            value=True,
            description='Update data?',
            disabled=False,
            style=self.style,
        )
        # default parameters
        self.defaults = copy.copy(kwargs['defaults'])

    def select_directory(self, b):
        """function for directory selection
        """
        IPython.display.clear_output()
        root = tkinter.Tk()
        root.withdraw()
        root.call('wm', 'attributes', '.', '-topmost', True)
        b.directory = tkinter.filedialog.askdirectory()
        self.directory.value = copy.copy(b.directory)

    def select_product(self):
        """
        Widgets for setting specific data product and months
        """
        # dropdown menu for setting processing center
        # CSR: University of Texas Center for Space Research
        # GFZ: German Research Centre for Geosciences (GeoForschungsZentrum)
        # JPL: Jet Propulsion Laboratory
        # CNES: French Centre National D'Etudes Spatiales
        self.center = ipywidgets.Dropdown(
            options=['CSR', 'GFZ', 'JPL', 'CNES'],
            value=self.defaults[0],
            description='Center:',
            disabled=False,
            style=self.style,
        )

        # dropdown menu for setting data release
        self.release = ipywidgets.Dropdown(
            description='Release:',
            options=['RL04', 'RL05', 'RL06'],
            value=self.defaults[1],
            disabled=False,
            style=self.style,
        )

        # dropdown menu for setting data product
        # GAA: non-tidal atmospheric correction
        # GAB: non-tidal oceanic correction
        # GAC: combined non-tidal atmospheric and oceanic correction
        # GAD: GRACE/GRACE-FO ocean bottom pressure product
        # GSM: corrected monthly GRACE/GRACE-FO static field product
        self.product = ipywidgets.Dropdown(
            description='Product:',
            options=['GAC', 'GAD', 'GSM'],
            value=self.defaults[2],
            disabled=False,
            style=self.style,
        )

        # find available months for data product
        total_months = grace_find_months(self.base_directory,
            self.center.value, self.release.value,
            DSET=self.product.value)
        # select months to run
        # https://tsutterley.github.io/data/GRACE-Months.html
        options=[str(m).zfill(3) for m in total_months['months']]
        self.months = ipywidgets.SelectMultiple(
            options=options,
            value=options,
            description='Months:',
            disabled=False,
            style=self.style,
        )

        # watch widgets for changes
        self.center.observe(self.set_release)
        self.release.observe(self.set_product)
        self.center.observe(self.set_product)
        self.center.observe(self.update_months)
        self.release.observe(self.update_months)


    # function for setting the data release
    def set_release(self, sender):
        """function for updating available releases
        """
        if (self.center.value == 'CNES'):
            releases = ['RL01','RL02','RL03', 'RL04', 'RL05']
        else:
            releases = ['RL04', 'RL05', 'RL06']
        self.release.options=releases
        self.release.value=releases[-1]

    # function for setting the data product
    def set_product(self, sender):
        """function for updating available products
        """
        if (self.center.value == 'CNES'):
            products = {}
            products['RL01'] = ['GAC', 'GSM']
            products['RL02'] = ['GAA', 'GAB', 'GSM']
            products['RL03'] = ['GAA', 'GAB', 'GSM']
            products['RL04'] = ['GSM']
            products['RL05'] = ['GAA', 'GAB', 'GSM']
            valid_products = products[self.release.value]
        elif (self.center.value == 'CSR'):
            valid_products = ['GAC', 'GAD', 'GSM']
        elif (self.center.value in ('GFZ','JPL')):
            valid_products = ['GAA', 'GAB', 'GAC', 'GAD', 'GSM']
        self.product.options=valid_products
        self.product.value=self.defaults[2]

    # function for updating the available months
    def update_months(self, sender):
        """function for updating available months
        """
        # https://tsutterley.github.io/data/GRACE-Months.html
        total_months = grace_find_months(self.base_directory,
            self.center.value, self.release.value,
            DSET=self.product.value)
        options=[str(m).zfill(3) for m in total_months['months']]
        self.months.options=options
        self.months.value=options

    def select_options(self, **kwargs):
        """
        Widgets for setting data truncation and harmonic replacements
        """
        # set default keyword arguments

        # set the spherical harmonic truncation parameters
        # text entry for spherical harmonic degree
        self.lmax = ipywidgets.BoundedIntText(
            min=0,
            max=self.defaults[3],
            value=self.defaults[3],
            step=1,
            description='<i>&#8467;</i><sub>max</sub>:',
            disabled=False,
            style=self.style,
        )

        # text entry for spherical harmonic order
        self.mmax = ipywidgets.BoundedIntText(
            min=0,
            max=self.defaults[3],
            value=self.defaults[3],
            step=1,
            description='<i>m</i><sub>max</sub>:',
            disabled=False,
            style=self.style,
        )

        # dropdown menu for setting geocenter
        # Tellus: GRACE/GRACE-FO TN-13 from PO.DAAC
        #    https://grace.jpl.nasa.gov/data/get-data/geocenter/
        # SLR: satellite laser ranging from CSR
        #    ftp://ftp.csr.utexas.edu/pub/slr/geocenter/
        # SLF: Sutterley and Velicogna, Remote Sensing (2019)
        #    https://www.mdpi.com/2072-4292/11/18/2108
        geocenter_default = 'SLF' if (self.product.value == 'GSM') else '[none]'
        self.geocenter = ipywidgets.Dropdown(
            options=['[none]', 'Tellus', 'SLR', 'SLF'],
            value=geocenter_default,
            description='Geocenter:',
            disabled=False,
            style=self.style,
        )

        # SLR C20
        C20_default = 'GSFC' if (self.product.value == 'GSM') else '[none]'
        self.C20 = ipywidgets.Dropdown(
            options=['[none]','CSR','GSFC'],
            value=C20_default,
            description='SLR C20:',
            disabled=False,
            style=self.style,
        )

        # SLR C21 and S21
        self.CS21 = ipywidgets.Dropdown(
            options=['[none]','CSR'],
            value='[none]',
            description='SLR CS21:',
            disabled=False,
            style=self.style,
        )

        # SLR C22 and S22
        self.CS22 = ipywidgets.Dropdown(
            options=['[none]','CSR'],
            value='[none]',
            description='SLR CS22:',
            disabled=False,
            style=self.style,
        )

        # SLR C30
        C30_default = 'GSFC' if (self.product.value == 'GSM') else '[none]'
        self.C30 = ipywidgets.Dropdown(
            options=['[none]','CSR','GSFC'],
            value=C30_default,
            description='SLR C30:',
            disabled=False,
            style=self.style,
        )

        # SLR C50
        self.C50 = ipywidgets.Dropdown(
            options=['[none]','CSR','GSFC'],
            value='[none]',
            description='SLR C50:',
            disabled=False,
            style=self.style,
        )

        # Pole Tide Drift (Wahr et al., 2015) for Release-5
        poletide_default = True if ((self.release.value == 'RL05')
            and (self.product.value == 'GSM')) else False
        self.pole_tide = ipywidgets.Checkbox(
            value=poletide_default,
            description='Pole Tide Corrections',
            disabled=False,
            style=self.style,
        )

        # ECMWF Atmospheric Jump Corrections for Release-5
        atm_default = True if (self.release.value == 'RL05') else False
        self.atm = ipywidgets.Checkbox(
            value=atm_default,
            description='ATM Corrections',
            disabled=False,
            style=self.style,
        )

        # watch processing center widget for changes
        self.product.observe(self.set_max_degree)
        # watch data release widget for changes
        self.release.observe(self.set_max_degree)
        self.release.observe(self.set_pole_tide)
        self.release.observe(self.set_atm_corr)
        # watch data product widget for changes
        self.release.observe(self.set_pole_tide)
        # watch spherical harmonic degree widget for changes
        self.lmax.observe(self.set_max_order)

    # function for setting the spherical harmonic degree
    def set_max_degree(self, sender):
        """function for setting max degree of a product
        """
        if (self.center == 'CNES'):
            LMAX = dict(RL01=50,RL02=50,RL03=80,RL04=90,RL05=90)
        elif (self.center in ('CSR','JPL')):
            # CSR RL04/5/6 at LMAX 60
            # JPL RL04/5/6 at LMAX 60
            LMAX = dict(RL04=60,RL05=60,RL06=60)
        elif (self.center == 'GFZ'):
            # GFZ RL04/5 at LMAX 90
            # GFZ RL06 at LMAX 60
            LMAX = dict(RL04=90,RL05=90,RL06=60)
        self.lmax.max=LMAX[self.release.value]
        self.lmax.value=LMAX[self.release.value]

    # function for setting the spherical harmonic order
    def set_max_order(self, sender):
        """function for setting default max order
        """
        self.mmax.max=self.lmax.value
        self.mmax.value=self.lmax.value

    # function for setting pole tide drift corrections for Release-5
    def set_pole_tide(self, sender):
        """function for setting default pole tide correction for a release
        """
        self.pole_tide.value = True if ((self.release.value == 'RL05')
            and (self.product.value == 'GSM')) else False

    # function for setting atmospheric jump corrections for Release-5
    def set_atm_corr(self, sender):
        """function for setting default ATM correction for a release
        """
        self.atm.value = True if (self.release.value == 'RL05') else False

    def select_corrections(self, **kwargs):
        """
        Widgets for setting data corrections and processing
        """
        # set default keyword arguments
        kwargs.setdefault('units', ['cmwe','mmGH','mmCU',u'\u03BCGal','mbar'])

        # set the GIA file
        # files come in different formats depending on the group
        self.GIA_file = ipywidgets.Text(
            value='',
            description='GIA File:',
            disabled=False,
            style=self.style,
        )
        # button and label for input file selection
        self.GIA_button = ipywidgets.Button(
            description="File select",
            width="30%",
        )
        # connect fileselect button with action
        self.GIA_button.on_click(self.select_GIA_file)

        # dropdown menu for setting GIA model
        # IJ05-R2: Ivins R2 GIA Models
        # W12a: Whitehouse GIA Models
        # SM09: Simpson/Milne GIA Models
        # ICE6G: ICE-6G GIA Models
        # Wu10: Wu (2010) GIA Correction
        # AW13-ICE6G: Geruo A ICE-6G GIA Models
        # Caron: Caron JPL GIA Assimilation
        # ICE6G-D: ICE-6G Version-D GIA Models
        # ascii: GIA reformatted to ascii
        # netCDF4: GIA reformatted to netCDF4
        # HDF5: GIA reformatted to HDF5
        gia_list = ['[None]','IJ05-R2','W12a','SM09','ICE6G',
            'Wu10','AW13-ICE6G','Caron','ICE6G-D',
            'ascii','netCDF4','HDF5']
        self.GIA = ipywidgets.Dropdown(
            options=gia_list,
            value='[None]',
            description='GIA Type:',
            disabled=False,
            style=self.style,
        )

        # set the files to be removed
        self.remove_files = []
        self.remove_file = ipywidgets.Text(
            value='',
            description='Rem. Files:',
            disabled=False,
            style=self.style,
        )

        # button and label for input file selection
        self.remove_button = ipywidgets.Button(
            description="File select",
        )
        # connect fileselect button with action
        self.remove_button.on_click(self.select_remove_file)
        self.remove_file.observe(self.set_removefile)

        # dropdown menu for setting remove file type
        # netCDF4: single netCDF4 file
        # HDF5: single HDF5 file
        # index (ascii): index of monthly ascii files
        # index (netCDF4): index of monthly netCDF4 files
        # index (HDF5): index of monthly HDF5 files
        remove_list = ['[None]','netCDF4','HDF5',
            'index (ascii)','index (netCDF4)','index (HDF5)']
        self.remove_format = ipywidgets.Dropdown(
            options=remove_list,
            value='[None]',
            description='Rem. Type:',
            disabled=False,
            style=self.style,
        )

        # redestribute removed file mass over the ocean
        self.redistribute_removed = ipywidgets.Checkbox(
            value=False,
            description='Redistribute Removed',
            disabled=False,
            style=self.style,
        )

        # path to land-sea mask for ocean redistribution
        self.mask = ipywidgets.Text(
            value='',
            description='Mask File:',
            disabled=False,
            style=self.style,
        )
        # button and label for input file selection
        self.mask_button = ipywidgets.Button(
            description="File select",
            width="30%",
        )
        # connect fileselect button with action
        self.mask_button.on_click(self.select_mask_file)

        # text entry for Gaussian Smoothing Radius in km
        self.gaussian = ipywidgets.BoundedFloatText(
            value=300,
            min=0,
            max=1000.0,
            step=50,
            description='Gaussian:',
            disabled=False,
            style=self.style,
        )

        # Destripe Spherical Harmonics
        self.destripe = ipywidgets.Checkbox(
            value=True,
            description='Destripe',
            disabled=False,
            style=self.style,
        )

        # text entry for output spatial degree spacing
        self.spacing = ipywidgets.BoundedFloatText(
            value=1.0,
            min=0,
            max=360.0,
            step=0.5,
            description='Spacing:',
            disabled=False,
            style=self.style,
        )

        # dropdown menu for setting output degree interval
        interval_list = ['(-180:180,90:-90)', '(Degree spacing)/2']
        self.interval = ipywidgets.Dropdown(
            options=interval_list,
            value='(Degree spacing)/2',
            description='Interval:',
            disabled=False,
            style=self.style,
        )

        # dropdown menu for setting units
        # 1: cm of water thickness
        # 2: mm of geoid height
        # 3: mm of elastic crustal deformation
        # 4: microGal gravitational perturbation
        # 5: millibar of equivalent surface pressure
        self.units = ipywidgets.Dropdown(
            options=kwargs['units'],
            value='cmwe',
            description='Units:',
            disabled=False,
            style=self.style,
        )

    def select_GIA_file(self, b):
        """function for GIA file selection
        """
        IPython.display.clear_output()
        root = tkinter.Tk()
        root.withdraw()
        root.call('wm', 'attributes', '.', '-topmost', True)
        filetypes = (("All Files", "*.*"))
        b.files = tkinter.filedialog.askopenfilename(
            filetypes=filetypes,
            multiple=False)
        self.GIA_file.value = copy.copy(b.files)

    def select_remove_file(self, b):
        """function for removed file selection
        """
        IPython.display.clear_output()
        root = tkinter.Tk()
        root.withdraw()
        root.call('wm', 'attributes', '.', '-topmost', True)
        filetypes = (("ascii file", "*.txt"),
            ("HDF5 file", "*.h5"),
            ("netCDF file", "*.nc"),
            ("All Files", "*.*"))
        b.files = tkinter.filedialog.askopenfilename(
            defaultextension='nc',
            filetypes=filetypes,
            multiple=True)
        self.remove_files.extend(b.files)
        self.set_removelabel()

    def set_removefile(self, sender):
        """function for updating removed file list
        """
        if self.remove_file.value:
            self.remove_files = self.remove_file.value.split(',')
        else:
            self.remove_files = []

    def set_removelabel(self):
        """function for updating removed file label
        """
        self.remove_file.value = ','.join(self.remove_files)

    def select_mask_file(self, b):
        """function for mask file selection
        """
        IPython.display.clear_output()
        root = tkinter.Tk()
        root.withdraw()
        root.call('wm', 'attributes', '.', '-topmost', True)
        filetypes = (("netCDF file", "*.nc"),
            ("All Files", "*.*"))
        b.files = tkinter.filedialog.askopenfilename(
            defaultextension='nc',
            filetypes=filetypes,
            multiple=False)
        self.mask.value = copy.copy(b.files)

    def select_output(self, **kwargs):
        """
        Widget for setting output data file format
        """
        # set default keyword arguments
        # dropdown menu for setting output data format
        self.output_format = ipywidgets.Dropdown(
            options=['[None]','netCDF4', 'HDF5'],
            value='[None]',
            description='Output:',
            disabled=False,
            style=self.style,
        )

    @property
    def base_directory(self):
        return os.path.expanduser(self.directory.value)

    @property
    def landmask(self):
        return os.path.expanduser(self.mask.value)

    @property
    def unit_index(self):
        return self.units.index + 1

    @property
    def format(self):
        return self.output_format.value

class colormap:
    def __init__(self, **kwargs):
        # set default keyword arguments
        kwargs.setdefault('vmin', None)
        kwargs.setdefault('vmax', None)
        kwargs.setdefault('steps', 10)
        kwargs.setdefault('cmaps_listed', {})
        kwargs.setdefault('style', {})
        # set style
        self.style = copy.copy(kwargs['style'])
        # vmin and vmax for normalization
        self.vmin = copy.copy(kwargs['vmin'])
        self.vmax = copy.copy(kwargs['vmax'])
        # slider for range of color bar
        self.range = ipywidgets.IntRangeSlider(
            value=[self.vmin,self.vmax],
            min=self.vmin,
            max=self.vmax,
            step=1,
            description='Plot Range:',
            disabled=False,
            continuous_update=False,
            orientation='horizontal',
            readout=True,
            style=self.style,
        )

        # slider for steps in color bar
        step = (self.vmax-self.vmin)//kwargs['steps']
        self.step = ipywidgets.IntSlider(
            value=step,
            min=0,
            max=self.vmax-self.vmin,
            step=1,
            description='Plot Step:',
            disabled=False,
            continuous_update=False,
            orientation='horizontal',
            readout=True,
            style=self.style,
        )

        # all listed colormaps in matplotlib version
        cmap_set = set(cm.datad.keys()) | set(cm.cmaps_listed.keys())
        # colormaps available in this program
        # (no reversed, qualitative or miscellaneous)
        self.cmaps_listed = copy.copy(kwargs['cmaps_listed'])
        self.cmaps_listed['Perceptually Uniform Sequential'] = [
            'viridis','plasma','inferno','magma','cividis']
        self.cmaps_listed['Sequential'] = ['Greys','Purples',
            'Blues','Greens','Oranges','Reds','YlOrBr','YlOrRd',
            'OrRd','PuRd','RdPu','BuPu','GnBu','PuBu','YlGnBu',
            'PuBuGn','BuGn','YlGn']
        self.cmaps_listed['Sequential (2)'] = ['binary','gist_yarg',
            'gist_gray','gray','bone','pink','spring','summer',
            'autumn','winter','cool','Wistia','hot','afmhot',
            'gist_heat','copper']
        self.cmaps_listed['Diverging'] = ['PiYG','PRGn','BrBG',
            'PuOr','RdGy','RdBu','RdYlBu','RdYlGn','Spectral',
            'coolwarm', 'bwr','seismic']
        self.cmaps_listed['Cyclic'] = ['twilight',
            'twilight_shifted','hsv']
        # create list of available colormaps in program
        cmap_list = []
        for val in self.cmaps_listed.values():
            cmap_list.extend(val)
        # reduce colormaps to available in program and matplotlib
        cmap_set &= set(cmap_list)
        # dropdown menu for setting colormap
        self.name = ipywidgets.Dropdown(
            options=sorted(cmap_set),
            value='viridis',
            description='Colormap:',
            disabled=False,
            style=self.style,
        )

        # Reverse the colormap
        self.reverse = ipywidgets.Checkbox(
            value=False,
            description='Reverse Colormap',
            disabled=False,
            style=self.style,
        )

    @property
    def _r(self):
        cmap_reverse_flag = '_r' if self.reverse.value else ''
        return cmap_reverse_flag

    @property
    def value(self):
        return copy.copy(cm.get_cmap(self.name.value + self._r))

    @property
    def norm(self):
        cmin,cmax = self.range.value
        return colors.Normalize(vmin=cmin,vmax=cmax)

    @property
    def levels(self):
        cmin,cmax = self.range.value
        return [l for l in range(cmin,cmax+self.step.value,self.step.value)]

    @property
    def label(self):
        return ['{0:0.0f}'.format(ct) for ct in self.levels]

def from_cpt(filename, use_extremes=True, **kwargs):
    """
    Reads a GMT color palette table
    Can import HSV (hue-saturation-value) or RGB values

    Keyword Arguments
    -----------------
    use_extremes: use the under, over and bad values from the cpt file
    """

    # read the cpt file and get contents
    with open(filename,'r') as f:
        file_contents = f.read().splitlines()
    # extract basename from cpt filename
    name = re.sub(r'\.cpt','',os.path.basename(filename),flags=re.I)

    # compile regular expression operator to find numerical instances
    rx = re.compile(r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?')

    # create list objects for x, r, g, b
    x,r,g,b = ([],[],[],[])
    # assume RGB color model
    colorModel = "RGB"
    # back, forward and no data flags
    flags = dict(B=None,F=None,N=None)
    for line in file_contents:
        # find back, forward and no-data flags
        model = re.search(r'COLOR_MODEL.*(HSV|RGB)',line,re.I)
        BFN = re.match(r'[BFN]',line,re.I)
        # parse non-color data lines
        if model:
            # find color model
            colorModel = model.group(1)
            continue
        elif BFN:
            flags[BFN.group(0)] = [float(i) for i in rx.findall(line)]
            continue
        elif re.search(r"#",line):
            # skip over commented header text
            continue
        # find numerical instances within line
        x1,r1,g1,b1,x2,r2,g2,b2 = rx.findall(line)
        # append colors and locations to lists
        x.append(float(x1))
        r.append(float(r1))
        g.append(float(g1))
        b.append(float(b1))
    # append end colors and locations to lists
    x.append(float(x2))
    r.append(float(r2))
    g.append(float(g2))
    b.append(float(b2))

    # convert input colormap to output
    xNorm = [None]*len(x)
    if (colorModel == "HSV"):
        # convert HSV (hue-saturation-value) to RGB
        # calculate normalized locations (0:1)
        for i,xi in enumerate(x):
            rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
            r[i] = rr
            g[i] = gg
            b[i] = bb
            xNorm[i] = (xi - x[0])/(x[-1] - x[0])
    elif (colorModel == "RGB"):
        # normalize hexadecimal RGB triple from (0:255) to (0:1)
        # calculate normalized locations (0:1)
        for i,xi in enumerate(x):
            r[i] /= 255.0
            g[i] /= 255.0
            b[i] /= 255.0
            xNorm[i] = (xi - x[0])/(x[-1] - x[0])

    # output RGB lists containing normalized location and colors
    cdict = dict(red=[None]*len(x),green=[None]*len(x),blue=[None]*len(x))
    for i,xi in enumerate(x):
        cdict['red'][i] = [xNorm[i],r[i],r[i]]
        cdict['green'][i] = [xNorm[i],g[i],g[i]]
        cdict['blue'][i] = [xNorm[i],b[i],b[i]]

    # create colormap for use in matplotlib
    cmap = colors.LinearSegmentedColormap(name, cdict, **kwargs)
    # set flags for under, over and bad values
    extremes = dict(under=None,over=None,bad=None)
    for key,attr in zip(['B','F','N'],['under','over','bad']):
        if flags[key] is not None:
            r,g,b = flags[key]
            if (colorModel == "HSV"):
                # convert HSV (hue-saturation-value) to RGB
                r,g,b = colorsys.hsv_to_rgb(r/360.,g,b)
            elif (colorModel == 'RGB'):
                # normalize hexadecimal RGB triple from (0:255) to (0:1)
                r,g,b = (r/255.0,g/255.0,b/255.0)
            # set attribute for under, over and bad values
            extremes[attr] = (r,g,b)
    # create copy of colormap with extremes
    if use_extremes:
        cmap = cmap.with_extremes(**extremes)
    # register colormap to be recognizable by cm.get_cmap()
    cm.register_cmap(name=name, cmap=cmap)
    # return the colormap
    return cmap
