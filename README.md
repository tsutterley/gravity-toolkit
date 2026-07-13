# gravity-toolkit

Python tools for obtaining and working with Level-2 spherical harmonic coefficients from the NASA/DLR Gravity Recovery and Climate Experiment (GRACE) and the NASA/GFZ Gravity Recovery and Climate Experiment Follow-On (GRACE-FO) missions

## About

<table>
  <tr>
    <td><b>Version:</b></td>
    <td>
        <a href="https://pypi.python.org/pypi/gravity-toolkit/" alt="PyPI"><img src="https://img.shields.io/pypi/v/gravity-toolkit.svg"></a>
        <a href="https://anaconda.org/conda-forge/gravity-toolkit" alt="conda-forge"><img src="https://img.shields.io/conda/vn/conda-forge/gravity-toolkit"></a>
        <a href="https://github.com/tsutterley/gravity-toolkit/releases/latest" alt="commits-since"><img src="https://img.shields.io/github/commits-since/tsutterley/gravity-toolkit/latest"></a>
    </td>
  </tr>
  <tr>
    <td><b>Citation:</b></td>
    <td>
        <a href="https://doi.org/10.5281/zenodo.5156864" alt="zenodo"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.5156864.svg"></a>
    </td>
  </tr>
  <tr>
    <td><b>Tests:</b></td>
    <td>
        <a href="https://gravity-toolkit.readthedocs.io/en/latest/?badge=latest" alt="Documentation Status"><img src="https://readthedocs.org/projects/gravity-toolkit/badge/?version=latest"></a>
        <a href="https://github.com/tsutterley/gravity-toolkit/actions/workflows/python-request.yml" alt="Build"><img src="https://github.com/tsutterley/gravity-toolkit/actions/workflows/python-request.yml/badge.svg"></a>
        <a href="https://github.com/tsutterley/gravity-toolkit/actions/workflows/ruff-format.yml" alt="Ruff"><img src="https://github.com/tsutterley/gravity-toolkit/actions/workflows/ruff-format.yml/badge.svg"></a>
    </td>
  </tr>
  <tr>
    <td><b>Data:</b></td>
    <td>
        <a href="https://doi.org/10.6084/m9.figshare.7388540" alt="figshare"><img src="https://img.shields.io/badge/figshare-geocenter-a60845?logo=figshare"></a>
        <a href="https://doi.org/10.6084/m9.figshare.9702338" alt="figshare"><img src="https://img.shields.io/badge/figshare-sealevel-71cac5?logo=figshare"></a>
    </td>
  </tr>
  <tr>
    <td><b>License:</b></td>
    <td>
        <a href="https://github.com/tsutterley/gravity-toolkit/blob/main/LICENSE" alt="License"><img src="https://img.shields.io/github/license/tsutterley/gravity-toolkit"></a>
    </td>
  </tr>
</table>

For more information: see the documentation at [gravity-toolkit.readthedocs.io](https://gravity-toolkit.readthedocs.io/)

## Installation

From PyPI:

```bash
python3 -m pip install gravity-toolkit
```

To include all optional dependencies:

```bash
python3 -m pip install gravity-toolkit[all]
```

Using `conda` or `mamba` from conda-forge:

```bash
conda install -c conda-forge gravity-toolkit
```

```bash
mamba install -c conda-forge gravity-toolkit
```

Development version from GitHub:

```bash
python3 -m pip install git+https://github.com/tsutterley/gravity-toolkit.git
```

### Running with Pixi

Alternatively, you can use [Pixi](https://pixi.sh/) for a streamlined workspace environment:

1. Install Pixi following the [installation instructions](https://pixi.sh/latest/#installation)
2. Clone the project repository:

```bash
git clone https://github.com/tsutterley/gravity-toolkit.git
```

3. Move into the `gravity-toolkit` directory

```bash
cd gravity-toolkit
```

4. Install dependencies and start JupyterLab:

```bash
pixi run start
```

This will automatically create the environment, install all dependencies, and launch JupyterLab in the [notebooks](./doc/source/notebooks/) directory.

## Resources

- [NASA GRACE mission site](https://www.nasa.gov/mission_pages/Grace/index.html)
- [NASA GRACE-FO mission site](https://www.nasa.gov/missions/grace-fo)
- [JPL GRACE Tellus site](https://grace.jpl.nasa.gov/)
- [JPL GRACE-FO site](https://gracefo.jpl.nasa.gov/)
- [UTCSR GRACE site](http://www.csr.utexas.edu/grace/)
- [GRACE at the NASA Physical Oceanography Distributed Active Archive Center (PO.DAAC)](https://podaac.jpl.nasa.gov/grace)
- [GRACE at the GFZ Information System and Data Center](http://isdc.gfz-potsdam.de/grace-isdc/)

## Dependencies

- [boto3: Amazon Web Services (AWS) SDK for Python](https://boto3.amazonaws.com/v1/documentation/api/latest/index.html)
- [future: Compatibility layer between Python 2 and Python 3](https://python-future.org/)
- [lxml: processing XML and HTML in Python](https://pypi.python.org/pypi/lxml)
- [matplotlib: Python 2D plotting library](https://matplotlib.org/)
- [netCDF4: Python interface to the netCDF C library](https://unidata.github.io/netcdf4-python/)
- [numpy: Scientific Computing Tools For Python](https://www.numpy.org)
- [platformdirs: Python module for determining platform-specific directories](https://pypi.org/project/platformdirs/)
- [python-dateutil: powerful extensions to datetime](https://dateutil.readthedocs.io/en/stable/)
- [PyYAML: YAML parser and emitter for Python](https://github.com/yaml/pyyaml)
- [scipy: Scientific Tools for Python](https://docs.scipy.org/doc/)

## Download

The program homepage is:  
<https://github.com/tsutterley/gravity-toolkit>

A zip archive of the latest version is available directly at:  
<https://github.com/tsutterley/gravity-toolkit/archive/main.zip>

## Disclaimer

This package includes software developed at the University of California at Irvine (UCI), the NASA Jet Propulsion Laboratory (JPL), NASA Goddard Space Flight Center (GSFC) and the University of Washington Applied Physics Laboratory (UW-APL).
This program is not sponsored or maintained by the Universities Space Research Association (USRA),
the Center for Space Research at the University of Texas (UTCSR), the Jet Propulsion Laboratory (JPL),
the German Research Centre for Geosciences (GeoForschungsZentrum, GFZ) or NASA.
The software is provided here for your convenience but *with no guarantees whatsoever*.

## Contributing

This project contains work and contributions from the [scientific community](./CONTRIBUTORS.md).
If you would like to contribute to the project, please have a look at the [contribution guidelines](./doc/source/getting_started/Contributing.rst), [open issues](https://github.com/tsutterley/gravity-toolkit/issues) and [discussions board](https://github.com/tsutterley/gravity-toolkit/discussions).

## References

> T. C. Sutterley, I. Velicogna, and C.-W. Hsu,
> "Self-Consistent Ice Mass Balance and Regional Sea Level From Time-Variable Gravity",
> *Earth and Space Science*, 7, (2020).
> [doi: 10.1029/2019EA000860](https://doi.org/10.1029/2019EA000860)
> 
> T. C. Sutterley and I. Velicogna,
> "Improved estimates of geocenter variability from time-variable gravity and ocean model outputs", 
> *Remote Sensing*, 11(18), 2108, (2019).
> [doi: 10.3390/rs11182108](https://doi.org/10.3390/rs11182108)
> 
> J. Wahr, S. C. Swenson, and I. Velicogna,
> "Accuracy of GRACE mass estimates",
> *Geophysical Research Letters*, 33(6), L06401, (2006).
> [doi: 10.1029/2005GL025305](https://doi.org/10.1029/2005GL025305)
> 
> J. Wahr, M. Molenaar, and F. Bryan,
> "Time variability of the Earth's gravity field: Hydrological and oceanic effects and their possible > detection using GRACE",
> *Journal of Geophysical Research: Solid Earth*, 103(B12), (1998).
> [doi: 10.1029/98JB02844](https://doi.org/10.1029/98JB02844)
> 
> D. Han and J. Wahr,
> "The viscoelastic relaxation of a realistically stratified earth, and a further analysis of postglacial rebound",
> *Geophysical Journal International*, 120(2), (1995).
> [doi: 10.1111/j.1365-246X.1995.tb01819.x](https://doi.org/10.1111/j.1365-246X.1995.tb01819.x)

## Data Repositories

> T. C. Sutterley, I. Velicogna, and C.-W. Hsu,
> "Ice Mass and Regional Sea Level Estimates from Time-Variable Gravity", (2020).
> [doi: 10.6084/m9.figshare.9702338](https://doi.org/10.6084/m9.figshare.9702338)
> 
> T. C. Sutterley and I. Velicogna,
> "Geocenter Estimates from Time-Variable Gravity and Ocean Model Outputs", (2019).
> [doi: 10.6084/m9.figshare.7388540](https://doi.org/10.6084/m9.figshare.7388540)

## License

The content of this project is licensed under the [Creative Commons Attribution 4.0 Attribution license](https://creativecommons.org/licenses/by/4.0/) and the source code is licensed under the [MIT license](LICENSE).
