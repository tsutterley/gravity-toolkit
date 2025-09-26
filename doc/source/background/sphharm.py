import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import gravity_toolkit as gravtk

# latitude and longitude
dlon, dlat = 0.625, 0.5
lat = np.arange(-90 + dlat/2.0, 90 + dlat/2.0, dlat)
lon = np.arange(0 + dlon/2.0, 360 + dlon/2.0, dlon)
gridlon, gridlat = np.meshgrid(lon, lat)
nlat, nlon = gridlat.shape
# colatitude and longitude in radians
theta = np.radians(90.0 - gridlat)
phi = np.radians(gridlon)

# minimum and maximum degree of spherical harmonics
lmin, lmax = (1, 4)
# number of rows and columns for subplots
nrows = (lmax - lmin) + 1
ncols = 2*lmax + 1
# compute associated Legendre functions 
Plm, dPlm = gravtk.associated_legendre(lmax, np.cos(theta))
# reshape to [l,m,lat,lon]
Plm = Plm.reshape((lmax+1, lmax+1, nlat, nlon))

# projection for the plots
projection = ccrs.Orthographic(central_longitude=0.0, central_latitude=0.0)

# plot spherical harmonics
fig = plt.figure(num=1, figsize=(12,7))
patch = mpatches.Rectangle((0, 0), 0.445, 1, color='0.975',
    zorder=0, transform=fig.transFigure)
fig.add_artist(patch)

for n, l in enumerate(range(lmin, lmax+1)):
    for m in range(-l, l+1):
        # setup subplot
        i = n*ncols + l + (lmax-l) + m + 1
        ax = fig.add_subplot(nrows, ncols, i, projection=projection)
        # spherical harmonics of degree l and order m
        Ylms = Plm[l,np.abs(m),:,:]*np.exp(1j*m*phi)
        Ylm = Ylms.imag if (m < 0) else Ylms.real
        # plot the surface
        ax.pcolormesh(lon, lat, Ylm,
            transform=ccrs.PlateCarree(),
            cmap='viridis', rasterized=True)
        # set the title
        ax.set_title(f'$l={l}, m={m}$')
        # add coastlines and set global
        ax.coastlines()
        ax.set_global()
        # turn off the axis
        ax.set_axis_off()

# add labels for cosine and sine terms
t1 = fig.text(0.05, 0.925, '$S_{lm}$', size=20,
    ha="center", va="center", transform=fig.transFigure)
t2 = fig.text(0.95, 0.925, '$C_{lm}$', size=20,
    ha="center", va="center", transform=fig.transFigure)

# adjust spacing and show
plt.tight_layout()
plt.show()