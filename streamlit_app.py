import streamlit as st
from sunpy import map as smap
from astropy.time import Time
from sunpy.util import MetaDict
import numpy as np
import astropy.constants as const
import astropy.units as u
from CoronalLoopBuilder.builder import CoronalLoopBuilder
import matplotlib.pyplot as plt
import io

st.title('Coronal Loop Builder')
st.write('Interactive visualization for coronal loops on solar images.')

# Sidebar for user inputs
st.sidebar.header('Loop Parameters')
radius = st.sidebar.slider('Radius (Mm)', min_value=10, max_value=100, value=50)
height = st.sidebar.slider('Height (Mm)', min_value=-50, max_value=50, value=-10)
phi0 = st.sidebar.slider('Phi0 (degrees)', min_value=0, max_value=360, value=45)
theta0 = st.sidebar.slider('Theta0 (degrees)', min_value=0, max_value=180, value=30)
el = st.sidebar.slider('El (degrees)', min_value=0, max_value=360, value=100)
az = st.sidebar.slider('Az (degrees)', min_value=-90, max_value=90, value=10)
samples_num = st.sidebar.slider('Samples Number', min_value=10, max_value=1000, value=100)

# Create dummy solar maps for visualization (similar to the example script)
time_now = Time.now()

# The following dummy map creation is repetitive for the two maps, could be functionized for brevity
data1 = np.ones((10, 10))
meta1 = MetaDict({
    'ctype1': 'HPLN-TAN', 'ctype2': 'HPLT-TAN',
    'cunit1': 'arcsec', 'cunit2': 'arcsec',
    'crpix1': (data1.shape[0] + 1) / 2., 'crpix2': (data1.shape[1] + 1) / 2.,
    'cdelt1': 1.0, 'cdelt2': 1.0, 'crval1': 0.0, 'crval2': 0.0,
    'hgln_obs': 0.0,  ## Stonyhurst heliographic longitude in degree
    'hglt_obs': 0.0,  ## Stonyhurst heliographic latitude in degree
    'dsun_obs': const.au.to(u.m).value, 'dsun_ref': const.au.to(u.m).value,
    'rsun_ref': const.R_sun.to(u.m).value,
    'rsun_obs': ((const.R_sun / const.au).decompose() * u.radian).to(u.arcsec).value,
    't_obs': time_now.iso, 'date-obs': time_now.iso,
})
dummy_map1 = smap.GenericMap(data1, meta1)

data2 = np.ones((10, 10))
meta2 = MetaDict({
    'ctype1': 'HPLN-TAN', 'ctype2': 'HPLT-TAN',
    'cunit1': 'arcsec', 'cunit2': 'arcsec',
    'crpix1': (data2.shape[0] + 1) / 2., 'crpix2': (data2.shape[1] + 1) / 2.,
    'cdelt1': 1.0, 'cdelt2': 1.0, 'crval1': 0.0, 'crval2': 0.0,
    'hgln_obs': 106.0,  ## Stonyhurst heliographic longitude in degree
    'hglt_obs': 5.0,  ## Stonyhurst heliographic latitude in degree
    'dsun_obs': const.au.to(u.m).value, 'dsun_ref': const.au.to(u.m).value,
    'rsun_ref': const.R_sun.to(u.m).value,
    'rsun_obs': 2 * ((const.R_sun / const.au).decompose() * u.radian).to(u.arcsec).value,
    't_obs': time_now.iso, 'date-obs': time_now.iso,
})
dummy_map2 = smap.GenericMap(data2, meta2)

dummy_maps = [dummy_map1, dummy_map2]
num_maps = len(dummy_maps)

# Visualization logic
fig = plt.figure(figsize=(6 * num_maps, 6))
axs = []
for midx, dummy_map in enumerate(dummy_maps):
    ax = fig.add_subplot(1, num_maps, midx + 1, projection=dummy_map)
    axs.append(ax)
    dummy_map.plot(alpha=0, extent=[-1000, 1000, -1000, 1000], title=False, axes=ax)
    dummy_map.draw_grid(axes=ax, grid_spacing=10 * u.deg, color='k')
    dummy_map.draw_limb(axes=ax, color='k')

# coronal_loop = CoronalLoopBuilder(
#     fig, axs, dummy_maps,
#     radius=radius * u.Mm, height=height * u.Mm, phi0=phi0 * u.deg, theta0=theta0 * u.deg,
#     el=el * u.deg, az=az * u.deg, samples_num=samples_num
# )

# Convert plot to image and display in Streamlit
buf = io.BytesIO()
plt.savefig(buf, format="png")
buf.seek(0)
st.image(buf)

# Additional functionality can be added, like saving the loop data or adjusting parameters further
