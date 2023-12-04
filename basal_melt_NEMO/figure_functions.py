"""
These are scripts to analyse NEMO output more or less rapidly

"""

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import cmocean
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature
from cartopy.util import add_cyclic_point

def maps_compare3(lon,lat,ref,modif,name,cmap,region,legend,lat_lim=-50):
    
    f = plt.figure(figsize=(8.25*3,8.25))
    #f = plt.figure()
    #f.suptitle(str(time_in.values)[0:16],fontsize=22)

    if region == 'global':
        proj = ccrs.Mollweide(central_longitude=0)
        wrap_ref = ref
        wrap_lon = lon
        wrap_modif = modif
    elif region == 'Ant':
        proj = ccrs.SouthPolarStereo(central_longitude=0)
        wrap_ref, wrap_lon = ref, lon #add_cyclic_point(ref.values,coord=lon,axis=1)
        wrap_modif, wrap_lon = modif, lon #add_cyclic_point(modif.values,coord=lon,axis=1)

        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpl.path.Path(verts * radius + center)        

    #### REFERENCE
    ax1 = plt.subplot(1, 3, 1, projection=proj)
    abso0 = ax1.pcolormesh(wrap_lon,lat,wrap_ref,transform=ccrs.PlateCarree(),cmap=cmap,rasterized=True)
    ax1.coastlines(resolution='110m', linewidth=0.5)
    
    if region == 'global':
        ax1.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='black')
    elif region == 'Ant':
        ax1.set_extent([-180, 180, -90, -50], crs=ccrs.PlateCarree())
        ax1.set_boundary(circle, transform=ax1.transAxes)

    if legend=='yes':
        cbar = f.colorbar(abso0, ax=ax1, shrink=0.3,orientation='vertical',extend='both')
        cbar.set_label(name,rotation=90)  

    #### MODIFICATION
    ax2 = plt.subplot(1, 3, 2, projection=proj)
    ax2.coastlines(resolution='110m', linewidth=0.5)
    abso = ax2.pcolormesh(wrap_lon,lat,wrap_modif,transform=ccrs.PlateCarree(),cmap=cmap,rasterized=True)

    if region == 'global':
        ax2.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='black')
    elif region == 'Ant':
        ax2.set_extent([-180, 180, -90, -50], crs=ccrs.PlateCarree())
        ax2.set_boundary(circle, transform=ax2.transAxes)

    if legend=='yes':
        cbar = f.colorbar(abso0, ax=ax2, shrink=0.3,orientation='vertical',extend='both')
        cbar.set_label(name,rotation=90)    

    #### DIFFERENCE
    abs_lim = abs(modif - ref).quantile(0.99)
    ax3 = plt.subplot(1, 3, 3, projection=proj)
    ax3.coastlines(resolution='110m', linewidth=0.5)
    diff1 = ax3.pcolormesh(wrap_lon,lat,wrap_modif - wrap_ref,transform=ccrs.PlateCarree(),cmap=mpl.cm.coolwarm,rasterized=True,vmin=-abs_lim,vmax=abs_lim)
    
    if region == 'global':
        ax3.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='black')
    elif region == 'Ant':
        ax3.set_extent([-180, 180, -90, -50], crs=ccrs.PlateCarree())
        ax3.set_boundary(circle, transform=ax3.transAxes)

    if legend=='yes':
        cbar = f.colorbar(diff1, ax=ax3, shrink=0.3,orientation='vertical',extend='both')
        cbar.set_label('Modif - ref '+name,rotation=90)    


    
    f.tight_layout()
    
    return f
    