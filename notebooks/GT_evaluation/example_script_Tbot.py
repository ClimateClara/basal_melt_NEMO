"""
Try to make a simple script for bottom T
"""

import xarray as xr

inputpath = '/thredds/tgcc/store/burgardc/FORMATTED/'
file_T_all = xr.open_dataset(inputpath + 'opencav-presc02_18500101_18891231_1M_grid_T_varofint.nc')

mask_ocean = file_T_all['so'].isel(time=0).drop('time_counter') > 0
vert_diff_minus_all = (mask_ocean - mask_ocean.shift(deptht=-1)).isel(deptht=range(1,len(mask_ocean.deptht)))
bot_depth_all = (mask_ocean.deptht * vert_diff_minus_all).where(vert_diff_minus_all > 0).sum('deptht').astype('float')

bot_depth_all_wo0 = bot_depth_all.where(bot_depth_all != 0, 5.057600e-01)
Tbot = file_T_all['thetao'].sel(deptht=bot_depth_all_wo0).where(bot_depth_all > 0)

Tbot.drop('deptht').to_dataset(name='T_bottom')