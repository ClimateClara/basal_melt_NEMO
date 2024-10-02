import xarray as xr
import numpy as np

import matplotlib.pyplot as plt

import seaborn as sns

from basal_melt_NEMO.constants import *
import basal_melt_NEMO.useful_functions as uf

def var_obs_mean_std(var_list_in):

    var_list = ['global_SST','wed_gyre','ross_gyre','ACC',
                'FRIS_melt','Ross_melt','LarsenC_melt','total_melt',
                'mar_sie_arc','sep_sie_arc','feb_sie_ant','sep_sie_ant',
                'mar_sia_arc','sep_sia_arc','feb_sia_ant','sep_sia_ant',
                'mar_siv_arc','sep_siv_arc','feb_siv_ant','sep_siv_ant',
                'Sbot_WWED','Sbot_EWED','Sbot_WROSS','Sbot_EROSS','Sbot_AMU',
                'Tbot_WWED','Tbot_EWED','Tbot_WROSS','Tbot_EROSS','Tbot_AMU',
               'AMOC_26N','AMOC_30S','OHC_tot','OHC_700','OHC_2000']
    
    var_obs_mean = xr.DataArray(data=np.array([
                20., 56.0, 20.0, 136.7, 
               155.4, 47.7, 20.7, 1500,
               15.3, 6.3, 3.3, 19.8, 
               14.4, 6.0, 2.25, 15.0,
               27., 9., 2.5, 20.,
               34.9, np.nan, 35.0, np.nan, np.nan, 
               np.nan, -1.95, -1.9, -1.7, np.nan,
               17.0, np.nan, np.nan, np.nan, np.nan]), dims='var').assign_coords({'var': var_list})

    var_obs_std = xr.DataArray(data=np.array([
                  np.nan, 8.0, 5.0, 7.8, 
                  45.0, 34.0, 67.0, 237,
                  1.0, 1.0, 0.5, 0.5, 
                  1.29, 1.49, 0.25, 0.4,
                  3, 3, 2, 5,
                  0.0, np.nan, 0.0, np.nan, np.nan,
                  np.nan, 0.2, 0.4, 0.4, np.nan,
                  10.0, np.nan, np.nan, np.nan, np.nan]), dims='var').assign_coords({'var': var_list})

    return var_obs_mean.sel(var=var_list_in), var_obs_std.sel(var=var_list_in)

def bottom_prop_reg(var,lon,lat,reg,bottom_cell):

    # mask regions for bottom properties
    mask_regions = xr.Dataset()
    mask_regions['AMU'] = (lon >= -109.64) & (lon <= -102.23) & (lat >= -75.80) & (lat <= -71.66)
    mask_regions['WROSS'] = (lon >= 157.100) & (lon <= 173.333) & (lat >= -78.130) & (lat <= -74.040)
    mask_regions['EROSS'] = (lon >= -176.790) & (lon <= -157.820) & (lat >= -78.870) & (lat <= -77.520)
    mask_regions['WWED'] = (lon >= -65.130) & (lon <= -53.020) & (lat >= -75.950) & (lat <= -72.340)
    mask_regions['EWED'] = (lon >= -45.647) & (lon <= -32.253) & (lat >= -78.632) & (lat <= -76.899)
    
    return var.isel(deptht=bottom_cell).where(mask_regions[reg], drop=True)

def compute_VALSO_ocean_metrics(var_list_in, run_name_list, vector_T_list, vector_U_list, vector_V_list, cellarea_list, ocean_masks, smallest_domain=True):

    var_ds_list = []

    if smallest_domain:
        cc_prev = 10**23
        for n,cc in enumerate(cellarea_list):
            if cc.sum() < cc_prev:
                cc_prev = cc.sum().load()
                n_small = n
                
        cell_area = cellarea_list[n_small]
        print('USING CELL AREA FROM RUN NUMBER '+str(n)+' IN THE LIST')
    else:
        print('USING CELL AREA FROM EACH RUN SEPARATELY')
        
    for n,rrun in enumerate(run_name_list):
        print(run_name_list[n])
        
        var_to_plot = xr.Dataset()

        file_T_n = vector_T_list[n]
        file_U_n = vector_U_list[n]
        file_V_n = vector_V_list[n]

        #### format cell_area
        if not smallest_domain:
            cell_area = cellarea_list[n]

        if 'deptht' in cell_area.dims:
            cell_area_here = cell_area.isel(deptht=0)
        else:
            cell_area_here = cell_area

        mask_land = file_T_n['so'].isel(time=0,deptht=0).squeeze().drop('time')

        lon = file_T_n.nav_lon
        lat = file_T_n.nav_lat


        ### prepare integrated horizontal mass transport
        if ('ACC' in var_list_in) or ('wed_gyre' in var_list_in) or ('ross_gyre' in var_list_in):
            
            u_vertsum = (file_U_n['uocetr_eff'] * (mask_land.rename({'deptht':'depthu'}) > 0)).sum('depthu').load()

        #######
        if 'global_SST' in var_list_in:
            print('Computing global SST')
            if 'deptht' in cell_area.dims:
                cell_area_here = cell_area.isel(deptht=0)
            else:
                cell_area_here = cell_area
                
            var_to_plot['global_SST'] = uf.weighted_mean(file_T_n['thetao'].isel(deptht=0),['x','y'],cell_area_here)
        
        #######
        if 'ACC' in var_list_in:
            print('Computing ACC')
            var_to_plot['ACC'] = u_vertsum.sel(x=220,y=range(79,107)).sum('y')
            var_to_plot['ACC'] = var_to_plot['ACC']/10**6


        ########
        print('Computing gyres')
        if 'wed_gyre' in var_list_in:
            uocetr_Wed = u_vertsum.where(uf.in_range(cell_area_here.nav_lat,[-66.50,-60.40]) & uf.in_range(cell_area_here.nav_lon,[-31.25,37.50]), drop=True)

            var_to_plot['wed_gyre'] = uocetr_Wed.cumsum('y').max(['y','x'])/10**6
        
        ########
        if 'ross_gyre' in var_list_in:
            uocetr_Ross = u_vertsum.where(uf.in_range(cell_area_here.nav_lat,[-72.650,-61.600]) & ((cell_area_here.nav_lon <= -135.75) | (cell_area_here.nav_lon >= 360-168.500)), drop=True)
            
            var_to_plot['ross_gyre'] = uocetr_Ross.cumsum('y').max(['y','x'])/10**6

        ########
        print('Computing bottom properties')
        mask_cells = (file_T_n['so'] > 0).isel(time=0).drop('time').drop('time_centered').load()
        bottom_cell = mask_cells.sum('deptht') - 1

        if 'Tbot_AMU' in var_list_in:
            reg = 'AMU'
            var_to_plot['Tbot_'+reg] = bottom_prop_reg(file_T_n['thetao'],lon,lat,reg,bottom_cell).mean(['x','y'])
            var_to_plot['Sbot_'+reg] = bottom_prop_reg(file_T_n['so'],lon,lat,reg,bottom_cell).mean(['x','y'])

        if 'Tbot_WROSS' in var_list_in:
            reg = 'WROSS'
            var_to_plot['Tbot_'+reg] = bottom_prop_reg(file_T_n['thetao'],lon,lat,reg,bottom_cell).mean(['x','y'])
            var_to_plot['Sbot_'+reg] = bottom_prop_reg(file_T_n['so'],lon,lat,reg,bottom_cell).mean(['x','y'])

        if 'Tbot_EROSS' in var_list_in:
            reg = 'EROSS'
            var_to_plot['Tbot_'+reg] = bottom_prop_reg(file_T_n['thetao'],lon,lat,reg,bottom_cell).mean(['x','y'])
            var_to_plot['Sbot_'+reg] = bottom_prop_reg(file_T_n['so'],lon,lat,reg,bottom_cell).mean(['x','y'])

        if 'Tbot_WWED' in var_list_in:
            reg = 'WWED'
            var_to_plot['Tbot_'+reg] = bottom_prop_reg(file_T_n['thetao'],lon,lat,reg,bottom_cell).mean(['x','y'])
            var_to_plot['Sbot_'+reg] = bottom_prop_reg(file_T_n['so'],lon,lat,reg,bottom_cell).mean(['x','y'])

        if 'Tbot_EWED' in var_list_in:
            reg = 'EWED'
            var_to_plot['Tbot_'+reg] = bottom_prop_reg(file_T_n['thetao'],lon,lat,reg,bottom_cell).mean(['x','y'])
            var_to_plot['Sbot_'+reg] = bottom_prop_reg(file_T_n['so'],lon,lat,reg,bottom_cell).mean(['x','y'])

        ########
        print('Computing ocean heat content')
        
        rho0 = 1020.
        c_p = 4000.
        
        if 'OHC_tot' in var_list_in:
            
            var_to_plot['OHC_tot'] = rho0 * c_p * (file_T_n['thetao'] * cell_area_here * file_T_n['e3t']).sum(['x','y','deptht']) /10**22

        if 'OHC_700' in var_list_in:
            
            deptht_700 = (file_T_n.deptht <= 700).sum()

            var_to_plot['OHC_700'] = rho0 * c_p * (file_T_n['thetao'] * cell_area_here* file_T_n['e3t']).isel(deptht = np.arange(deptht_700)).sum(['x','y','deptht']) /10**22

        if 'OHC_2000' in var_list_in:
            
            deptht_2000 = (file_T_n.deptht <= 2000).sum()

            var_to_plot['OHC_2000'] = rho0 * c_p * (file_T_n['thetao'] * cell_area_here * file_T_n['e3t']).isel(deptht = np.arange(deptht_2000)).sum(['x','y','deptht']) /10**22

        if 'AMOC_26N' in var_list_in:
            
            Atl_26N_mask = (lat<=26.5) & (lat>=25.5) & np.isfinite(ocean_masks['atlantic'])
            v_atlsum_26N = file_V_n['vocetr_eff'].where(Atl_26N_mask,drop=True).sum('x')
            var_to_plot['AMOC_26N'] = (v_atlsum_26N.assign_coords({'depthv': -1*file_V_n['vocetr_eff'].depthv}).cumsum('depthv') / 10**6).max(['depthv','y'])


        if 'AMOC_30S' in var_list_in:

            Atl_30S_mask = (lat<=-29.5) & (lat>=-30.5) & np.isfinite(ocean_masks['atlantic'])
            v_atlsum_30S = file_V_n['vocetr_eff'].where(Atl_30S_mask,drop=True).sum('x')
            var_to_plot['AMOC_30S'] = (v_atlsum_30S.assign_coords({'depthv': -1*file_V_n['vocetr_eff'].depthv}).cumsum('depthv') / 10**6).min(['depthv','y'])

        
        var_ds_list.append(var_to_plot.assign_coords({'run': rrun}))

    var_ds_all = xr.concat(var_ds_list, dim='run')
    return var_ds_all

def compute_VALSO_seaice_metrics(var_seaice_list_in, run_name_list, vector_T_list, vector_seaice_list, cellarea_list, ocean_masks, smallest_domain=True):

    var_ds_list = []

    if smallest_domain:
        cc_prev = 10**23
        for n,cc in enumerate(cellarea_list):
            if cc.sum() < cc_prev:
                cc_prev = cc.sum().load()
                n_small = n
                
        cell_area = cellarea_list[n_small]
        print('USING CELL AREA FROM RUN NUMBER '+str(n)+' IN THE LIST')
    else:
        print('USING CELL AREA FROM EACH RUN SEPARATELY')
        
    for n,rrun in enumerate(run_name_list):
        print(run_name_list[n])
        
        var_to_plot = xr.Dataset()

        file_T_n = vector_T_list[n]
        file_seaice_n = vector_seaice_list[n]


        #### format cell_area
        if not smallest_domain:
            cell_area = cellarea_list[n]
    
        if 'deptht' in cell_area.dims:
            cell_area_here = cell_area.isel(deptht=0).drop('deptht').load()
        else:
            cell_area_here = cell_area.load()

        mask_land = file_T_n['so'].isel(time=0,deptht=0).squeeze().drop('time')

        ########
        
        lon = cell_area_here.nav_lon
        lat = cell_area_here.nav_lat
        
        mask_Arc = (lat >= 50) 
        mask_Ant = (lat <= -50) 
        
        ########

        if ('mar_sie_arc' in var_seaice_list_in) or ('sep_sie_arc' in var_seaice_list_in) or ('feb_sie_ant' in var_seaice_list_in) or ('sep_sie_ant' in var_seaice_list_in):
            print('Computing sea-ice extent')
            
            file_ice_15 = file_seaice_n['siconc'] > 0.15

            sie_Ant = (file_ice_15.where(mask_Ant) * cell_area_here).sum(['x','y']).rename({'time_counter':'time'}).load()
            sie_Arc = (file_ice_15.where(mask_Arc) * cell_area_here).sum(['x','y']).rename({'time_counter':'time'}).load()

        if 'mar_sie_arc' in var_seaice_list_in:
            vvar = sie_Arc.where(sie_Arc['time.month'] == 3, drop=True).squeeze()/10**12
            var_to_plot['mar_sie_arc'] = vvar.assign_coords({'time': vvar['time.year']})
            
        if 'sep_sie_arc' in var_seaice_list_in:
            vvar = sie_Arc.where(sie_Arc['time.month'] == 9, drop=True).squeeze()/10**12
            var_to_plot['sep_sie_arc'] = vvar.assign_coords({'time': vvar['time.year']})

        if 'feb_sie_ant' in var_seaice_list_in:
            vvar = sie_Ant.where(sie_Ant['time.month'] == 2, drop=True).squeeze()/10**12
            var_to_plot['feb_sie_ant'] = vvar.assign_coords({'time': vvar['time.year']})

        if 'sep_sie_ant' in var_seaice_list_in:
            vvar = sie_Ant.where(sie_Ant['time.month'] == 9, drop=True).squeeze()/10**12
            var_to_plot['sep_sie_ant'] = vvar.assign_coords({'time': vvar['time.year']})

        ########

        sia_Ant = (file_seaice_n['siconc'].where(mask_Ant) * cell_area_here).sum(['x','y']).rename({'time_counter':'time'}).load()
        sia_Arc = (file_seaice_n['siconc'].where(mask_Arc) * cell_area_here).sum(['x','y']).rename({'time_counter':'time'}).load()

        print('Computing sea-ice area')
        
        if 'mar_sia_arc' in var_seaice_list_in:
            vvar = sia_Arc.where(sia_Arc['time.month'] == 3, drop=True).squeeze()/10**12
            var_to_plot['mar_sia_arc'] = vvar.assign_coords({'time': vvar['time.year']})

        if 'sep_sia_arc' in var_seaice_list_in:
            vvar = sia_Arc.where(sia_Arc['time.month'] == 9, drop=True).squeeze()/10**12
            var_to_plot['sep_sia_arc'] = vvar.assign_coords({'time': vvar['time.year']})

        if 'feb_sia_ant' in var_seaice_list_in:
            vvar = sia_Ant.where(sia_Ant['time.month'] == 2, drop=True).squeeze()/10**12
            var_to_plot['feb_sia_ant'] = vvar.assign_coords({'time': vvar['time.year']})


        if 'sep_sia_ant' in var_seaice_list_in:
            vvar = sia_Ant.where(sia_Ant['time.month'] == 9, drop=True).squeeze()/10**12
            var_to_plot['sep_sia_ant'] = vvar.assign_coords({'time': vvar['time.year']})

        ########

        siv_Ant = (file_seaice_n['sivolu'].where(mask_Ant) * cell_area_here).sum(['x','y']).rename({'time_counter':'time'}).load()
        siv_Arc = (file_seaice_n['sivolu'].where(mask_Arc) * cell_area_here).sum(['x','y']).rename({'time_counter':'time'}).load()

        print('Computing sea-ice volume')
        
        if 'mar_siv_arc' in var_seaice_list_in:
            vvar = siv_Arc.where(siv_Arc['time.month'] == 3, drop=True).squeeze()/10**12
            var_to_plot['mar_siv_arc'] = vvar.assign_coords({'time': vvar['time.year']})

        if 'sep_siv_arc' in var_seaice_list_in:
            vvar = siv_Arc.where(siv_Arc['time.month'] == 9, drop=True).squeeze()/10**12
            var_to_plot['sep_siv_arc'] = vvar.assign_coords({'time': vvar['time.year']})

        if 'feb_siv_ant' in var_seaice_list_in:
            vvar = siv_Ant.where(siv_Ant['time.month'] == 2, drop=True).squeeze()/10**12
            var_to_plot['feb_siv_ant'] = vvar.assign_coords({'time': vvar['time.year']})
            
        if 'sep_siv_ant' in var_seaice_list_in:
            vvar = siv_Ant.where(siv_Ant['time.month'] == 9, drop=True).squeeze()/10**12
            var_to_plot['sep_siv_ant'] = vvar.assign_coords({'time': vvar['time.year']})
            
        ########
        
        var_ds_list.append(var_to_plot.assign_coords({'run': rrun}))

    var_ds_all = xr.concat(var_ds_list, dim='run')
    return var_ds_all

def VALSO_plot_ocean_comparison(var_list_in, var_ds_all, run_name_list, color_list):
    
    var_obs_mean, var_obs_std = var_obs_mean_std(var_list_in)

    var_to_plot = var_ds_all.load()
    
    f = plt.figure()
    f.set_size_inches(8.25*2, 8.25*2.2)
    
    ax={}
    
    leg_hdl = []
    
    i = 0
    
    for vv in var_list_in:
        
        ax[i] = f.add_subplot(5,6,i+1)
    
        if vv in list(var_to_plot.keys()):
            for n,nname in enumerate(run_name_list):
                ax[i].plot(var_to_plot[vv].sel(run=nname), color=color_list[n])
        
        ax[i].axhline(y=var_obs_mean.sel(var=vv), color='black', linewidth=2)
        ax[i].fill_between(x=np.arange(0,100),y1=var_obs_mean.sel(var=vv)-var_obs_std.sel(var=vv), y2=var_obs_mean.sel(var=vv)+var_obs_std.sel(var=vv), color='grey',alpha=0.2)
    
        if vv[0:3] == 'OHC':
            ax[i].set_title(vv+' x 10$^{22}$ J')
        else:
            ax[i].set_title(vv)
    
        i = i+1
    #f.legend()
    #f.subplots_adjust(bottom=0.05, wspace=0.1)
    
    f.tight_layout()
    sns.despine()

    return f

def VALSO_plot_seaice_comparison(var_list_in, var_ds_all, run_name_list, color_list):
    
    var_obs_mean, var_obs_std = var_obs_mean_std(var_list_in)

    #return var_ds_all
    var_to_plot = var_ds_all.load()
    
    f = plt.figure()
    f.set_size_inches(8.25*2, 8.25*2.2)
    
    ax={}
    
    leg_hdl = []
    
    i = 0
    
    for vv in var_list_in:
        
        ax[i] = f.add_subplot(5,6,i+1)
    
        if vv in list(var_to_plot.keys()):
            for n,nname in enumerate(run_name_list):
                ax[i].plot(var_to_plot[vv].sel(run=nname), color=color_list[n])
        
        ax[i].axhline(y=var_obs_mean.sel(var=vv), color='black', linewidth=2)
        ax[i].fill_between(x=np.arange(0,100),y1=var_obs_mean.sel(var=vv)-var_obs_std.sel(var=vv), y2=var_obs_mean.sel(var=vv)+var_obs_std.sel(var=vv), color='grey',alpha=0.2)
    
        if vv[0:3] == 'OHC':
            ax[i].set_title(vv+' x 10$^22$ J')
        else:
            ax[i].set_title(vv)
    
        i = i+1
    #f.legend()
    #f.subplots_adjust(bottom=0.05, wspace=0.1)
    
    f.tight_layout()
    sns.despine()

    return f, var_ds_all

