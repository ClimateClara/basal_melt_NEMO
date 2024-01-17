#!/bin/bash

cas_path=/data/cdelaver
clara_path=/data/cburgard/CASIMIR_SIMU/interim/CDO_PROCESSED

##### PREPARE PATHS
path_closed=$cas_path/n42tm21
path_open=$cas_path/n42openc
path_closed2=$clara_path/n42tm21
path_open2=$clara_path/n42openc

### GLOBAL THETA
for yy in {0..9}
do

yy1=$(expr $yy + 1)
cdo fldmean -selvar,thetao $path_open/n42openc_00"$yy"10101_00"$yy1"01231_1Y_grid_T.nc $path_open2/n42openc_00"$yy"10101_00"$yy1"01231_1Y_thetao_fldmean.nc
cdo fldmean -selvar,thetao $path_closed/n42tm21_00"$yy"10101_00"$yy1"01231_1Y_grid_T.nc $path_closed2/n42tm21_00"$yy"10101_00"$yy1"01231_1Y_thetao_fldmean.nc

done

cdo mergetime $path_open2/n42openc_*_1Y_thetao_fldmean.nc $path_open2/n42openc_01-100_1Y_thetao_fldmean.nc
cdo mergetime $path_closed2/n42tm21_*_1Y_thetao_fldmean.nc $path_closed2/n42tm21_01-100_1Y_thetao_fldmean.nc

### 

### SEA ICE EXTENT

for yy in {0..9}
do

yy1=$(expr $yy + 1)
cdo sellonlatbox,0,360,-90,-50 $path_open/n42openc_00"$yy"10101_00"$yy1"01231_1M_icemod.nc $path_open2/n42openc_00"$yy"10101_00"$yy1"01231_1M_icemod_Ant.nc
cdo sellonlatbox,0,360,50,90 $path_open/n42openc_00"$yy"10101_00"$yy1"01231_1M_icemod.nc $path_open2/n42openc_00"$yy"10101_00"$yy1"01231_1M_icemod_Arc.nc
cdo sellonlatbox,0,360,-90,-50 $path_closed/n42tm21_00"$yy"10101_00"$yy1"01231_1M_icemod.nc $path_closed2/n42tm21_00"$yy"10101_00"$yy1"01231_1M_icemod_Ant.nc
cdo sellonlatbox,0,360,50,90 $path_closed/n42tm21_00"$yy"10101_00"$yy1"01231_1M_icemod.nc $path_closed2/n42tm21_00"$yy"10101_00"$yy1"01231_1M_icemod_Arc.nc

done

cdo mergetime $path_open2/n42openc_*_1M_icemod_Ant.nc $path_open2/n42openc_01-100_1M_icemod_Ant.nc
cdo mergetime $path_open2/n42openc_*_1M_icemod_Arc.nc $path_open2/n42openc_01-100_1M_icemod_Arc.nc
cdo mergetime $path_closed2/n42tm21_*_1M_icemod_Ant.nc $path_closed2/n42tm21_01-100_1M_icemod_Ant.nc
cdo mergetime $path_closed2/n42tm21_*_1M_icemod_Arc.nc $path_closed2/n42tm21_01-100_1M_icemod_Arc.nc

cdo gtc,0.15 -selvar,siconc $path_open2/n42openc_01-100_1M_icemod_Ant.nc $path_open2/n42openc_01-100_1M_icemod_Ant_SIEmask.nc
cdo gtc,0.15 -selvar,siconc $path_open2/n42openc_01-100_1M_icemod_Arc.nc $path_open2/n42openc_01-100_1M_icemod_Arc_SIEmask.nc
cdo gtc,0.15 -selvar,siconc $path_closed2/n42tm21_01-100_1M_icemod_Ant.nc $path_closed2/n42tm21_01-100_1M_icemod_Ant_SIEmask.nc
cdo gtc,0.15 -selvar,siconc $path_closed2/n42tm21_01-100_1M_icemod_Arc.nc $path_closed2/n42tm21_01-100_1M_icemod_Arc_SIEmask.nc

cdo fldsum -ifthen $path_closed2/n42tm21_01-100_1M_icemod_Ant_SIEmask.nc -selvar,cell_area $path_closed2/n42tm21_01-100_1M_icemod_Ant.nc $path_closed2/n42tm21_01-100_1M_icemod_Ant_SIE.nc
cdo fldsum -ifthen $path_closed2/n42tm21_01-100_1M_icemod_Arc_SIEmask.nc -selvar,cell_area $path_closed2/n42tm21_01-100_1M_icemod_Arc.nc $path_closed2/n42tm21_01-100_1M_icemod_Arc_SIE.nc

cdo fldsum -ifthen $path_open2/n42openc_01-100_1M_icemod_Ant_SIEmask.nc -selvar,cell_area $path_open2/n42openc_01-100_1M_icemod_Ant.nc $path_open2/n42openc_01-100_1M_icemod_Ant_SIE.nc
cdo fldsum -ifthen $path_open2/n42openc_01-100_1M_icemod_Arc_SIEmask.nc -selvar,cell_area $path_open2/n42openc_01-100_1M_icemod_Arc.nc $path_open2/n42openc_01-100_1M_icemod_Arc_SIE.nc


### CUT OUT REGIONS
for yy in {0..9}
do

yy1=$(expr $yy + 1)

# Amundsen 
cdo sellonlatbox,-109.64,-102.23,-75.80,-71.66 -selvar,thetao,so $path_open/n42openc_00"$yy"10101_00"$yy1"01231_1Y_grid_T.nc $path_open2/n42openc_00"$yy"10101_00"$yy1"01231_1Y_T_AMU.nc
cdo sellonlatbox,-109.64,-102.23,-75.80,-71.66 -selvar,thetao,so $path_closed/n42tm21_00"$yy"10101_00"$yy1"01231_1Y_grid_T.nc $path_closed2/n42tm21_00"$yy"10101_00"$yy1"01231_1Y_T_AMU.nc

# WRoss
cdo sellonlatbox,157.100,173.333,-78.130,-74.040 -selvar,thetao,so $path_open/n42openc_00"$yy"10101_00"$yy1"01231_1Y_grid_T.nc $path_open2/n42openc_00"$yy"10101_00"$yy1"01231_1Y_T_WROSS.nc
cdo sellonlatbox,157.100,173.333,-78.130,-74.040 -selvar,thetao,so $path_closed/n42tm21_00"$yy"10101_00"$yy1"01231_1Y_grid_T.nc $path_closed2/n42tm21_00"$yy"10101_00"$yy1"01231_1Y_T_WROSS.nc

# ERoss
cdo sellonlatbox,-176.790,-157.820,-78.870,-77.520 -selvar,thetao,so $path_open/n42openc_00"$yy"10101_00"$yy1"01231_1Y_grid_T.nc $path_open2/n42openc_00"$yy"10101_00"$yy1"01231_1Y_T_EROSS.nc
cdo sellonlatbox,-176.790,-157.820,-78.870,-77.520 -selvar,thetao,so $path_closed/n42tm21_00"$yy"10101_00"$yy1"01231_1Y_grid_T.nc $path_closed2/n42tm21_00"$yy"10101_00"$yy1"01231_1Y_T_EROSS.nc

# Weddell
cdo sellonlatbox,-65.130,-53.020,-75.950,-72.340  -selvar,thetao,so $path_open/n42openc_00"$yy"10101_00"$yy1"01231_1Y_grid_T.nc $path_open2/n42openc_00"$yy"10101_00"$yy1"01231_1Y_T_WED.nc
cdo sellonlatbox,-65.130,-53.020,-75.950,-72.340  -selvar,thetao,so $path_closed/n42tm21_00"$yy"10101_00"$yy1"01231_1Y_grid_T.nc $path_closed2/n42tm21_00"$yy"10101_00"$yy1"01231_1Y_T_WED.nc

# EWeddell
cdo sellonlatbox,-45.647,-32.253,-78.632,-76.899  -selvar,thetao,so $path_open/n42openc_00"$yy"10101_00"$yy1"01231_1Y_grid_T.nc $path_open2/n42openc_00"$yy"10101_00"$yy1"01231_1Y_T_EWED.nc
cdo sellonlatbox,-45.647,-32.253,-78.632,-76.899  -selvar,thetao,so $path_closed/n42tm21_00"$yy"10101_00"$yy1"01231_1Y_grid_T.nc $path_closed2/n42tm21_00"$yy"10101_00"$yy1"01231_1Y_T_EWED.nc

done

for reg in {AMU,WROSS,EROSS,WED,EWED}
do
cdo mergetime $path_open2/n42openc_*_1Y_T_"$reg".nc $path_open2/n42openc_01-100_1Y_T_"$reg".nc
cdo mergetime $path_closed2/n42tm21_*_1Y_T_"$reg".nc $path_closed2/n42tm21_01-100_1Y_T_"$reg".nc

done

### ACC

for yy in {0..9}
do

yy1=$(expr $yy + 1)

cdo vertsum -selvar,uocetr_eff $path_open/n42openc_00"$yy"10101_00"$yy1"01231_1Y_grid_U.nc $path_open2/n42openc_00"$yy"10101_00"$yy1"01231_1Y_uo_vertsum.nc
cdo vertsum -selvar,uocetr_eff $path_closed/n42tm21_00"$yy"10101_00"$yy1"01231_1Y_grid_U.nc $path_closed2/n42tm21_00"$yy"10101_00"$yy1"01231_1Y_uo_vertsum.nc

done 

cdo mergetime $path_open2/n42openc_*_1Y_uo_vertsum.nc $path_open2/n42openc_0-100_1Y_uo_vertsum.nc
cdo mergetime $path_closed2/n42tm21_*_1Y_uo_vertsum.nc $path_closed2/n42tm21_0-100_1Y_uo_vertsum.nc

### ICE SHELF MELT
for yy in {0..9}
do

yy1=$(expr $yy + 1)
cdo -selvar,iceshelf_cav,iceshelf $path_open/n42openc_00"$yy"10101_00"$yy1"01231_1Y_grid_T.nc $path_open2/n42openc_00"$yy"10101_00"$yy1"01231_1Y_grid_isfvars.nc

done

cdo mergetime $path_open2/n42openc_*_1Y_grid_isfvars.nc $path_open2/n42openc_0-100_1Y_grid_isfvars.nc

### WEDDELL AND ROSS GYRE

#Weddell
cdo -31.250 37.500 -66.500 -60.400