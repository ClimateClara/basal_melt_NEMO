# for the folder BLKDIR
path1=/ccc/work/cont003/igcmg/igcmg/IGCM/OCE/NEMO/FORCINGS/MRI-JRA55-do-1-4-0

ORCA_choice=eORCA1.4.2

for vv in {u_10,v_10,t_10,q_10,ncar_rad,ncar_precip,slp}
do
cp $path1/"$vv".15JUNE2009_fill.nc "$vv"_CORE2_fill.nc
done

cp $path1/weights_coreII_2_${ORCA_choice}_bilinear.nc weights_bilinear.nc
cp $path1/weights_coreII_2_${ORCA_choice}_bicubic.nc weights_bicubic.nc

# for the folder INITDIR
path2=/ccc/work/cont003/igcmg/igcmg/IGCM/OCE/NEMO/${ORCA_choice}/OPA

for vv in {maskMFO,resto,subbasins,sali_ref_clim_monthly}
do
cp $path2/${ORCA_choice}_${vv}.nc ${vv}.nc
done

cp $path2/${ORCA_choice}_sss_absolute_PHC2_salx_2004_08_03_clim.nc sss_absolute_salinity.nc
cp $path2/Lucazeau_ghflux.nc .
cp $path2/merged_ESACCI_BIOMER4V1R1_CHL_REG05.nc .
cp $path2/weights_WOA13d1_2_${ORCA_choice}_bilinear.nc weights_WOA13d1_2_bilinear.nc
cp $path2/weights_3D_WOA13d1_2_${ORCA_choice}_bilinear.nc weights_3D_WOA13d1_2_bilinear.nc
cp $path2/weights_reg05_2_${ORCA_choice}_bilinear.nc weights_reg05_2_bilinear.nc
cp $path2/weights_Lucazeau1_2_${ORCA_choice}_bilinear.nc weights_Lucazeau1_2_bilinear.nc


