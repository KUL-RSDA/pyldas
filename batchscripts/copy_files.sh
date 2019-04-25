#!/bin/bash

vscroot=/scratch/leuven/$(echo "$USER" |cut -c 4-6)/$USER/output/TEST_RUNS
macroot=/data_sets/LDAS_runs

exp=US_M36_SMOS40_DA_cal_scaled
dom=SMAP_EASEv2_M36_US
yr=2010
mo=01

obsparamf=$exp.ldas_obsparam.${yr}${mo}01_0000z.txt
catparamf=$exp.ldas_catparam.${yr}${mo}01_0000z.bin
rtmparamf=$exp.ldas_mwRTMparam.${yr}${mo}01_0000z.bin
drvparamf=$exp.ldas_driver_inputs.${yr}${mo}01_0000z.nml

mkdir $exp/output/$dom/ana/ens_avg
mkdir $exp/output/$dom/cat/ens_avg
mkdir $exp/output/$dom/rc_out/Y$yr/M01/

# tilegrids and tilecoords file
scp vsc:$vscroot/$exp/output/$dom/rc_out/* $macroot/$exp/output/$dom/rc_out/

# parameter files
scp vsc:$vscroot/$exp/output/$dom/rc_out/Y$yr/M01/$obsparamf $macroot/$exp/output/$dom/rc_out/Y$yr/M01
scp vsc:$vscroot/$exp/output/$dom/rc_out/Y$yr/M01/$catparamf $macroot/$exp/output/$dom/rc_out/Y$yr/M01
scp vsc:$vscroot/$exp/output/$dom/rc_out/Y$yr/M01/$rtmparamf $macroot/$exp/output/$dom/rc_out/Y$yr/M01
scp vsc:$vscroot/$exp/output/$dom/rc_out/Y$yr/M01/$drvparamf $macroot/$exp/output/$dom/rc_out/Y$yr/M01

# catchment model output
scp vsc:$vscroot/$exp/output/$dom/cat/ens_avg/xhourly_timeseries.nc $macroot/$exp/output/$dom/cat/ens_avg
scp vsc:$vscroot/$exp/output/$dom/cat/ens_avg/xhourly_images.nc $macroot/$exp/output/$dom/cat/ens_avg

# analysis output
scp vsc:$vscroot/$exp/output/$dom/ana/ens_avg/ObsFcstAna_timeseries.nc $macroot/$exp/output/$dom/ana/ens_avg
scp vsc:$vscroot/$exp/output/$dom/ana/ens_avg/ObsFcstAna_images.nc $macroot/$exp/output/$dom/ana/ens_avg
