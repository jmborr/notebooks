#!/bin/bash

module load lammps/5-Sep-2014  #lammps version used for these project
PROJD=/SNSlocal/projects/jbq/pStyrene_Toluene

##############################
##  Semi-equilibration of r0.5
##############################
ratio="0.5"
# generation of lammps file
mkdir -p $PROJD/8styrene32/r${ratio}/msi2lmp
cd $PROJD/8styrene32/r${ratio}/msi2lmp/
ln -s ../8styrene32.car
ln -s ../8styrene32.mdf
msi2lmp 8styrene32 -class I -frc cvff.frc  # we need force-field file cvff.frc in this directory
# minimization
mkdir -p $PROJD/8styrene32/r${ratio}/minimize
cd $PROJD/8styrene32/r${ratio}/minimize
ln -s ../msi2lmp/8styrene32.data
mpirun -np 12 lmp_camm2 -in minimize.in
# compression
mkdir -p $PROJD/8styrene32/r${ratio}/compress
cd $PROJD/8styrene32/r${ratio}/compress
ln -s ../minimize/minimize.rst compress.0.rst
mpirun -np 12 lmp_camm2 -in compress.1.in
mpirun -np 12 lmp_camm2 -in compress.2.in
mpirun -np 12 lmp_camm2 -in compress.3.in
mpirun -np 12 lmp_camm2 -in compress.4.in
mpirun -np 12 lmp_camm2 -in compress.5.in
# heatup
mkdir -p $PROJD/8styrene32/r${ratio}/heatup
cd $PROJD/8styrene32/r${ratio}/heatup
ln -s ../compress/compress.5.rst compress.rst
time mpirun -np 12 lmp_camm2 -in heatup.in

#for T in T150 T170 T190 T210 T220 T230 T250 T270 T290 T310 T330 T350 T360 T370 T380 T390 T400 T405 T410 T415 T420 T430 T440 T450 T460 T470;do
for T in T270 T280 T290 T310 T330 T350 T370 T390 T400 T410 T420 T430 T440 T450;do
  cd $PROJD/8styrene32/r0.5/$T
  echo -e "\n\n PROCESSING T=$T\n"; sleep 2s
  #python $PROJD/python/merge_trajectories.py ../8styrene32.car ../8styrene32.mdf 3 equil.dcd --name equil.X_$T.dcd
  #python $PROJD/python/8styrene_diffusion_polymerCM.py ../8styrene32.pdb equil.dcd styCMH_r.xmgr # MSD(t) for CM of polymers
  #ptraj ../8styrene32.pdb < $PROJD/cpptraj/8styrene32_atomicfluct_styH.in #use module amber12
  #ptraj ../8styrene32.pdb < $PROJD/cpptraj/average_structure.in
  ptraj ../8styrene32.pdb < $PROJD/cpptraj/rms2avg_styH.in  #RMS of styrene hydrogens to average structure of the whole system
  #ptraj ../8styrene32.pdb < $PROJD/cpptraj/8styrene32_diffusion_styH.in
done
popwindow.py 0s 'finished pStyrene/bash_commands'


# Coalesce the 30ns trajectories
PROJD=/SNSlocal/projects/jbq/pStyrene_Toluene
for T in T380 T390 T400 T405 T410 T420 T430 T440 T450 T460 T470;do
  cd $PROJD/8styrene32/r0/$T/
  vmd -dispdev text -eofexit -e $PROJD/vmd/coalesce_trajectories.tcl
done


# Remove global rotations and translations
PROJD=/SNSlocal/projects/jbq/pStyrene_Toluene
for T in T400 T420 T430;do
 cd $PROJD/8styrene32/r1/$T/
 #ln -s equil.1.dcd equil.dcd
 #vmd -dispdev win -eofexit -e $PROJD/vmd/rms2first.tcl
 echo -n "$T "; dumpdcd equil_rms2first.dcd |grep 1000
done

# Find centroid and average structure by doing RMS over styH hydrogens using 500 conformations
PROJD=/SNSlocal/projects/jbq/pStyrene_Toluene
#for T in T150 T170 T190 T210 T220 T230 T250 T270 T290 T310 T330 T360 T370 T380 T390 T400 T405 T410 T420 T430 T440 T450 T460 T470;do
  #cd $PROJD/8styrene32/r0/$T/
#for T in T270 T280 T290 T310 T330 T350 T370 T380 T400 T410 T420 T430 T440 T450;do   #for r0.5
for T in T400 T410 T420 T430 T440 T450;do   #for r0.5
  cd $PROJD/8styrene32/r0.5/$T/
  python $PROJD/python/cluster_trajectory.py ./8styrene32.pdb ./equil_rms2first.dcd 500 '@1-4133&@H*' centroid_styH.pdb average_styH.pdb
do


         #############                                  
         ##   r0.0  ##
         #############
r='r0'; WD=$PROJD/8styrene32/r0
for temp in `seq 90 10 500`;do
  T=T$temp; #echo -e "\n##############\n T = $T\n##########"
  cd $WD/$T
  #ln -s ../8styrene32.pdb pdb
  #ln -s ../8styrene32.pdb
  if [[ "T380 T390 T400 T410 T420 T430 T440 T450 T460 T470 T480 T490 T500" == *${T}* ]]
  then
    nframes=100000
    #vmd -dispdev text -eofexit -e $PROJD/vmd/coalesce_trajectories.tcl; sleep 9s
    #cpptraj -p pdb -i $PROJD/amber/radial_8styrene32_r0_105001.cpptraj
  else
    nframes=10000
    #ln -s equil.1.dcd equil.dcd
    #cpptraj -p pdb -i $PROJD/amber/radial_8styrene32_r0_15001.cpptraj
  fi
  #echo -n "$T "; dumpdcd equil.dcd |head -1
  #vmd -dispdev win -eofexit -e $PROJD/vmd/rms2first.tcl; sleep 9s
  #echo -n "$T "; dumpdcd equil_rms2first.dcd | head -1
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-256)&(@H*)' atomicfluct_4.56ns_styH.dat --seriesfile atomicfluct_4.56ns_styH_series.dat --msd yes #atomic fluctuations of all hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styH.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-256)&(@H1,@H2,@H3,@H4,@H5)' atomicfluct_4.56ns_styd3.dat --seriesfile atomicfluct_4.56ns_styd3_series.dat --msd yes #atomic fluctuations of all styrene ring hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styd3.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-256)&(@H62,@H71,@H72,@H73)' atomicfluct_4.56ns_styd5.dat --seriesfile atomicfluct_4.56ns_styd5_series.dat --msd yes #atomic fluctuations of all styrene backbone hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styd5.dat
done

  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 1220 100 '(:1-256)&(@H*)' atomicfluct_1.22ns_styH.dat --seriesfile atomicfluct_1.22ns_styH_series.dat --msd yes &  #atomic fluctuations of all hydrogens in the time-scale of 1.22ns
  #echo -n "$T "; cat atomicfluct_1.22ns_styH.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 7640 100 '(:1-256)&(@H*)' atomicfluct_7.64ns_styH.dat --seriesfile atomicfluct_7.64ns_styH_series.dat --msd yes &  #atomic fluctuations of all hydrogens in the time-scale of 7.64ns
  #echo -n "$T "; cat atomicfluct_7.64ns_styH.dat
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 1073 200 '(:1-256)&(@H*)' atomicfluct_BASIS_sty.dat --seriesfile  atomicfluct_BASIS_sty_series.dat --msd yes & #atomic fluctuations of all hydrogens in the time-scale of BASIS (1.073ns)
  #echo -n "$T "; cat atomicfluct_BASIS_sty.dat
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 1073 200 '(:1-256)&(@H61,H62,H71,H72,H73)' atomicfluct_BASIS_bk.dat --seriesfile  atomicfluct_BASIS_bk_series.dat --msd yes & #atomic fluctuations of bakbone hydrogens in the time-scale of BASIS (1.073ns)
  #echo -n "$T "; cat atomicfluct_BASIS_bk.dat
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 1073 200 '(:1-256)&(@H1,H2,H3,H4,H5)' atomicfluct_BASIS_sc.dat --seriesfile  atomicfluct_BASIS_sc_series.dat --msd yes & #atomic fluctuations of side-chain hydrogens in the time-scale of BASIS (1.073ns)
  #echo -n "$T "; cat atomicfluct_BASIS_sc.dat
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4562 200 '(:1-256)&(@H*)' atomicfluct_HFBS_sty.dat --seriesfile  atomicfluct_HFBS_sty_series.dat --msd yes & #atomic fluctuations of all hydrogens in the time-scale of HFBS (4.562ns)
  #echo -n "$T "; cat atomicfluct_HFBS_sty.dat
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 10 200 '(:1-256)&(@H*)' atomicfluct_10ps_sty.dat --seriesfile  atomicfluct_10ps_sty_series.dat --msd yes &
  #sleep 9s
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 100 200 '(:1-256)&(@H*)' atomicfluct_100ps_sty.dat --seriesfile  atomicfluct_100ps_sty_series.dat --msd yes &
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 5000 200 '(:1-256)&(@H*)' atomicfluct_5ns_sty.dat --seriesfile  atomicfluct_5ns_sty_series.dat --msd yes &
  #python $PROJD/python/diffusion_t0average.py pdb equil_rms2first.dcd $nframes 1.0 5000 100 '(:1-256)&(@H*)' diffusion_sty.dat &
  #python $PROJD/python/residueCMtrajectory.py pdb equil_rms2first.dcd 256 residueCM.pdb equil_rms2first_residueCM.crd & #trajectory of the center of mass of the styrene monomers
  #python $PROJD/python/diffusion_t0average.py residueCM.pdb equil_rms2first_residueCM.crd $nframes 1.0 5000 100 ':1-256' diffusion_styCM.dat &  
  #time python $PROJD/python/atomic_fluct_t0average.py residueCM.pdb equil_rms2first_residueCM.crd $nframes 1073 60 ':1-256' atomicfluct_BASIS_styCM.dat --seriesfile atomicfluct_BASIS_styCM_series.dat --msd yes & #atomic fluctuations of the center of mass of the styrene monomers
  #python $PROJD/python/CMtrajectory.py pdb equil_rms2first.dcd 32 8 CM.pdb equil_rms2first_CM.crd
  #time python $PROJD/python/atomic_fluct_t0average.py CM.pdb equil_rms2first_CM.crd $nframes 1073 99 ':1-8' atomicfluct_BASIS_CM.dat --seriesfile atomicfluct_BASIS_CM_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_BASIS_CM.dat
  #python $PROJD/python/diffusion_t0average.py CM.pdb equil_rms2first_CM.crd $nframes 1.0 5000 100 ':1-8' diffusion_CM.dat &
  #python $PROJD/python/thermo_averages.py  equil.log
  #for i in `seq 1 8`;do for f in backbone sidechain;do python $PROJD/python/generate_mwcovar.py $i $f;done;done
done

         #############                                  
         ##  r0.2   ##
         #############
#REDOING SIMULATIONS FOR 410. REMOVE LINE if [[ "T410" != *${T}* ]];then WHEN FINISHED
r='r0.2'; WD=$PROJD/8styrene32/$r
for temp in `seq 80 10 490`;do
  T=T$temp; #echo -e "\n##############\n T = $T\n##########"
if [[ "T410" != *${T}* ]];then
  cd $WD/$T
  #ln -s ../8styrene32.pdb
  #ln -s ../8styrene32.pdb pdb
  if [[ "T260 T270 T280 T290 T300 T310 T320 T330 T340 T350 T360 T370 T380 T390 T400 T410 T420 T430 T440 T450 T460 T470 T480 T490" == *${T}* ]]
  then
    nframes=100000
    #vmd -dispdev text -eofexit -e $PROJD/vmd/coalesce_trajectories.tcl
    #cpptraj -p pdb -i $PROJD/amber/radial_8styrene32_105001.cpptraj
    #python $PROJD/python/CMAnalysis.py 8styrene32.pdb equil.dcd 5001 105001 100 'not resnum 1:256 and name C' 6.11 contacts_tolC_tolC.h5
    #python $PROJD/python/CMAnalysis.py 8styrene32.pdb equil.dcd 5001 105001 100 'not resnum 1:256 and name C' 7.0 contacts_tolC_tolC_7.0.h5
  else
    nframes=10000
    #ln -s equil.1.dcd equil.dcd
    #cpptraj -p pdb -i $PROJD/amber/radial_8styrene32_15001.cpptraj
    #python $PROJD/python/CMAnalysis.py 8styrene32.pdb equil.dcd 5001 15001 10 'not resnum 1:256 and name C' 6.11 contacts_tolC_tolC.h5
    #python $PROJD/python/CMAnalysis.py 8styrene32.pdb equil.dcd 5001 15001 10 'not resnum 1:256 and name C' 7.0 contacts_tolC_tolC_7.0.h5
  fi
  #vmd -dispdev win -eofexit -e $PROJD/vmd/rms2first.tcl
  #echo -n "$T "; dumpdcd equil_rms2first.dcd | head -1
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-256)&(@H*)' atomicfluct_4.56ns_styH.dat --seriesfile atomicfluct_4.56ns_styH_series.dat --msd yes #atomic fluctuations of all polymer hydrogens in the time-scale of HFBS
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-256)&(@H1,@H2,@H3,@H4,@H5)' atomicfluct_4.56ns_styd3.dat --seriesfile atomicfluct_4.56ns_styd3_series.dat --msd yes #atomic fluctuations of all styrene ring hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styd3.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-256)&(@H62,@H71,@H72,@H73)' atomicfluct_4.56ns_styd5.dat --seriesfile atomicfluct_4.56ns_styd5_series.dat --msd yes #atomic fluctuations of all styrene backbone hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styd5.dat
  #echo -n "$T "; cat atomicfluct_4.56ns_styH.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(!(:1-256))&(@H=)' atomicfluct_4.56ns_tol.dat --seriesfile  atomicfluct_4.56ns_tol_series.dat --msd yes #atomic fluctuations of all toluene hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_tol.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(!(:1-256))&(@H1,@H2,@H3,@H4,@H5)' atomicfluct_4.56ns_told3.dat --seriesfile  atomicfluct_4.56ns_told3_series.dat --msd yes #atomic fluctuations of all toluene hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_tol.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(!(:1-256))&(@C)' atomicfluct_4.56ns_tolC.dat --seriesfile  atomicfluct_4.56ns_tolC_series.dat --msd yes #atomic fluctuations of central carobon of toluene in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_tolC.dat
  #mkdir atomicfluct_4.56ns_told3
  #for ires in `seq 257 307`;do
  #  atoms="':${ires}&(@H1,@H2,@H3,@H4,@H5)'"
  #  python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 $atoms atomicfluct_4.56ns_told3/atomicfluct_4.56ns_tol${ires}d3.dat --seriesfile  atomicfluct_4.56ns_told3/atomicfluct_4.56ns_tol${ires}d3_series.dat --msd yes #atomic fluctuations of all toluene hydrogens in the time-scale of HFBS
  #done
  #mkdir atomicfluct_4.56ns_styd3
  #for ires in `seq 1 256`;do
  #  atoms="':${ires}&(@H1,@H2,@H3,@H4,@H5)'"
  #  python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 $atoms atomicfluct_4.56ns_styd3/atomicfluct_4.56ns_sty${ires}d3.dat --seriesfile  atomicfluct_4.56ns_styd3/atomicfluct_4.56ns_sty${ires}d3_series.dat --msd yes #atomic fluctuations of all toluene hydrogens in the time-scale of HFBS
  #done
  #mkdir atomicfluct_4.56ns_tolC
  #for ires in `seq 257 307`;do
  #  atoms="':${ires}&(@C)'"
  #  python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 $atoms atomicfluct_4.56ns_tolC/atomicfluct_4.56ns_tol${ires}C.dat --seriesfile  atomicfluct_4.56ns_tolC/atomicfluct_4.56ns_tol${ires}C_series.dat --msd yes #atomic fluctuations of central carbon of tol in the time-scale of HFBS
  #done
  #watershell (:257)&(@C) sty_around_tol257.dat lower 6.0 upper 8.0 (:1-256)&(@C)
  #python $PROJD/python/cluster_contacts.py contacts_tolC_tolC.h5 257-307 contacts_tolC_tolC_cluster_sizes.dat
  #echo -n "$temp "; head -2 contacts_tolC_tolC_cluster_sizes.dat|tail --lines=1|tr '#' ' '
  #python $PROJD/python/cluster_contacts.py contacts_tolC_tolC_7.0.h5 257-307 contacts_tolC_tolC_7.0_cluster_sizes.dat
  #echo -n "$temp "; head -2 contacts_tolC_tolC_7.0_cluster_sizes.dat|tail --lines=1|tr '#' ' '
fi
done
popwindow.py 1s 'Finished r0.2 jobs'

         #############                                  
         ##  r0.3   ##
         #############
#REDOING SIMULATIONS FOR 370. REMOVE LINE if [[ "T370" != *${T}* ]];then WHEN FINISHED
r='r0.3'; WD=$PROJD/8styrene32/r0.3
for temp in `seq 80 10 490`;do
  T=T$temp; #echo -e "\n##############\n T = $T\n##########"
if [[ "T370" != *${T}* ]];then
  cd $WD/$T
  #ln -s ../8styrene32.pdb
  #ln -s ../8styrene32.pdb pdb
  if [[ "T260 T270 T280 T290 T300 T310 T320 T330 T340 T350 T360 T370 T380 T390 T400 T410 T420 T430 T440 T450 T460 T470 T480 T490" == *${T}* ]]
  then
    nframes=100000
    #vmd -dispdev text -eofexit -e $PROJD/vmd/coalesce_trajectories.tcl
    #cpptraj -p pdb -i $PROJD/amber/radial_8styrene32_105001.cpptraj
    #python $PROJD/python/CMAnalysis.py 8styrene32.pdb equil.dcd 5001 105001 100 'not resnum 1:256 and name C' 6.11 contacts_tolC_tolC.h5
    #python $PROJD/python/CMAnalysis.py 8styrene32.pdb equil.dcd 5001 105001 100 'not resnum 1:256 and name C' 7.0 contacts_tolC_tolC_7.0.h5
  else
    nframes=10000
    #ln -s equil.1.dcd equil.dcd
    #cpptraj -p pdb -i $PROJD/amber/radial_8styrene32_15001.cpptraj
    #python $PROJD/python/CMAnalysis.py 8styrene32.pdb equil.dcd 5001 15001 10 'not resnum 1:256 and name C' 6.11 contacts_tolC_tolC.h5
    #python $PROJD/python/CMAnalysis.py 8styrene32.pdb equil.dcd 5001 15001 10 'not resnum 1:256 and name C' 7.0 contacts_tolC_tolC_7.0.h5
  fi
  #vmd -dispdev win -eofexit -e $PROJD/vmd/rms2first.tcl; sleep 9s
  #echo -n "$T "; dumpdcd equil.dcd|head -1
  #echo -n "$T "; dumpdcd equil_rms2first.dcd | head -1
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-256)&(@H*)' atomicfluct_4.56ns_styH.dat --seriesfile atomicfluct_4.56ns_styH_series.dat --msd yes #atomic fluctuations of all hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styH.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-256)&(@H1,@H2,@H3,@H4,@H5)' atomicfluct_4.56ns_styd3.dat --seriesfile atomicfluct_4.56ns_styd3_series.dat --msd yes #atomic fluctuations of all styrene ring hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styd3.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-256)&(@H62,@H71,@H72,@H73)' atomicfluct_4.56ns_styd5.dat --seriesfile atomicfluct_4.56ns_styd5_series.dat --msd yes #atomic fluctuations of all styrene backbone hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styd5.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(!(:1-256))&(@H=)' atomicfluct_4.56ns_tol.dat --seriesfile  atomicfluct_4.56ns_tol_series.dat --msd yes #atomic fluctuations of all toluene hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_tol.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(!(:1-256))&(@C)' atomicfluct_4.56ns_tolC.dat --seriesfile  atomicfluct_4.56ns_tolC_series.dat --msd yes #atomic fluctuations of central carobon of toluene in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_tolC.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(!(:1-256))&(@H1,@H2,@H3,@H4,@H5)' atomicfluct_4.56ns_told3.dat --seriesfile  atomicfluct_4.56ns_told3_series.dat --msd yes #atomic fluctuations of all toluene hydrogens in the time-scale of HFBS
  #python $PROJD/python/cluster_contacts.py contacts_tolC_tolC.h5 257-333 contacts_tolC_tolC_cluster_sizes.dat 
  #echo -n "$temp "; head -2 contacts_tolC_tolC_cluster_sizes.dat|tail --lines=1|tr '#' ' '
  #python $PROJD/python/cluster_contacts.py contacts_tolC_tolC_7.0.h5 257-333 contacts_tolC_tolC_7.0_cluster_sizes.dat
  #echo -n "$temp "; head -2 contacts_tolC_tolC_7.0_cluster_sizes.dat|tail --lines=1|tr '#' ' '
fi
done
popwindow.py 1s 'Finished r0.3 jobs'

         #############                                  
         ##  r0.4   ##
         #############
r='r0.4'; WD=$PROJD/8styrene32/r0.4
for temp in `seq 80 10 490`;do
  T=T$temp;
  #echo -e "\n##############\n T = $T\n##########"
  cd $WD/$T
  #ln -s ../8styrene32.pdb
  #ln -s ../8styrene32.pdb pdb
  if [[ "T260 T270 T280 T290 T300 T310 T320 T330 T340 T350 T360 T370 T380 T390 T400 T410 T420 T430 T440 T450 T460 T470 T480 T490" == *${T}* ]]
  then
    nframes=100000
    #vmd -dispdev text -eofexit -e $PROJD/vmd/coalesce_trajectories.tcl; sleep 9s
    #cpptraj -p pdb -i $PROJD/amber/radial_8styrene32_105001.cpptraj
    #python $PROJD/python/CMAnalysis.py 8styrene32.pdb equil.dcd 5001 105001 100 'not resnum 1:256 and name C' 6.11 contacts_tolC_tolC.h5
    #python $PROJD/python/CMAnalysis.py 8styrene32.pdb equil.dcd 5001 105001 100 'not resnum 1:256 and name C' 7.0 contacts_tolC_tolC_7.0.h5
  else
    nframes=10000
    #ln -s equil.1.dcd equil.dcd
    #cpptraj -p pdb -i $PROJD/amber/radial_8styrene32_15001.cpptraj
    #python $PROJD/python/CMAnalysis.py 8styrene32.pdb equil.dcd 5001 15001 10 'not resnum 1:256 and name C' 6.11 contacts_tolC_tolC.h5
    #python $PROJD/python/CMAnalysis.py 8styrene32.pdb equil.dcd 5001 15001 10 'not resnum 1:256 and name C' 7.0 contacts_tolC_tolC_7.0.h5
  fi
  #echo -n "$T "; dumpdcd equil.dcd|head -1
  #python $PROJD/python/CMsystemTrajectory.py 8styrene32.pdb equil.dcd cm.pdb cm.crd & #get trajectory of the CM of the whole system
  #echo -n "$T ";   head -2 cm.crd;   tail --lines=2 cm.crd
  #vmd -dispdev text -e $PROJD/vmd/crop_equil.tcl # Remove first 5ns, creating equil_cropped.dcd
  #echo -n "$T "; dumpdcd equil_cropped.dcd | head -1
  #python $PROJD/python/diffusion_t0average.py pdb equil_cropped.dcd $nframes 1.0 8000 100 '(!:1-256)&(@H1,@H2,@H3,@H4,@H5)' diffusion_styd3.dat --rms2t0 no & # MSD(t) of the hydrogens in toluene rings
  #echo -n "$T "; python $PROJD/python/benedettoMSD.py diffusion_styd3.dat BASIS #BASIS MSD according to Benedetto12 reference
  echo -n "$T "; python $PROJD/python/benedettoMSD.py diffusion_styd3.dat HFBS #HFBS MSD according to Benedetto12 reference
  #python $PROJD/python/CMtraj.py pdb equil_cropped.dcd '(!:1-256)&(@H1,@H2,@H3,@H4,@H5)' styd3CoM.pdb styd3CoM.dcd --aggregate byres & #trajectory and topology of Center of Masses of the hydrogens of toluene rings
  
  #vmd -dispdev win -eofexit -e $PROJD/vmd/rms2first.tcl; sleep 9s
  #echo -n "$T "; dumpdcd equil_rms2first.dcd | head -1
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-256)&(@H*)' atomicfluct_4.56ns_styH.dat --seriesfile atomicfluct_4.56ns_styH_series.dat --msd yes #atomic fluctuations of all hydrogens in the time-scale of HFBS
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-256)&(@H1,@H2,@H3,@H4,@H5)' atomicfluct_4.56ns_styd3.dat --seriesfile atomicfluct_4.56ns_styd3_series.dat --msd yes #atomic fluctuations of all styrene ring hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styd3.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-256)&(@H62,@H71,@H72,@H73)' atomicfluct_4.56ns_styd5.dat --seriesfile atomicfluct_4.56ns_styd5_series.dat --msd yes #atomic fluctuations of all styrene backbone hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styd5.dat
  #echo -n "$T "; cat atomicfluct_4.56ns_styH.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(!(:1-256))&(@H=)' atomicfluct_4.56ns_tol.dat --seriesfile  atomicfluct_4.56ns_tol_series.dat --msd yes #atomic fluctuations of all toluene hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_tol.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(!(:1-256))&(@C)' atomicfluct_4.56ns_tolC.dat --seriesfile  atomicfluct_4.56ns_tolC_series.dat --msd yes #atomic fluctuations of central carobon of toluene in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_tolC.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(!(:1-256))&(@H1,@H2,@H3,@H4,@H5)' atomicfluct_4.56ns_told3.dat --seriesfile  atomicfluct_4.56ns_told3_series.dat --msd yes #atomic fluctuations of all toluene hydrogens in the time-scale of HFBS
  #python $PROJD/python/cluster_contacts.py contacts_tolC_tolC.h5 257-358 contacts_tolC_tolC_cluster_sizes.dat
  #echo -n "$temp "; head -2 contacts_tolC_tolC_cluster_sizes.dat|tail --lines=1|tr '#' ' '
  #python $PROJD/python/cluster_contacts.py contacts_tolC_tolC_7.0.h5 257-358 contacts_tolC_tolC_7.0_cluster_sizes.dat
  #echo -n "$temp "; head -2 contacts_tolC_tolC_7.0_cluster_sizes.dat|tail --lines=1|tr '#' ' '
done
popwindow.py 1s 'Finished r0.4 jobs'


         #############                                  
         ##  r0.5   ##
         #############
#REDOING SIMULATIONS FOR 160 and 390. REMOVE LINES if [[ "T160 T390" != *${T}* ]];then WHEN FINISHED
r='r0.5'; WD=$PROJD/8styrene32/r0.5
for temp in `seq 90 10 490`;do
  T=T$temp; #echo -e "\n##############\n T = $T\n##########"
if [[ "T160 T390" != *${T}* ]];then
  cd $WD/$T
  #ln -s ../8styrene32.pdb pdb
  #ln -s ../8styrene32.pdb .
  if [[ "T320 T330 T340 T350 T360 T370 T380 T390 T400 T410 T420 T430 T440 T450 T460 T470 T480 T490" == *${T}* ]];then
    nframes=100000
    #vmd -dispdev text -eofexit -e $PROJD/vmd/coalesce_trajectories.tcl; sleep 9s
    #cpptraj -p pdb -i $PROJD/amber/radial_8styrene32_105001.cpptraj
    #python $PROJD/python/CMAnalysis.py 8styrene32.pdb equil.dcd 5001 105001 100 'not resnum 1:256 and name C' 6.11 contacts_tolC_tolC.h5
    #python $PROJD/python/CMAnalysis.py 8styrene32.pdb equil.dcd 5001 105001 100 'not resnum 1:256 and name C' 7.0 contacts_tolC_tolC_7.0.h5
  else
    nframes=10000
    #ln -s equil.1.dcd equil.dcd
    #cpptraj -p pdb -i $PROJD/amber/radial_8styrene32_15001.cpptraj
    #python $PROJD/python/CMAnalysis.py 8styrene32.pdb equil.dcd 5001 15001 10 'not resnum 1:256 and name C' 6.11 contacts_tolC_tolC.h5
    #python $PROJD/python/CMAnalysis.py 8styrene2.pdb equil.dcd 5001 15001 10 'not resnum 1:256 and name C' 7.0 contacts_tolC_tolC_7.0.h5
  fi
  #echo -n "$T "; dumpdcd equil.dcd|head -1
  
  python $PROJD/python/CMsystemTrajectory.py 8styrene32.pdb equil.dcd cm.pdb cm.crd  #get trajectory of the CM of the whole system
  vmd -dispdev text -e $PROJD/vmd/crop_equil.tcl # Remove first 5ns, creating equil_cropped.dcd
  #vmd -dispdev win -eofexit -e $PROJD/vmd/rms2first.tcl; sleep 9s
  #echo -n "$T "; dumpdcd equil_rms2first.dcd | head -1
  #python $PROJD/python/diffusion_t0average.py pdb equil_cropped.dcd $nframes 1.0 8000 100 '(!:1-256)&(@H1,@H2,@H3,@H4,@H5)' diffusion_styd3.dat --rms2t0 no & # MSD(t) of the hydrogens in toluene rings
  #echo -n "$T "; python $PROJD/python/benedettoMSD.py diffusion_styd3.dat BASIS #BASIS MSD according to Benedetto12 reference
  #echo -n "$T "; python $PROJD/python/benedettoMSD.py diffusion_styd3.dat HFBS #HFBS MSD according to Benedetto12 reference
  #python $PROJD/python/CMtraj.py pdb equil_cropped.dcd '(!:1-256)&(@H1,@H2,@H3,@H4,@H5)' styd3CoM.pdb styd3CoM.dcd --aggregate byres & #trajectory and topology of Center of Masses of the hydrogens of toluene rings
  #python $PROJD/python/diffusion_t0average.py styd3CoM.pdb styd3CoM.dcd $nframes 1.0 8000 100 '@*' diffusion_styd3CoM.dat --rms2t0 no & #MSD(t) of Center of Masses of the hydrogens in toluene rings
  #echo -n "$T "; python $PROJD/python/benedettoMSD.py diffusion_styd3CoM.dat HFBS #BASIS MSD according to Benedetto12 reference
  #echo -n "$T "; python $PROJD/python/benedettoMSD.py diffusion_styd3CoM.dat BASIS #BASIS MSD according to Benedetto12 reference
  #python $PROJD/python/internalMSD.py diffusion_styd3.dat diffusion_styd3CoM.dat diffusion_styd3int.dat
  #echo -n "$T "; python $PROJD/python/benedettoMSD.py diffusion_styd3int.dat HFBS
  #echo -n "$T "; python $PROJD/python/benedettoMSD.py diffusion_styd3int.dat BASIS
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-256)&(@H*)' atomicfluct_4.56ns_styH.dat --seriesfile atomicfluct_4.56ns_styH_series.dat --msd yes #atomic fluctuations of all hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styH.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-256)&(@H1,@H2,@H3,@H4,@H5)' atomicfluct_4.56ns_styd3.dat --seriesfile atomicfluct_4.56ns_styd3_series.dat --msd yes #atomic fluctuations of all styrene ring hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styd3.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-256)&(@H62,@H71,@H72,@H73)' atomicfluct_4.56ns_styd5.dat --seriesfile atomicfluct_4.56ns_styd5_series.dat --msd yes #atomic fluctuations of all styrene backbone hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styd5.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(!(:1-256))&(@H=)' atomicfluct_4.56ns_tol.dat --seriesfile  atomicfluct_4.56ns_tol_series.dat --msd yes #atomic fluctuations of all toluene hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_tol.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(!(:1-256))&(@C)' atomicfluct_4.56ns_tolC.dat --seriesfile  atomicfluct_4.56ns_tolC_series.dat --msd yes #atomic fluctuations of central carobon of toluene in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_tolC.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(!(:1-256))&(@H1,@H2,@H3,@H4,@H5)' atomicfluct_4.56ns_told3.dat --seriesfile  atomicfluct_4.56ns_told3_series.dat --msd yes #atomic fluctuations of all toluene hydrogens in the time-scale of HFBS
  #python $PROJD/python/cluster_contacts.py contacts_tolC_tolC.h5 257-384 contacts_tolC_tolC_cluster_sizes.dat 
  #echo -n "$temp "; head -2 contacts_tolC_tolC_cluster_sizes.dat|tail --lines=1|tr '#' ' '
  #python $PROJD/python/cluster_contacts.py contacts_tolC_tolC_7.0.h5 257-384 contacts_tolC_tolC_7.0_cluster_sizes.dat
  #echo -n "$temp "; head -2 contacts_tolC_tolC_7.0_cluster_sizes.dat|tail --lines=1|tr '#' ' '
fi
done
popwindow.py 1s 'Finished r0.5 jobs'

  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 1220 100 '(:1-256)&(@H*)' atomicfluct_1.22ns_styH.dat --seriesfile atomicfluct_1.22ns_styH_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_1.22ns_styH.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 7640 100 '(:1-256)&(@H*)' atomicfluct_7.64ns_styH.dat --seriesfile atomicfluct_7.64ns_styH_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_7.64ns_styH.dat
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 1073 60 '(:1-256)&(@H*)' atomicfluct_BASIS_sty.dat --seriesfile  atomicfluct_BASIS_sty_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_BASIS_sty.dat
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 1073 60 '(:1-256)&(@H61,H62,H71,H72,H73)' atomicfluct_BASIS_bk.dat --seriesfile  atomicfluct_BASIS_bk_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_BASIS_bk.dat
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 1073 60 '(:1-256)&(@H1,H2,H3,H4,H5)' atomicfluct_BASIS_sc.dat --seriesfile  atomicfluct_BASIS_sc_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_BASIS_sc.dat
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 1073 60 '(!(:1-256))&(@H=)' atomicfluct_BASIS_tol.dat --seriesfile  atomicfluct_BASIS_tol_series.dat --msd yes & 
  #echo -n "$T "; cat atomicfluct_BASIS_tol.dat
  #python $PROJD/python/diffusion_t0average.py pdb equil_rms2first.dcd $nframes 1.0 5000 100 '(:1-256)&(@H*)' diffusion_sty.dat &
  #python $PROJD/python/diffusion_t0average.py pdb equil_rms2first.dcd $nframes 1.0 5000 100 '(!(:1-256))&(@H*)' diffusion_tol.dat &
  #python $PROJD/python/residueCMtrajectory.py pdb equil_rms2first.dcd 384 residueCM.pdb equil_rms2first_residueCM.crd &
  #python $PROJD/python/diffusion_t0average.py residueCM.pdb equil_rms2first_residueCM.crd $nframes 1.0 5000 100 ':1-256' diffusion_styCM.dat &  
  #time python $PROJD/python/atomic_fluct_t0average.py residueCM.pdb equil_rms2first_residueCM.crd $nframes 1073 60 ':1-256' atomicfluct_BASIS_styCM.dat --seriesfile atomicfluct_BASIS_styCM_series.dat --msd yes &
  #python $PROJD/python/CMtrajectory.py pdb equil_rms2first.dcd 32 8 CM.pdb equil_rms2first_CM.crd
  #time python $PROJD/python/atomic_fluct_t0average.py CM.pdb equil_rms2first_CM.crd $nframes 1073 99 ':1-8' atomicfluct_BASIS_CM.dat --seriesfile atomicfluct_BASIS_CM_series.dat --msd yes &
#  python $PROJD/python/diffusion_t0average.py CM.pdb equil_rms2first_CM.crd $nframes 1.0 5000 100 ':1-8' diffusion_CM.dat &

  #python $PROJD/python/rdf.py pdb equil_rms2first.dcd C3 C --npol 8 --spacing 0.1 --maximum 30.0 &
  #python $PROJD/python/rdf.py pdb equil_rms2first.dcd C7 C --npol 8 --spacing 0.1 --maximum 30.0 &
  #python $PROJD/python/toluene_shell.py pdb equil_rms2first.dcd C3 C toluene_shell_around_C3.dat --npol 8 --spacing 0.1 --minimum 3.0 --maximum 9.0 &
  #python $PROJD/python/toluene_shell.py pdb equil_rms2first.dcd C7 C toluene_shell_around_C7.dat --npol 8 --spacing 0.1 --minimum 3.0 --maximum 9.0 &
  #python $PROJD/python/thermo_averages.py  equil.log
  #for i in `seq 1 8`;do for f in backbone sidechain;do python $PROJD/python/generate_mwcovar.py $i $f;done;done
done


         #############                                  
         ##   r1    ##
         #############
r='r1'; WD=$PROJD/8styrene32/r1
for temp in `seq 90 10 460`;do
  T=T$temp; echo -e "\n##############\n T = $T\n##########"
if [[ "T100 T160 T410" != *${T}* ]];then
  cd $WD/$T
  #ln -s ../8styrene32.pdb pdb
  #ln -s ../8styrene32.pdb .
  nframes=10000
  #echo -n "$T "; dumpdcd equil.1.dcd | head -1
  #ln -s equil.1.dcd equil.dcd
 #vmd -dispdev win -eofexit -e $PROJD/vmd/rms2first.tcl; sleep 9s
  #echo -n "$T "; dumpdcd equil_rms2first.dcd | head -1
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-256)&(@H*)' atomicfluct_4.56ns_styH.dat --seriesfile atomicfluct_4.56ns_styH_series.dat --msd yes #atomic fluctuations of all hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styH.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-256)&(@H1,@H2,@H3,@H4,@H5)' atomicfluct_4.56ns_styd3.dat --seriesfile atomicfluct_4.56ns_styd3_series.dat --msd yes #atomic fluctuations of all styrene ring hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styd3.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-256)&(@H62,@H71,@H72,@H73)' atomicfluct_4.56ns_styd5.dat --seriesfile atomicfluct_4.56ns_styd5_series.dat --msd yes #atomic fluctuations of all styrene backbone hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styd5.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(!(:1-256))&(@H=)' atomicfluct_4.56ns_tol.dat --seriesfile  atomicfluct_4.56ns_tol_series.dat --msd yes #atomic fluctuations of all toluene hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_tol.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(!(:1-256))&(@C)' atomicfluct_4.56ns_tolC.dat --seriesfile  atomicfluct_4.56ns_tolC_series.dat --msd yes #atomic fluctuations of central carobon of toluene in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_tolC.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(!(:1-256))&(@H1,@H2,@H3,@H4,@H5)' atomicfluct_4.56ns_told3.dat --seriesfile  atomicfluct_4.56ns_told3_series.dat --msd yes #atomic fluctuations of all toluene hydrogens in the time-scale of HFBS
  #cpptraj -p pdb -i $PROJD/amber/radial_8styrene32_15001.cpptraj
  #python $PROJD/python/CMAnalysis.py 8styrene32.pdb equil.dcd 5001 15001 10 'not resnum 1:256 and name C' 6.11 contacts_tolC_tolC.h5
    #python $PROJD/python/CMAnalysis.py 8styrene32.pdb equil.dcd 5001 15001 10 'not resnum 1:256 and name C' 7.0 contacts_tolC_tolC_7.0.h5
  #python $PROJD/python/cluster_contacts.py contacts_tolC_tolC.h5 257-512 contacts_tolC_tolC_cluster_sizes.dat 
  #echo -n "$temp "; head -2 contacts_tolC_tolC_cluster_sizes.dat|tail --lines=1|tr '#' ' '
  #python $PROJD/python/cluster_contacts.py contacts_tolC_tolC_7.0.h5 257-512 contacts_tolC_tolC_7.0_cluster_sizes.dat
  #echo -n "$temp "; head -2 contacts_tolC_tolC_7.0_cluster_sizes.dat|tail --lines=1|tr '#' ' '
fi
done
  #ln -s ../8styrene32.pdb pdb
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1220 100 '(:1-256)&(@H*)' atomicfluct_1.22ns_styH.dat --seriesfile atomicfluct_1.22ns_styH_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_1.22ns_styH.dat #;  echo ''
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 7640 100 '(:1-256)&(@H*)' atomicfluct_7.64ns_styH.dat --seriesfile atomicfluct_7.64ns_styH_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_7.64ns_styH.dat #;  echo ''
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1073 30 '(:1-256)&(@H*)' atomicfluct_BASIS_sty.dat --seriesfile  atomicfluct_BASIS_sty_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_BASIS_sty.dat
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1073 30 '(:1-256)&(@H61,H62,H71,H72,H73)' atomicfluct_BASIS_bk.dat --seriesfile  atomicfluct_BASIS_bk_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_BASIS_bk.dat
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1073 30 '(:1-256)&(@H1,H2,H3,H4,H5)' atomicfluct_BASIS_sc.dat --seriesfile  atomicfluct_BASIS_sc_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_BASIS_sc.dat
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1073 30 '(!(:1-256))&(@H=)' atomicfluct_BASIS_tol.dat --seriesfile  atomicfluct_BASIS_tol_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_BASIS_tol.dat
  #python $PROJD/python/diffusion_t0average.py pdb equil_rms2first.dcd 10000 1.0 5000 100 '(:1-256)&(@H*)' diffusion_sty.dat &
  #python $PROJD/python/diffusion_t0average.py pdb equil_rms2first.dcd 10000 1.0 5000 100 '(!(:1-256))&(@H*)' diffusion_tol.dat &
  #python $PROJD/python/residueCMtrajectory.py pdb equil_rms2first.dcd 512 residueCM.pdb equil_rms2first_residueCM.crd &
  #python $PROJD/python/diffusion_t0average.py residueCM.pdb equil_rms2first_residueCM.crd 10000 1.0 5000 100 ':1-256' diffusion_styCM.dat &
  #time python $PROJD/python/atomic_fluct_t0average.py residueCM.pdb equil_rms2first_residueCM.crd 10000 1073 60 ':1-256' atomicfluct_BASIS_styCM.dat --seriesfile atomicfluct_BASIS_styCM_series.dat --msd yes &
  #python $PROJD/python/CMtrajectory.py pdb equil_rms2first.dcd 32 8 CM.pdb equil_rms2first_CM.crd
  #time python $PROJD/python/atomic_fluct_t0average.py CM.pdb equil_rms2first_CM.crd 10000 1073 99 ':1-8' atomicfluct_BASIS_CM.dat --seriesfile atomicfluct_BASIS_CM_series.dat --msd yes &
  #python $PROJD/python/diffusion_t0average.py CM.pdb equil_rms2first_CM.crd 10000 1.0 5000 100 ':1-8' diffusion_CM.dat &

  #python $PROJD/python/rdf.py pdb equil_rms2first.dcd C3 C --npol 8 --spacing 0.1 --maximum 30.0 &
  #python $PROJD/python/rdf.py pdb equil_rms2first.dcd C7 C --npol 8 --spacing 0.1 --maximum 30.0 &
  #python $PROJD/python/toluene_shell.py pdb equil_rms2first.dcd C3 C toluene_shell_around_C3.dat --npol 8 --spacing 0.1 --minimum 3.0 --maximum 9.0 &
  #python $PROJD/python/toluene_shell.py pdb equil_rms2first.dcd C7 C toluene_shell_around_C7.dat --npol 8 --spacing 0.1 --minimum 3.0 --maximum 9.0 &
  #python $PROJD/python/thermo_averages.py  equil.log
  #for i in `seq 1 8`;do for f in backbone sidechain;do python $PROJD/python/generate_mwcovar.py $i $f;done;done
done


         #############                                  
         ##   r2    ##
         #############
r='r2'; WD=$PROJD/8styrene32/r2
for temp in `seq 90 10 460`;do
  T=T$temp; #echo -e "\n##############\n T = $T\n##########"
if [[ "T100 T110 T300" != *${T}* ]];then
  cd $WD/$T
  nframes=10000
  #ln -s ../8styrene32.pdb .;   ln -s ../8styrene32.pdb pdb;   ln -s equil.1.dcd equil.dcd
  #echo -n "$T "; dumpdcd equil.dcd|head -1
  #vmd -dispdev win -eofexit -e $PROJD/vmd/rms2first.tcl; sleep 9s
  #echo -n "$T "; dumpdcd equil_rms2first.dcd | head -1
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-256)&(@H*)' atomicfluct_4.56ns_styH.dat --seriesfile atomicfluct_4.56ns_styH_series.dat --msd yes #atomic fluctuations of all hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styH.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-256)&(@H1,@H2,@H3,@H4,@H5)' atomicfluct_4.56ns_styd3.dat --seriesfile atomicfluct_4.56ns_styd3_series.dat --msd yes #atomic fluctuations of all styrene ring hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styd3.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-256)&(@H62,@H71,@H72,@H73)' atomicfluct_4.56ns_styd5.dat --seriesfile atomicfluct_4.56ns_styd5_series.dat --msd yes #atomic fluctuations of all styrene backbone hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styd5.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(!(:1-256))&(@H=)' atomicfluct_4.56ns_tol.dat --seriesfile  atomicfluct_4.56ns_tol_series.dat --msd yes #atomic fluctuations of all toluene hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_tol.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(!(:1-256))&(@C)' atomicfluct_4.56ns_tolC.dat --seriesfile  atomicfluct_4.56ns_tolC_series.dat --msd yes #atomic fluctuations of central carobon of toluene in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_tolC.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(!(:1-256))&(@H1,@H2,@H3,@H4,@H5)' atomicfluct_4.56ns_told3.dat --seriesfile  atomicfluct_4.56ns_told3_series.dat --msd yes #atomic fluctuations of all toluene hydrogens in the time-scale of HFBS
  #cpptraj -p pdb -i $PROJD/amber/radial_8styrene32_15001.cpptraj
  #python $PROJD/python/CMAnalysis.py 8styrene32.pdb equil.dcd 5001 15001 10 'not resnum 1:256 and name C' 6.11 contacts_tolC_tolC.h5
  #python $PROJD/python/CMAnalysis.py 8styrene32.pdb equil.dcd 5001 15001 10 'not resnum 1:256 and name C' 7.0 contacts_tolC_tolC_7.0.h5
  #python $PROJD/python/cluster_contacts.py contacts_tolC_tolC.h5 257-768 contacts_tolC_tolC_cluster_sizes.dat 
  #echo -n "$temp "; head -2 contacts_tolC_tolC_cluster_sizes.dat|tail --lines=1|tr '#' ' '
  #python $PROJD/python/cluster_contacts.py contacts_tolC_tolC_7.0.h5 257-768 contacts_tolC_tolC_7.0_cluster_sizes.dat
fi
done
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1220 100 '(:1-256)&(@H*)' atomicfluct_1.22ns_styH.dat --seriesfile atomicfluct_1.22ns_styH_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_1.22ns_styH.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 7640 100 '(:1-256)&(@H*)' atomicfluct_7.64ns_styH.dat --seriesfile atomicfluct_7.64ns_styH_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_7.64ns_styH.dat
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1073 30 '(:1-256)&(@H*)' atomicfluct_BASIS_sty.dat --seriesfile  atomicfluct_BASIS_sty_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_BASIS_sty.dat
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1073 30 '(:1-256)&(@H61,H62,H71,H72,H73)' atomicfluct_BASIS_bk.dat --seriesfile  atomicfluct_BASIS_bk_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_BASIS_bk.dat
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1073 30 '(:1-256)&(@H1,H2,H3,H4,H5)' atomicfluct_BASIS_sc.dat --seriesfile  atomicfluct_BASIS_sc_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_BASIS_sc.dat
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1073 30 '(!(:1-256))&(@H=)' atomicfluct_BASIS_tol.dat --seriesfile  atomicfluct_BASIS_tol_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_BASIS_tol.dat
  #python $PROJD/python/diffusion_t0average.py pdb equil_rms2first.dcd 10000 1.0 5000 100 '(:1-256)&(@H*)' diffusion_sty.dat &
  #python $PROJD/python/diffusion_t0average.py pdb equil_rms2first.dcd 10000 1.0 5000 100 '(!(:1-256))&(@H*)' diffusion_tol.dat &
  #python $PROJD/python/residueCMtrajectory.py pdb equil_rms2first.dcd 768 residueCM.pdb equil_rms2first_residueCM.crd &
  #python $PROJD/python/diffusion_t0average.py residueCM.pdb equil_rms2first_residueCM.crd 10000 1.0 5000 100 ':1-256' diffusion_styCM.dat &
  #time python $PROJD/python/atomic_fluct_t0average.py residueCM.pdb equil_rms2first_residueCM.crd 10000 1073 60 ':1-256' atomicfluct_BASIS_styCM.dat --seriesfile atomicfluct_BASIS_styCM_series.dat --msd yes &
  #python $PROJD/python/CMtrajectory.py pdb equil_rms2first.dcd 32 8 CM.pdb equil_rms2first_CM.crd
  #time python $PROJD/python/atomic_fluct_t0average.py CM.pdb equil_rms2first_CM.crd 10000 1073 99 ':1-8' atomicfluct_BASIS_CM.dat --seriesfile atomicfluct_BASIS_CM_series.dat --msd yes &
  #python $PROJD/python/diffusion_t0average.py CM.pdb equil_rms2first_CM.crd 10000 1.0 5000 100 ':1-8' diffusion_CM.dat &
  #python $PROJD/python/rdf.py pdb equil_rms2first.dcd C3 C --npol 8 --spacing 0.1 --maximum 30.0 &
  #python $PROJD/python/rdf.py pdb equil_rms2first.dcd C7 C --npol 8 --spacing 0.1 --maximum 30.0 &
  #python $PROJD/python/toluene_shell.py pdb equil_rms2first.dcd C3 C toluene_shell_around_C3.dat --npol 8 --spacing 0.1 --minimum 3.0 --maximum 9.0 &
  #python $PROJD/python/toluene_shell.py pdb equil_rms2first.dcd C7 C toluene_shell_around_C7.dat --npol 8 --spacing 0.1 --minimum 3.0 --maximum 9.0 &
  #python $PROJD/python/thermo_averages.py  equil.log
  #for i in `seq 1 8`;do for f in backbone sidechain;do python $PROJD/python/generate_mwcovar.py $i $f;done;done
done


         #########################                                  
         ##   SINGLE POLYMER    ##
         #########################
sleep 5h
WD=$PROJD/1styrene32/1styrene32_4017toluene
for temp in 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 230 250 270 290 310 330 340 350 360;do
  T="T$temp"; cd $WD/$T;
  #ln -s ../1styrene32.pdb pdb; ln -s ../1styrene32.pdb .
  nframes=10000
  #echo -n "$T "; dumpdcd equil.dcd|head -1
  #vmd -dispdev win -eofexit -e $PROJD/vmd/rms2first.tcl #coalesce trajectories
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-32)&(@H*)' atomicfluct_4.56ns_styH.dat --seriesfile atomicfluct_4.56ns_styH_series.dat --msd yes #atomic fluctuations of all hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styH.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-32)&(@H1,@H2,@H3,@H4,@H5)' atomicfluct_4.56ns_styd3.dat --seriesfile atomicfluct_4.56ns_styd3_series.dat --msd yes #atomic fluctuations of all styrene ring hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styd3.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 100 '(:1-32)&(@H62,@H71,@H72,@H73)' atomicfluct_4.56ns_styd5.dat --seriesfile atomicfluct_4.56ns_styd5_series.dat --msd yes #atomic fluctuations of all styrene backbone hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_styd5.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 30 '(!(:1-32))&(@H=)' atomicfluct_4.56ns_tol.dat --seriesfile  atomicfluct_4.56ns_tol_series.dat --msd yes #atomic fluctuations of all toluene hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_tol.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 30 '(!(:1-32))&(@C)' atomicfluct_4.56ns_tolC.dat --seriesfile  atomicfluct_4.56ns_tolC_series.dat --msd yes &#atomic fluctuations of central carobon of toluene in the time-scale of HFBS
  echo -n "$T "; cat atomicfluct_4.56ns_tolC.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 4560 30 '(!(:1-32))&(@H1,@H2,@H3,@H4,@H5)' atomicfluct_4.56ns_told3.dat --seriesfile  atomicfluct_4.56ns_told3_series.dat --msd yes & #atomic fluctuations of all toluene hydrogens in the time-scale of HFBS
  #echo -n "$T "; cat atomicfluct_4.56ns_told3.dat
  #cpptraj -p pdb -i $PROJD/amber/radial_1styrene32_15001.cpptraj
  #python $PROJD/python/CMAnalysis.py 1styrene32.pdb equil.dcd 5001 15001 10 'not resnum 1:32 and name C' 6.11 contacts_tolC_tolC.h5 &
  #python $PROJD/python/CMAnalysis.py 1styrene32.pdb equil.dcd 5001 15001 10 'not resnum 1:32 and name C' 7.0 contacts_tolC_tolC_7.0.h5
  #python $PROJD/python/cluster_contacts.py contacts_tolC_tolC.h5 33-4139 contacts_tolC_tolC_cluster_sizes.dat 
done
  #python $PROJD/python/cluster_trajectory.py ./pdb ./equil_rms2first.dcd 500 '(:1-32)&(@H*)' centroid_styH.pdb average_styH.pdb
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1220 100 '(:1-32)&(@H*)' atomicfluct_1.22ns_styH.dat --seriesfile atomicfluct_1.22ns_styH_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_1.22ns_styH.dat;
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 7640 100 '(:1-32)&(@H*)' atomicfluct_7.64ns_styH.dat --seriesfile atomicfluct_7.64ns_styH_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_7.64ns_styH.dat;
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1073 30 '(:1-32)&(@H*)' atomicfluct_BASIS_sty.dat --seriesfile atomicfluct_BASIS_sty_series.dat --msd yes &
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1073 30 '(:1-32)&(@H61,H62,H71,H72,H73)' atomicfluct_BASIS_bk.dat --seriesfile  atomicfluct_BASIS_bk_series.dat --msd yes &
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1073 30 '(:1-32)&(@H1,H2,H3,H4,H5)' atomicfluct_BASIS_sc.dat --seriesfile  atomicfluct_BASIS_sc_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_BASIS_sc.dat
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1073 30 '(!(:1-32))&(@H=)' atomicfluct_BASIS_tol.dat --seriesfile  atomicfluct_BASIS_tol_series.dat --msd yes &
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 10 60 '(:1-32)&(@H*)' atomicfluct_10ps_sty.dat --seriesfile atomicfluct_10ps_sty_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_10ps_sty.dat
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 100 60 '(:1-32)&(@H*)' atomicfluct_100ps_sty.dat --seriesfile atomicfluct_100ps_sty_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_100ps_sty.dat
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 5000 60 '(:1-32)&(@H*)' atomicfluct_5ns_sty.dat --seriesfile atomicfluct_5ns_sty_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_5ns_sty.dat
  #python $PROJD/python/diffusion_t0average.py pdb equil_rms2first.dcd 10000 1.0 5000 100 '(:1-32)&(@H*)' diffusion_sty.dat &
  #python $PROJD/python/diffusion_t0average.py pdb equil_rms2first.dcd 10000 1.0 5000 100 '(!(:1-32))&(@H*)' diffusion_tol.dat &
  #python $PROJD/python/residueCMtrajectory.py pdb equil_rms2first.dcd 4139 residueCM.pdb equil_rms2first_residueCM.crd &
  #python $PROJD/python/diffusion_t0average.py residueCM.pdb equil_rms2first_residueCM.crd 10000 1.0 5000 100 ':1-32' diffusion_styCM.dat &
  #time python $PROJD/python/atomic_fluct_t0average.py residueCM.pdb equil_rms2first_residueCM.crd 10000 1073 60 ':1-32' atomicfluct_BASIS_styCM.dat --seriesfile atomicfluct_BASIS_styCM_series.dat --msd yes &
  #python $PROJD/python/CMtrajectory.py pdb equil_rms2first.dcd 32 2 CM.pdb equil_rms2first_CM.crd &#Note: we take the first 32 toluenes as another polymer. The only reason is that cpptraj will not load a pdb file containing only one atom. Thus, CM.pdb need to contain two atoms
  #time python $PROJD/python/atomic_fluct_t0average.py CM.pdb equil_rms2first_CM.crd 10000 1073 99 ':1' atomicfluct_BASIS_CM.dat --seriesfile atomicfluct_BASIS_CM_series.dat --msd yes &
  #python $PROJD/python/diffusion_t0average.py CM.pdb equil_rms2first_CM.crd 10000 1.0 5000 100 ':1' diffusion_CM.dat --rms2t0 no & # we can't do rsm to first snapshot of each chunk because our system is made up of only one atom.
  #python $PROJD/python/rdf.py pdb equil_rms2first.dcd C3 C --npol 1 --spacing 0.1 --maximum 30.0 &
  #python $PROJD/python/rdf.py pdb equil_rms2first.dcd C7 C --npol 1 --spacing 0.1 --maximum 30.0 &
  #time python $PROJD/python/toluene_shell.py pdb equil_rms2first.dcd C3 C toluene_shell_around_C3.dat --npol 1 --spacing 0.1 --minimum 3.0 --maximum 9.0 &
  #time python $PROJD/python/toluene_shell.py pdb equil_rms2first.dcd C7 C toluene_shell_around_C7.dat --npol 1 --spacing 0.1 --minimum 3.0 --maximum 9.0 &
  #python $PROJD/python/thermo_averages.py  equil.log
  #for flavor in backbone sidechain;do python $PROJD/python/generate_mwcovar.py 1 $flavor;done
done
popwindow.py 1s 'Finished pStyrene_Toluene bash_commands.sh'
