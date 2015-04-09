#!/bin/bash

PROJD=/SNSlocal/projects/jbq/pStyrene_Toluene
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
#                                              #############
#                                              ###  MSD  ###
#                                              #############
                                          
# r0
WD=$PROJD/8styrene32/r0
for T in T150 T170 T190 T210 T220 T230 T250 T270 T290 T310 T330 T360 T370 T380 T390 T400 T405 T410 T420 T430 T440 T450 T460 T470;do
  cd $WD/$T
  #ln -s ../8styrene32.pdb pdb
  #nframes=10000
  #if [[ "T380 T390 T400 T405 T410 T420 T430 T440 T450 T460 T470" == *${T}* ]];then nframes=30000;fi
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 1220 100 '(:1-256)&(@H*)' atomicfluct_1.22ns_styH.dat --seriesfile atomicfluct_1.22ns_styH_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_1.22ns_styH.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 7640 100 '(:1-256)&(@H*)' atomicfluct_7.64ns_styH.dat --seriesfile atomicfluct_7.64ns_styH_series.dat --msd yes &
  echo -n "$T "; cat atomicfluct_7.64ns_styH.dat
  #python $PROJD/python/thermo_averages.py  equil.log
  #for i in `seq 1 8`;do for f in backbone sidechain;do python $PROJD/python/generate_mwcovar.py $i $f;done;done
done

# r0.5
WD=$PROJD/8styrene32/r0.5
for  T in T270 T280 T290 T310 T330 T350 T370 T380 T400 T410 T420 T430 T440 T450;do
  cd $WD/$T
  #ln -s ../8styrene32.pdb pdb
  #nframes=10000
  #if [[ "T400 T410 T420 T430 T440 T450" == *${T}* ]];then nframes=30000;fi
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 1220 100 '(:1-256)&(@H*)' atomicfluct_1.22ns_styH.dat --seriesfile atomicfluct_1.22ns_styH_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_1.22ns_styH.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd $nframes 7640 100 '(:1-256)&(@H*)' atomicfluct_7.64ns_styH.dat --seriesfile atomicfluct_7.64ns_styH_series.dat --msd yes &
  echo -n "$T "; cat atomicfluct_7.64ns_styH.dat
  #python $PROJD/python/rdf.py pdb equil_rms2first.dcd C3 C --npol 8 --spacing 0.1 --maximum 30.0 &
  #python $PROJD/python/rdf.py pdb equil_rms2first.dcd C7 C --npol 8 --spacing 0.1 --maximum 30.0 &
  #python $PROJD/python/toluene_shell.py pdb equil_rms2first.dcd C3 C toluene_shell_around_C3.dat --npol 8 --spacing 0.1 --minimum 3.0 --maximum 9.0 &
  #python $PROJD/python/toluene_shell.py pdb equil_rms2first.dcd C7 C toluene_shell_around_C7.dat --npol 8 --spacing 0.1 --minimum 3.0 --maximum 9.0 &
  #python $PROJD/python/thermo_averages.py  equil.log
  #for i in `seq 1 8`;do for f in backbone sidechain;do python $PROJD/python/generate_mwcovar.py $i $f;done;done
done

# r1
WD=$PROJD/8styrene32/r1
for  T in T230 T240 T250 T260 T270 T270 T280 T290 T300 T310 T320 T330 T340 T350 T360 T370 T380 T390 T400 T420 T430 T440;do
  cd $WD/$T
  #ln -s ../8styrene32.pdb pdb
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1220 100 '(:1-256)&(@H*)' atomicfluct_1.22ns_styH.dat --seriesfile atomicfluct_1.22ns_styH_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_1.22ns_styH.dat #;  echo ''
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 7640 100 '(:1-256)&(@H*)' atomicfluct_7.64ns_styH.dat --seriesfile atomicfluct_7.64ns_styH_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_7.64ns_styH.dat #;  echo ''
  #python $PROJD/python/rdf.py pdb equil_rms2first.dcd C3 C --npol 8 --spacing 0.1 --maximum 30.0 &
  #python $PROJD/python/rdf.py pdb equil_rms2first.dcd C7 C --npol 8 --spacing 0.1 --maximum 30.0 &
  #python $PROJD/python/toluene_shell.py pdb equil_rms2first.dcd C3 C toluene_shell_around_C3.dat --npol 8 --spacing 0.1 --minimum 3.0 --maximum 9.0 &
  #python $PROJD/python/toluene_shell.py pdb equil_rms2first.dcd C7 C toluene_shell_around_C7.dat --npol 8 --spacing 0.1 --minimum 3.0 --maximum 9.0 &
  #python $PROJD/python/thermo_averages.py  equil.log
  #for i in `seq 1 8`;do for f in backbone sidechain;do python $PROJD/python/generate_mwcovar.py $i $f;done;done
done

# r2
WD=$PROJD/8styrene32/r2
for T in T170 T190 T210 T230 T250 T270 T290 T310 T330 T350 T370 T390 T410;do
  cd $WD/$T
  #ln -s ../8styrene32.pdb pdb
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1220 100 '(:1-256)&(@H*)' atomicfluct_1.22ns_styH.dat --seriesfile atomicfluct_1.22ns_styH_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_1.22ns_styH.dat
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 7640 100 '(:1-256)&(@H*)' atomicfluct_7.64ns_styH.dat --seriesfile atomicfluct_7.64ns_styH_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_7.64ns_styH.dat
  time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1073 30 '(:1-32)&(@H61,H62,H71,H72,H73)' atomicfluct_BASIS_bk.dat --seriesfile  atomicfluct_BASIS_bk_series.dat --msd yes &
  sleep 9s
  time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1073 300 '(:1-32)&(@H1,H2,H3,H4,H5)' atomicfluct_BASIS_sc.dat --seriesfile  atomicfluct_BASIS_sc_series.dat --msd yes &
  sleep 9s
  time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1073 30 '(!(:1-32))&(@H=)' atomicfluct_BASIS_tol.dat --seriesfile  atomicfluct_BASIS_tol_series.dat --msd yes &
  sleep 9s
  #python $PROJD/python/rdf.py pdb equil_rms2first.dcd C3 C --npol 8 --spacing 0.1 --maximum 30.0 &
  #python $PROJD/python/rdf.py pdb equil_rms2first.dcd C7 C --npol 8 --spacing 0.1 --maximum 30.0 &
  #python $PROJD/python/toluene_shell.py pdb equil_rms2first.dcd C3 C toluene_shell_around_C3.dat --npol 8 --spacing 0.1 --minimum 3.0 --maximum 9.0 &
  #python $PROJD/python/toluene_shell.py pdb equil_rms2first.dcd C7 C toluene_shell_around_C7.dat --npol 8 --spacing 0.1 --minimum 3.0 --maximum 9.0 &
  #python $PROJD/python/thermo_averages.py  equil.log
  #for i in `seq 1 8`;do for f in backbone sidechain;do python $PROJD/python/generate_mwcovar.py $i $f;done;done
done

# SINGLE POLYMER
WD=$PROJD/1styrene32/1styrene32_4017toluene
for T in T90 T100 T110 T120 T130 T140 T150 T160 T170 T180 T190 T200 T210 T230 T250 T270 T290 T310;do
  cd $WD/$T; #ln -s ../1styrene32.pdb pdb
  #vmd -dispdev win -eofexit -e $PROJD/vmd/rms2first.tcl #coalesce trajectories
  #python $PROJD/python/cluster_trajectory.py ./pdb ./equil_rms2first.dcd 500 '(:1-32)&(@H*)' centroid_styH.pdb average_styH.pdb
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1220 100 '(:1-32)&(@H*)' atomicfluct_1.22ns_styH.dat --seriesfile atomicfluct_1.22ns_styH_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_1.22ns_styH.dat;
  #python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 7640 100 '(:1-32)&(@H*)' atomicfluct_7.64ns_styH.dat --seriesfile atomicfluct_7.64ns_styH_series.dat --msd yes &
  #echo -n "$T "; cat atomicfluct_7.64ns_styH.dat;
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1073 30 '(:1-32)&(@H*)' atomicfluct_BASIS_styH.dat --seriesfile atomicfluct_BASIS_styH_series.dat --msd yes &
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1073 30 '(:1-32)&(@H61,H62,H71,H72,H73)' atomicfluct_BASIS_bk.dat --seriesfile  atomicfluct_BASIS_bk_series.dat --msd yes &
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1073 300 '(:1-32)&(@H1,H2,H3,H4,H5)' atomicfluct_BASIS_sc.dat --seriesfile  atomicfluct_BASIS_sc_series.dat --msd yes &
  #time python $PROJD/python/atomic_fluct_t0average.py pdb equil_rms2first.dcd 10000 1073 30 '(!(:1-32))&(@H=)' atomicfluct_BASIS_tol.dat --seriesfile  atomicfluct_BASIS_tol_series.dat --msd yes &
  #python $PROJD/python/rdf.py pdb equil_rms2first.dcd C3 C --npol 1 --spacing 0.1 --maximum 30.0 &
  #python $PROJD/python/rdf.py pdb equil_rms2first.dcd C7 C --npol 1 --spacing 0.1 --maximum 30.0 &
  #time python $PROJD/python/toluene_shell.py pdb equil_rms2first.dcd C3 C toluene_shell_around_C3.dat --npol 1 --spacing 0.1 --minimum 3.0 --maximum 9.0 &
  #time python $PROJD/python/toluene_shell.py pdb equil_rms2first.dcd C7 C toluene_shell_around_C7.dat --npol 1 --spacing 0.1 --minimum 3.0 --maximum 9.0 &
  #python $PROJD/python/thermo_averages.py  equil.log
  #for flavor in backbone sidechain;do python $PROJD/python/generate_mwcovar.py 1 $flavor;done
done
popwindow.py 1s 'Finished pStyrene_Toluene bash_commands.sh'
