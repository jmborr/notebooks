#!/bin/bash

PROJD=/SNSlocal/projects/jbq/POSS/dms8t8
for T in `seq 20 10 300`;do
  cd $PROJD/${T}K/
  #ln -s ../dms8t8.pdb
  #ln -s  ${T}K_rms2first.dcd dms8t8_rms2first.dcd
  #ln -s  ${T}K_rms2first.dcd rms2first.dcd
  cpptraj -i ../cpptraj/atomicfluct.cpptraj
  #python $PROJD/python/atomic_fluct_t0average.py dms8t8.pdb ${T}K_rms2first.dcd 10000 1073 100 '@H=' atomicfluct_BASIS_styH.dat --seriesfile  atomicfluct_BASIS_styH_series.dat --msd yes
  #DIHEDRAL O-Si-O-Si
  #mkdir -p dihedrals/O-Si-O-Si
  #for dihedral in O11-Si10-O19-Si15 O13-Si4-O10-Si8 O14-Si11-O15-Si12 O3-Si1-O7-Si5 O4-Si9-O16-Si13 O5-Si2-O8-Si6 O6-Si3-O9-Si7 O8-Si14-O20-Si16;do
  #  cpptraj -p dms8t8_nd.pdb -i ../cpptraj/dihedrals/O-Si-O-Si/${dihedral}.ptraj
     #python $PROJD/python/dihedral_autocorrelation.py --infile ./dihedrals/O-Si-O-Si/${dihedral}.dat --outfile ./dihedrals/O-Si-O-Si/${dihedral}.corr.dat
  #done

  #DIHEDRAL Si-O-Si-C
  #mkdir -p dihedrals/Si-O-Si-C
  #for dihedral in Si10-O19-Si15-C13 Si11-O15-Si12-C9 Si14-O20-Si16-C16 Si1-O7-Si5-C2 Si2-O8-Si6-C4 Si3-O9-Si7-C6 Si4-O10-Si8-C8 Si9-O16-Si13-C11;do
  #  cpptraj -p dms8t8_nd.pdb -i ../cpptraj/dihedrals/Si-O-Si-C/${dihedral}.ptraj
  #   python $PROJD/python/dihedral_autocorrelation.py --infile ./dihedrals/Si-O-Si-C/${dihedral}.dat --outfile ./dihedrals/Si-O-Si-C/${dihedral}.corr.dat
  #done

  #DIHEDRAL O-Si-C-H
  #mkdir -p dihedrals/O-Si-C-H
  #for dihedral in O10-Si8-C7-H7 O10-Si8-C8-H12 O15-Si12-C10-H47 O15-Si12-C9-H44 O16-Si13-C11-H14 O16-Si13-C12-H17 O19-Si15-C13-H38 O19-Si15-C14-H41 O20-Si16-C15-H21 O20-Si16-C16-H22 O7-Si5-C1-H25 O7-Si5-C2-H28 O8-Si6-C3-H31 O8-Si6-C4-H36 O9-Si7-C5-H1 O9-Si7-C6-H6;do
  #  cpptraj -p dms8t8_nd.pdb -i ../cpptraj/dihedrals/O-Si-C-H/${dihedral}.ptraj
    #python $PROJD/python/dihedral_autocorrelation.py --infile ./dihedrals/O-Si-C-H/${dihedral}.dat --outfile ./dihedrals/O-Si-C-H/${dihedral}.corr.dat
  #done

  #python $PROJD/python/average_dihedral_autocorrelation.py
  
done

for T in `seq 20 10 300`;do
  echo -n "${T}K/dihedrals/O-Si-O-Si/O-Si-O-Si.corr.dat "
done

