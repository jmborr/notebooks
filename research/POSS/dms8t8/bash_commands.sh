#!/bin/bash

PROJD=/SNSlocal/projects/jbq/POSS/dms8t8
for T in `seq 20 10 300`;do
  cd $PROJD/${T}K/
  #ln -s ../dms8t8.pdb
  #ln -s  ${T}K_rms2first.dcd dms8t8_rms2first.dcd
  mkdir -p dihedrals/O-Si-O-Si
  for dihedral in O11-Si10-O19-Si15 O13-Si4-O10-Si8 O14-Si11-O15-Si12 O3-Si1-O7-Si5 O4-Si9-O16-Si13 O5-Si2-O8-Si6 O6-Si3-O9-Si7 O8-Si14-O20-Si16;do
    cpptraj -p dms8t8_nd.pdb -i ../cpptraj/dihedrals/O-Si-O-Si/${dihedral}.ptraj
  done
done

