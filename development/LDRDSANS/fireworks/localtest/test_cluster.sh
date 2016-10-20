#!/bin/bash

# Load environment for the clustering program
source $MODULESHOME/init/bash
module load fast_protein_cluster

# Set up the scratch directory
WD="/tmp/test_cluster"
/bin/mkdir -p $WD

# Transfer input files
DATADIR=$FASTPROTEINCLUSTERHOME/test
cd $DATADIR
echo -ne "\nTransferring input configurations to working directory..."
/bin/cp -r list 1af7_ test_rmsd_sse3.sh $WD/
echo -e "finished"
sleep 5s

# Run the cluster script
echo -e "\nStarting clustering..."
cd $WD/
./test_rmsd_sse3.sh
echo -e "..finished. Output files created in $WD"
