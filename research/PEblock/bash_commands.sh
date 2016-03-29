PROJD=/SNSlocal/projects/jbq/PEblock

############################
# FIST BATCH OF SIMULATIONS
############################
ljdir1s=("Low_surf_conc" "Philic_BCP2" "Phobic_BCP2" "Strong_Phobic_BCP2")
ljdir2s=("lj10_prod" "lj12_prod" "lj13_prod" "lj11_prod")

for i in `seq 0 3`;do
  ljdir1=${ljdir1s[$i]}
  ljdir2=${ljdir2s[$i]}
  cd $PROJD/$ljdir1/$ljdir2
  #Use VMD: load dump.xyz and save the first frame as dump.pdb
  #/bin/cp /SNSlocal/scratch/zf4/Block_PE_Surf/$ljdir1/$ljdir2/?cions.dat .
  #grep ' H ' dump.pdb > poly1.pdb        # neutral monomers of PE
  #grep ' He ' dump.pdb > poly2.pdb       # phylic or phobic block of PE
  #grep ' Li ' dump.pdb > polyCharge.pdb  # cation monomers of PE
  #grep ' B ' dump.pdb >  tail.pdb        # surfactant tail
  #grep ' Be ' dump.pdb >  head.pdb       # surfactant head
  #grep ' C ' dump.pdb >  Ncions.pdb       # ? counter-ions
  #grep ' N ' dump.pdb >  Pcions.pdb       # ? counter-ions
  #perl -p -i -e 's/ He/ H /g' poly2.pdb
  #perl -p -i -e 's/ Li/ H /g' polyCharge.pdb
  #perl -p -i -e 's/ B / H /g' tail.pdb
  #perl -p -i -e 's/ Be/ H /g' head.pdb
  #perl -p -i -e 's/ C / H /g' Ncions.pdb
  #perl -p -i -e 's/ N / H /g' Pcions.pdb
  #perl -p -i -e 's/ C / H /g' Ncions.pdb
  #perl -p -i -e 's/ N / H /g' Pcions.pdb
  #FIND OUT DIMENSIONS OF THE PERIODIC BOX
  #create trajectory of the first Ncion in text format (file junk.txt(
  #cpptraj -p Ncions.pdb -i $PROJD/cpptraj/Ncions100.cpptraj
  #Plot junk.txt with xmgrace to figure out the XYZ dimensions of the periodic box,and save to file dimensions.dat
  #ln -s ../../dimensions.dat #it turns out dimensions are same for all ljdirs (X=Y=Z=100)
  for name in poly1 poly2 polyCharge tail head Ncions Pcions;do
    #CREATE DCD FILES
    #/bin/cp $PROJD/vmd/dat2dcd.tcl junk.tcl
    #perl -pi -e "s/_NAME_/$name/g" junk.tcl
    #vmd -dispdev text -eofexit -e junk.tcl
    #CREATE UNWRAPED DCD FILES
    #/bin/cp $PROJD/vmd/dat2unwrapeddcd.tcl junk.tcl
    #perl -pi -e "s/_NAME_/$name/g" junk.tcl
    #vmd -dispdev text -eofexit -e junk.tcl
    #REMOVE TRANSLATIONS AND ROTATIONS AT EACH ELEMENTARY STEP
    #/bin/cp $PROJD/cpptraj/rms2prev.cpptraj junk.cpptraj
    #perl -pi -e "s/_NAME_/$name/g" junk.cpptraj
    #cpptraj -p $name.pdb -i junk.cpptraj &> /dev/null
    #CALCULATE MSD FOR rms2prev TRAJECTORIES
    #python $PROJD/python/diffusion_t0average.py $name.pdb ${name}_rms2prev.dcd 20000 500 10000 20 '@*' diffusion_${name}_rms2prev.dat --rms2t0=no &> /dev/null
    #Open MantiPlot and run $PROJD/python/h5todat.py to save Sassena output in ASCII format.
    #python $PROJD/python/nonnegative.py fqt_inc_$name.dat fqt_incNN_$name.dat
    #CALCULATE STATIC STRUCTURE FACTOR
    #/bin/cp ../../sassena_fq.xml sassena_fq_$name.xml
    #perl -pi -e "s/_NAME_/$name/g" sassena_fq_$name.xml
    #/bin/cp ../../sassena_fq.pbs sassena_fq_$name.pbs
    #perl -pi -e "s/_NAME_/$name/g" sassena_fq_$name.pbs
  done
done

##############################
# SECOND BATCH OF SIMULATIONS
##############################
zf4Root="/SNSlocal/projects/zf4/Block_PE_surf"
subDirs=("Philic_BCP2/300NC/lj7_prod" "Philic_BCP2/200NC/lj7_prod" "Philic_BCP2/50NC/lj7_prod" "Phobic_BCP2/300NC/lj7_prod" "Phobic_BCP2/200NC/lj7_prod" "Phobic_BCP2/50NC/lj7_prod")
lastIndex=5 # six subdirectories to work with, this is the index of last directory if start counting from zero
for iDir in `seq 0 "$lastIndex"`;do
    #mkdir -p $PROJD/${subDirs[$iDir]}
    #Use VMD: load dump.xyz and save the first frame as dump.pdb
    cd $PROJD/${subDirs[$iDir]}/
    #grep ' H '  dump.pdb > poly1.pdb       # neutral monomers of PE
    #grep ' He ' dump.pdb > poly2.pdb       # phylic or phobic block of PE
    #grep ' Li ' dump.pdb > polyCharge.pdb  # cation monomers of PE
    #grep ' B '  dump.pdb > tail.pdb        # surfactant tail
    #grep ' Be ' dump.pdb > head.pdb        # surfactant head
    #grep ' C '  dump.pdb > Ncions.pdb      # counter-anions
    #grep ' N '  dump.pdb > Pcions.pdb      # counter-cations
    #perl -p -i -e 's/ He/ H /g' poly2.pdb       #change element type to hydrogen
    #perl -p -i -e 's/ Li/ H /g' polyCharge.pdb
    #perl -p -i -e 's/ B / H /g' tail.pdb
    #perl -p -i -e 's/ Be/ H /g' head.pdb
    #perl -p -i -e 's/ C / H /g' Ncions.pdb
    #perl -p -i -e 's/ N / H /g' Pcions.pdb
    #perl -p -i -e 's/ C / H /g' Ncions.pdb
    #perl -p -i -e 's/ N / H /g' Pcions.pdb
    for name in poly1 poly2 polyCharge tail head Ncions Pcions;do
	#CREATE DCD FILES
	  #/bin/cp $PROJD/vmd/xyz2unwrappeddcd.tcl junk.tcl
          # XYZDIR is Monojoy's directory where the xyz file is located
	  #python $PROJD/python/replace_keyword.py junk.tcl "_XYZDIR_" $zf4Root/${subDirs[$iDir]}
  	  #python $PROJD/python/replace_keyword.py junk.tcl "_NAME_" $name
	  #vmd -dispdev text -eofexit -e junk.tcl &
	  #/bin/rm junk.tcl
        #REMOVE TRANSLATIONS AND ROTATIONS AT EACH ELEMENTARY STEP
	  #/bin/cp $PROJD/cpptraj/rms2prev.cpptraj junk.cpptraj
	  #perl -pi -e "s/_NAME_/$name/g" junk.cpptraj
	  #cpptraj -p $name.pdb -i junk.cpptraj &> /dev/null  # Be sure you do "module load amber" in the terminal
	  #/bin/rm junk.cpptraj
	#CALCULATE MSD FOR rms2prev TRAJECTORIES
	  #python $PROJD/python/diffusion_t0average.py $name.pdb ${name}_rms2prev.dcd 10000 500 7500 20 '@*' diffusion_${name}_rms2prev.dat --rms2t0=no &> /dev/null
	#
    done
done

##########################################################
# SECOND BATCH OF SIMULATIONS ASSUMING WRAPPED COORDINATES
##########################################################
zf4Root="/SNSlocal/projects/zf4/Block_PE_surf"
subDirs=("Philic_BCP2/300NC/lj7_prod" "Philic_BCP2/200NC/lj7_prod" "Philic_BCP2/50NC/lj7_prod" "Phobic_BCP2/300NC/lj7_prod" "Phobic_BCP2/200NC/lj7_prod" "Phobic_BCP2/50NC/lj7_prod")
lastIndex=5 # six subdirectories to work with, this is the index of last directory if start counting from zero
for iDir in `seq 0 "$lastIndex"`;do
    #mkdir -p $PROJD/secondBatchWrapped/${subDirs[$iDir]}
    #cd $PROJD/secondBatchWrapped/${subDirs[$iDir]}
    for name in poly1 poly2 polyCharge tail head Ncions Pcions;do
	#ln -s $PROJD/${subDirs[$iDir]}/${name}.pdb .
	#CREATE DCD FILES
         #/bin/cp $PROJD/vmd/xyz2dcd.tcl junk.tcl
         # XYZDIR is Monojoy's directory where the xyz file is located
	 #python $PROJD/python/replace_keyword.py junk.tcl "_XYZDIR_" $zf4Root/${subDirs[$iDir]}
  	 #python $PROJD/python/replace_keyword.py junk.tcl "_NAME_" $name
         #vmd -dispdev text -eofexit -e junk.tcl  #produces file ${name}.dcd
	#CREATE UNWRAPED DCD FILES
	 #/bin/cp $PROJD/vmd/xyz2unwrappeddcd.v2.tcl junk.tcl
	 #python $PROJD/python/replace_keyword.py junk.tcl "_XYZDIR_" $zf4Root/${subDirs[$iDir]}
  	 #python $PROJD/python/replace_keyword.py junk.tcl "_NAME_" $name
	 #vmd -dispdev text -eofexit -e junk.tcl
	#REMOVE TRANSLATIONS AND ROTATIONS AT EACH ELEMENTARY STEP
	 #/bin/mv ${name}_unwraped.dcd ${name}_unwrapped.dcd
	 #/bin/cp $PROJD/cpptraj/rms2prev.cpptraj junk.cpptraj
	 #perl -pi -e "s/_NAME_/$name/g" junk.cpptraj
	 #cpptraj -p $name.pdb -i junk.cpptraj &> /dev/null
	#CALCULATE MSD FOR rms2prev TRAJECTORIES
	 # python $PROJD/python/diffusion_t0average.py $name.pdb ${name}_rms2prev.dcd 10000 500 7500 20 '@*' diffusion_${name}_rms2prev.dat --rms2t0=no &> /dev/null
    done
done

##############################
# THIRD BATCH OF SIMULATIONS
##############################
zf4Root="/SNSlocal/projects/zf4/Block_PE_surf"
subDirs=("Philic_BCP2/400NC/lj8_prod" "Philic_BCP2/300NC/lj8_prod" "Philic_BCP2/200NC/lj8_prod" "Philic_BCP2/50NC/lj8_prod" "Phobic_BCP2/300NC/lj8_prod" "Phobic_BCP2/200NC/lj8_prod" "Phobic_BCP2/50NC/lj8_prod")
lastIndex=6 # seven subdirectories to work with, this is the index of last directory if start counting from zero
for iDir in `seq 0 "$lastIndex"`;do
    mkdir -p $PROJD/${subDirs[$iDir]}
    #Use VMD: load dump.xyz and save the first frame as dump.pdb
    cd $PROJD/${subDirs[$iDir]}/
    #grep ' H '  dump.pdb > poly1.pdb       # neutral monomers of PE
    #grep ' He ' dump.pdb > poly2.pdb       # phylic or phobic block of PE
    #grep ' Li ' dump.pdb > polyCharge.pdb  # cation monomers of PE
    #grep ' B '  dump.pdb > tail.pdb        # surfactant tail
    #grep ' Be ' dump.pdb > head.pdb        # surfactant head
    #grep ' C '  dump.pdb > Ncions.pdb      # ? counter-ions
    #grep ' N '  dump.pdb > Pcions.pdb      # ? counter-ions
    #perl -p -i -e 's/ He/ H /g' poly2.pdb       #change element type to hydrogen
    #perl -p -i -e 's/ Li/ H /g' polyCharge.pdb
    #perl -p -i -e 's/ B / H /g' tail.pdb
    #perl -p -i -e 's/ Be/ H /g' head.pdb
    #perl -p -i -e 's/ C / H /g' Ncions.pdb
    #perl -p -i -e 's/ N / H /g' Pcions.pdb
    #perl -p -i -e 's/    H/     /g' poly1.pdb     # remove trailing element name
    #perl -p -i -e 's/   HE/     /g' poly2.pdb
    #perl -p -i -e 's/   LI/     /g' polyCharge.pdb
    #perl -p -i -e 's/   BE/     /g' head.pdb    
    #perl -p -i -e 's/    B/     /g' tail.pdb
    #perl -p -i -e 's/    C/     /g' Ncions.pdb
    #perl -p -i -e 's/    N/     /g' Pcions.pdb
    for name in poly1 poly2 polyCharge tail head Ncions Pcions;do
	#CREATE DCD FILES
	 #/bin/cp $PROJD/vmd/dat2dcd_thirdBatch.tcl junk.tcl
         # XYZDIR is Monojoy's directory where the xyz file is located
	 #python $PROJD/python/replace_keyword.py junk.tcl "_XYZDIR_" $zf4Root/${subDirs[$iDir]}
  	 #python $PROJD/python/replace_keyword.py junk.tcl "_NAME_" $name
	 #vmd -dispdev text -eofexit -e junk.tcl  #generates _NAME_dcd
	 #/bin/rm junk.tcl
	#CREATE UNWRAPED DCD FILES
	 #/bin/cp $PROJD/vmd/dat2unwrapeddcd_thirdBatch.tcl junk.tcl
         # XYZDIR is Monojoy's directory where the xyz file is located
	 #python $PROJD/python/replace_keyword.py junk.tcl "_XYZDIR_" $zf4Root/${subDirs[$iDir]}
  	 #python $PROJD/python/replace_keyword.py junk.tcl "_NAME_" $name
	 #vmd -dispdev text -eofexit -e junk.tcl
	 #/bin/rm junk.tcl
	#REMOVE TRANSLATIONS AND ROTATIONS AT EACH ELEMENTARY STEP (do "module load amber")
	 #/bin/cp $PROJD/cpptraj/rms2prev.cpptraj junk.cpptraj
	 #perl -pi -e "s/_NAME_/$name/g" junk.cpptraj
	 #cpptraj -p $name.pdb -i junk.cpptraj &> /dev/null
	#CALCULATE MSD FOR rms2prev TRAJECTORIES
	 #python $PROJD/python/diffusion_t0average.py $name.pdb ${name}_rms2prev.dcd 20000 500 10000 20 '@*' diffusion_${name}_rms2prev.dat --rms2t0=no
    done
done

#####################################
# THIRD BATCH OF SIMULATIONS, MERGED
#####################################
#MERGE XYZ TRAJECTORIES (note some source directories are in Monojoy's area, and some in mine)
python $PROJD/python/merge_trajectories.py /SNSlocal/projects/zf4/Block_PE_surf/Philic_BCP2/50NC/lj8_prod initial.dat 50 peblock.dump
python $PROJD/python/merge_trajectories.py /SNSlocal/projects/jbq/PEblock/Philic_BCP2/100NC/lj8_prod initial.dat 100 peblock.dump
python $PROJD/python/merge_trajectories.py /SNSlocal/projects/zf4/Block_PE_surf/Philic_BCP2/200NC/lj8_prod initial.dat 200 peblock.dump
python $PROJD/python/merge_trajectories.py /SNSlocal/projects/zf4/Block_PE_surf/Philic_BCP2/300NC/lj8_prod initial.dat 300 peblock.dump
python $PROJD/python/merge_trajectories.py /SNSlocal/projects/zf4/Block_PE_surf/Philic_BCP2/400NC/lj8_prod initial.dat 400 peblock.dump
python $PROJD/python/merge_trajectories.py /SNSlocal/projects/zf4/Block_PE_surf/Phobic_BCP2/50NC/lj8_prod initial.dat 50 peblock.dump
python $PROJD/python/merge_trajectories.py /SNSlocal/projects/jbq/PEblock/Phobic_BCP2/100NC/lj8_prod initial.dat 100 peblock.dump
python $PROJD/python/merge_trajectories.py /SNSlocal/projects/zf4/Block_PE_surf/Phobic_BCP2/200NC/lj8_prod initial.dat 200 peblock.dump
python $PROJD/python/merge_trajectories.py /SNSlocal/projects/zf4/Block_PE_surf/Phobic_BCP2/300NC/lj8_prod initial.dat 300 peblock.dump
python $PROJD/python/merge_trajectories.py /SNSlocal/projects/jbq/PEblock/Phobic_BCP2/400NC/lj8_prod initial.dat 400 peblock.dump

subDirs=("Philic_BCP2/50NC/lj8_prod" "Philic_BCP2/100NC/lj8_prod" "Philic_BCP2/200NC/lj8_prod" "Philic_BCP2/300NC/lj8_prod" "Philic_BCP2/400NC/lj8_prod" "Phobic_BCP2/50NC/lj8_prod" "Phobic_BCP2/100NC/lj8_prod" "Phobic_BCP2/200NC/lj8_prod" "Phobic_BCP2/300NC/lj8_prod" "Phobic_BCP2/400NC/lj8_prod")
nchains=(50 100 200 300 400 50 100 200 300 400)
lastIndex=9 # ten subdirectories to work with, this is the index of last directory if start counting from zero
for iDir in `seq 0 "$lastIndex"`;do
    #mkdir -p $PROJD/${subDirs[$iDir]}/merged
    cd $PROJD/${subDirs[$iDir]}/merged
    #vmd -dispdev text -e $PROJD/vmd/create_psf.tcl  # create PSF topologies from the LAMMPS data file.
    #vmd -dispdev text -e $PROJD/vmd/unwrap_lammpstrj.tcl # unwrap trajectories and save in DCD format
    # CREATE animated GIFs by manually opening vmd and following commands from $PROJD/vmd/create_animated.tcl
    #vmd -dispdev win -e $PROJD/vmd/create_animated.tcl
    #vmd -dispdev text -e $PROJD/vmd/generate_pdbfile.tcl  # CREATE PDB file of the whole system
    #python $PROJD/python/rename_allH.py # create PDB file for creation of center of mass trajectories
    #cpptraj -i $PROJD/amber/center_of_mass_position.cpptraj # Create trajectories for the center of mass
    #python $PROJD/python/diffusion_t0average.py peblock_cm.parm peblock_cm.dcd 20000 500 10000 100 '@2' diffusion_peblock_cm.dat --rms2t0=no
    #python $PROJD/python/rename_atomnames.py # Create PDB file suitable for amber calculations
done
