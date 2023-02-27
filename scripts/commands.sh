#extraction of desinred ligand pose
function ligand_pose_extraction() {
cd $cwd
obabel $ligand.pdbqt -O $dirpath/$ligand.pdb -m
cd $dirpath                            
for i in $(ls $ligand*.pdb | grep -v $posenumber ); do rm $i ; done 
obabel $ligand"$posenumber".pdb -O $ligand"$posenumber"_h.pdb -h            
grep ATOM $ligand"$posenumber"_h.pdb > $ligand.pdb                          
for i in $(ls $ligand* |grep -v $ligand.pdbqt | grep -v $ligand.pdb); do rm $i ; done 
}

#parameterization of extracted ligand pose
function ligand_parameterization() {
source /home/software/ambermpi4py/amber22/amber.sh
antechamber -i $ligand.pdb -fi pdb -o $ligand.mol2 -fo mol2 -c bcc -s 2 -nc $lignad_charge -pf yes 
parmchk2 -i $ligand.mol2 -f mol2 -o $ligand.frcmod
antechamber -i $ligand.mol2 -fi mol2 -o $ligand.ac -fo ac -pf yes
antechamber -i $ligand.ac -fi ac -o $ligand-test.pdb -fo pdb  
}

# 2 positional parametes 1 funtion name 2 file for appending funtion commands
function append_commands() {
echo $(declare -f $1 | grep -v "$1 ()"| grep -v "{" | grep -v "}") >> $2
}

#frcmod generationo of ligand for bound protein
function tleap_ligand() {
ligand_pose_extraction
slurm_sript slurm tleap 
append_commands ligand_parameterization slurm.sh
sbatch slurm.sh
jobid1=$(squeue -u jrfnsm | grep tleap | awk '{print $1}')
echo $jobid1

}


#positioanl parameter 1 directory 
function parameters_files_generation() {
cd $cwd/$dirpath
mkdir $1
cd $1
mkdir $oxddir
mkdir $reddir
cd $oxddir
oxidized_state_files
cd $cwd/$dirpath/$1/$reddir
reduced_state_files
}

# comenting out ligand
function hashing_tleap() {
# 1 positinal parameter 1 directory
function commenting_out_ligand_from_tleap() {
cd $cwd/$dirpath/apo/$1
sed -i '16,17s/^/#/' tleap.in
cp $cwd/$protein.pdb protein.pdb
}

commenting_out_ligand_from_tleap $oxddir
commenting_out_ligand_from_tleap $reddir 
}

#$function 
function waiting_tleap() {
while true
do
  status=$(squeue -u jrfnsm | grep $jobid1 | wc -l)
  if [ $status -eq 0 ]
  then
    echo "Job $jobid1 has finished."
    break
  else
    sleep 30
    echo "Waiting for job $jobid1 to finish..."
  fi
done
cd $cwd/$dirpath
for i in $(ls *.pdb sqm* *.ac | grep -v $protein.pdb|grep -v $ligand.pdb ); do rm $i ; done
cp $cwd/$protein.pdb $cwd/$dirpath/bound/$protein.pdb
cp $cwd/$dirpath/$ligand.pdb $cwd/$dirpath/bound/$ligand.pdb
cd $cwd/$dirpath/bound
cat $protein.pdb $ligand.pdb | awk '$1!="END"' > $protein-$ligand-complex.pdb
pdb4amber -i $protein-$ligand-complex.pdb -o protein.pdb
rm *renum.txt *sslink *nonprot.pdb
cp protein.pdb $oxddir/protein.pdb
cp protein.pdb $reddir/protein.pdb
rm *.pdb

#$function 2 positional parameters 1 file name with extension 2 directory
function copying_files() {
cp $cwd/$dirpath/$1 $2/$1
}

#$function 1 positional parameters 1 directory
function both_dir() {
copying_files $ligand.mol2 $1
copying_files $ligand.frcmod $1
}

both_dir $oxddir
both_dir $reddir

}

#run tleap 2 positional parameters 1 dirctory_bond/apo 2 dir_oxy/reduced
function run_tleap() {
cd $cwd/$dirpath/$1/$2
tleap -f tleap.in
amb2gro_top_gro.py -p protein_solv.prmtop -c protein_solv.inpcrd -t complex.top -g solv_ions.gro -b complex.pdb
mkdir presimulation
for i in $(ls *.pdb *.mol2 *.inpcrd *.prmtop *.frcmod *.log *.in);do mv $i ./presimulation/ ;done
MD_files
}

# 1 positional argument dir bound/apo
function run_tleap_dir() {
run_tleap $1 $oxddir
run_tleap $1 $reddir
}

##apo protein
function apo_state() {
#apo protein
parameters_files_generation apo
hashing_tleap
run_tleap_dir apo
}

## protein - ligand complex
function bound_state() {
tleap_ligand
parameters_files_generation bound
waiting_tleap
run_tleap_dir bound
}

#function index_complex 1 directory positional parmameters (oxy/reduced)
function index_complex
{
cd $cwd/$dirpath/bound/$1
printf "1|13|14|15|16|17\nq\n"|gmx -quiet make_ndx -f solv_ions.gro -o index.ndx
}

#function 1 directory positional parmameters (oxy/reduced)
function index_apo
{
cd $cwd/$dirpath/apo/$1
printf "1|13|14|15|16\nq\n"|gmx -quiet make_ndx -f solv_ions.gro -o index.ndx
}

## complex indexging
function indexing_bound() {
index_complex $oxddir
index_complex $reddir
}

## complex indexging
function indexing_apo() {
index_apo $oxddir
index_apo $reddir
}

## final md running
function run_MD() {
gmx_mpi -quiet grompp -f em.mdp -c solv_ions.gro -n index.ndx -p complex.top -po em_mdout.mdp -pp em_processed.top -o em.tpr -maxwarn 3 
gmx_mpi -quiet mdrun -s em.tpr -mp em_processed.top -mn index.ndx -o em_traj.trr -x em_traj_comp.xtc -cpo em_state.cpt -c em_confout.gro -e em_ener.edr  -g em_md.log -xvg xmgrace -nb gpu &>> all_output.txt
echo "11 0" | gmx_mpi -quiet energy -f em_ener.edr -s em.tpr -o potential_energy.xvg -xvg xmgrace -fee -dp -mutot &>> all_output.txt
gmx_mpi -quiet grompp -f nvt.mdp -c em_confout.gro -n index.ndx -p em_processed.top -po nvt_mdout.mdp -pp nvt_processed.top -o nvt.tpr -maxwarn 3
gmx_mpi -quiet mdrun -s nvt.tpr -mp nvt_processed.top -mn index.ndx -o nvt_traj.trr -x nvt_traj_comp.xtc -cpo nvt_state.cpt -c nvt_confout.gro -e nvt_ener.edr  -g nvt_md.log -xvg xmgrace -nb gpu &>> all_output.txt
echo "15 0" | gmx_mpi -quiet energy -f nvt_ener.edr -s nvt.tpr -o temperature.xvg -xvg xmgrace -fee -dp -mutot &>> all_output.txt
gmx_mpi -quiet grompp -f npt.mdp -c nvt_confout.gro -n index.ndx -p nvt_processed.top -t nvt_state.cpt -e nvt_ener.edr -po npt_mdout.mdp -pp npt_processed.top -o npt.tpr -maxwarn 3
gmx_mpi -quiet mdrun -s npt.tpr -mp npt_processed.top -mn index.ndx -o npt_traj.trr -x npt_traj_comp.xtc -cpo npt_state.cpt -c npt_confout.gro -e npt_ener.edr  -g npt_md.log -xvg xmgrace -nb gpu &>> all_output.txt 
echo "16 0" | gmx_mpi -quiet energy -f npt_ener.edr -s npt.tpr -o pressure.xvg -xvg xmgrace -fee -dp -mutot &>> all_output.txt
echo "22 0" | gmx_mpi -quiet energy -f npt_ener.edr -s npt.tpr -o density.xvg -xvg xmgrace -fee -dp -mutot &>> all_output.txt
gmx_mpi -quiet grompp -f md.mdp -c npt_confout.gro -n index.ndx -p npt_processed.top -t npt_state.cpt -po md_mdout.mdp -pp md_processed.top -o md.tpr -maxwarn 3
gmx_mpi -quiet mdrun -s md.tpr -mp md_processed.top -mn index.ndx -o md_traj.trr -x md_traj_comp.xtc -cpo md_state.cpt -c md_confout.gro -e md_ener.edr  -g md_md.log -xvg xmgrace -nb gpu &>> all_output.txt
gmx_mpi -quiet check -f md_traj_comp.xtc -s1 md.tpr -c md_confout.gro -e md_ener.edr -n index.ndx -m doc.tex &>> check.txt
echo "1 0" |gmx_mpi -quiet trjconv -f md_traj_comp.xtc -s md.tpr -n index.ndx -o md_traj_center.xtc -xvg xmgrace -center -pbc mol -ur compact -force yes -conect yes &>> all_output.txt
echo "0" | gmx_mpi -quiet trjconv -s md.tpr -f md_traj_center.xtc -o start.pdb -dump 0 &>> all_output.txt
echo "4 0" |gmx_mpi -quiet trjconv -s md.tpr -f md_traj_center.xtc -o md_fit.xtc -fit rot+trans &>> all_output.txt
gmx_mpi -quiet grompp -f ie.mdp -c md_confout.gro -n index.ndx -p md_processed.top -t md_state.cpt -po ie_mdout.mdp -pp ie_processed.top -o ie.tpr -maxwarn 3
gmx_mpi -quiet mdrun -s ie.tpr -mp ie_processed.top -mn index.ndx -rerun md_traj_comp.xtc -o ie_traj.trr -x ie_traj_comp.xtc -cpo ie_state.cpt -c ie_confout.gro -e ie_ener.edr  -g ie_md.log -xvg xmgrace -nb gpu &>> all_output.txt
echo "11 0" | gmx_mpi -quiet energy -f ie_ener.edr -s ie.tpr -o single_point_energy.xvg -xvg xmgrace -fee -dp -mutot &>> spe.txt
gmx_mpi -quiet distance -f md_traj_comp.xtc -s ie.tpr -oav distave.xvg -oall dist.xvg -oxyz distxyz.xvg -oh disthist.xvg -oallstat diststat.xvg -xvg xmgrace -select 'resname "FE1"  plus resname "FMN" and name P ' &>> distance.txt
mkdir postsimulation_files
for i in $(ls | grep -v all_output.txt | grep -v md_fit.xtc |grep -v md_traj_center.xtc |grep -v md_traj_comp.xtc| grep -v postsimulation_files |grep -v presimulation);do mv $i ./postsimulation_files/ ;done

}


## all dir  md running 1 positional parameter single letter job name
function all_dir_md() {
slurm_sript $1-md $1-ms
append_commands run_MD $1-md.sh
sbatch $1-md.sh
jobid1=$(squeue -u jrfnsm | grep $1-ms | awk '{print $1}')
echo $jobid1
}

## final md running 3 parameters 1 apo or bound 2 oxi or reduced dir 3 single or 2 letter job name for md
function fs_md() {
cd $cwd/$dirpath/$1/$2
all_dir_md $3
}

## apo all dir md 
function apo_all_md_submit() {
fs_md apo $oxddir ao
oxyapojobid=$jobid1
echo "oxyapojobid: $oxyapojobid"
fs_md apo $reddir ar
reducedapojobid=$jobid1
echo "reducedapojobid: $reducedapojobid"
}

## bound all dir md 
function bound_all_md_submit() {
fs_md bound $oxddir bo
oxyboundjobid=$jobid1
echo "oxyboundjobid: $oxyboundjobid"
fs_md bound $reddir br
reducedboundjobid=$jobid1
echo "reducedboundjobid: $reducedboundjobid"
}



function copy_files_spec() {
##1 apo/bound 2 oxy/red 3 filename with extension
cp $cwd/$dirpath/$1/$2/postsimulation_files/$3 copied-$3
}


## 1 oxy/red 2 red/oxy
function apo_copy_files_spec() {
cd  $cwd/$dirpath/apo/$1/
copy_files_spec apo $2 ie_ener.edr
copy_files_spec apo $2 ie_processed.top
copy_files_spec apo $2 index.ndx
copy_files_spec apo $2 ie.tpr

}

## no argument
function apo_copy_files_spec_done() {
apo_copy_files_spec $oxddir $reddir
apo_copy_files_spec $reddir $oxddir 
}

## 1 oxy/red 2 red/oxy
function bound_copy_files_spec() {
cd  $cwd/$dirpath/bound/$1/
copy_files_spec bound $2 ie_ener.edr
copy_files_spec bound $2 ie_processed.top
copy_files_spec bound $2 index.ndx
copy_files_spec bound $2 ie.tpr
}

## no argument
function bound_copy_files_spec_done() {
bound_copy_files_spec $oxddir $reddir
bound_copy_files_spec $reddir $oxddir 
}

## 2 positional argument 1 bound/apo 2 oxy/reduced directory
function recalculating_spe() {
cd $cwd/$dirpath/$1/$2
gmx_mpi -quiet mdrun -s copied-ie.tpr -mp copied-ie_processed.top -mn copied-index.ndx -rerun md_traj_comp.xtc -o other_state-ie_traj.trr -x other_state-ie_traj_comp.xtc -cpo other_state-ie_state.cpt -c other_state-ie_confout.gro -e other_state-ie_ener.edr  -g other_state-ie_md.log -xvg xmgrace -nb gpu &>> all_output.txt
echo "11 0" | gmx_mpi -quiet energy -f other_state-ie_ener.edr -s copied-ie.tpr -o other_state-single_point_energy.xvg -xvg xmgrace -fee -dp -mutot &>> os-spe.txt
v1=$(grep "Energy                      Average   Err.Est.       RMSD  Tot-Drift  -kT ln<e^(E/kT)>" os-spe.txt -A 2 | grep Potential | awk '{print $2}')
v2=$(grep "Energy                      Average   Err.Est.       RMSD  Tot-Drift  -kT ln<e^(E/kT)>" ./postsimulation_files/spe.txt -A 2 | grep Potential | awk '{print $2}')
echo "v1 : $v1" &>> final_spe_file.txt
echo "v2 : $v2" &>> final_spe_file.txt
abs_v1=${v1#-}
abs_v2=${v2#-}
echo "abs_v1 : $abs_v1" &>> final_spe_file.txt
echo "abs_v2 : $abs_v2" &>> final_spe_file.txt
sp_energy=$(echo "($abs_v1) - ($abs_v2)" | bc)
abs_sp_energy=${sp_energy#-}
echo "sp_energy : $sp_energy" &>> final_spe_file.txt
echo "abs_sp_energy : $abs_sp_energy" &>> final_spe_file.txt
distance=$(grep "Average distance:" ./postsimulation_files/distance.txt | awk '{print $3}')
a_distance=$(expr $distance*10 | bc)
echo "distance : $a_distance" &>> final_spe_file.txt
}

## apo_distance_and_energy no positional argument 
function apo_distance_and_energy() {
recalculating_spe apo $oxddir 
recalculating_spe apo $reddir
}

## bound_distance_and_energy no positional argument 
function bound_distance_and_energy() {
recalculating_spe bound $oxddir 
recalculating_spe bound $reddir
}

#calculte parmmetrs 1 positional parameter apo or bound
function final_calculation() {
##1 positional argument apo or bound
function fun_grep_energy_all() {
Ea=$(grep "abs_sp_energy :" $cwd/$dirpath/$1/$reddir/final_spe_file.txt | awk '{print $3}')
Eb=$(grep "abs_sp_energy :" $cwd/$dirpath/$1/$oxddir/final_spe_file.txt | awk '{print $3}')
D1=$(grep "distance :" $cwd/$dirpath/$1/$reddir/final_spe_file.txt | awk '{print $3}')
D2=$(grep "distance :" $cwd/$dirpath/$1/$oxddir/final_spe_file.txt | awk '{print $3}')
}
fun_grep_energy_all $1

# calculate lambda and DeltaG
lambda=$(echo "scale=2; ($Ea - $Eb) / 2" | bc)
DeltaG=$(echo "scale=2; ($Ea + $Eb) / 2" | bc)

# Calculate the average of D1 and D2
V=$(echo "scale=5; ($D1+$D2)/2" | bc)

# Print the result
echo "Ea : $Ea" &>> $cwd/$dirpath/$1/$1-FCETP.txt
echo "Eb : $Eb" &>> $cwd/$dirpath/$1/$1-FCETP.txt
echo "D1 : $D1" &>> $cwd/$dirpath/$1/$1-FCETP.txt
echo "D2 : $D2" &>> $cwd/$dirpath/$1/$1-FCETP.txt
echo "V : $V" &>> $cwd/$dirpath/$1/$1-FCETP.txt
echo "lambda : $lambda" &>> $cwd/$dirpath/$1/$1-FCETP.txt
echo "DeltaG : $DeltaG" &>> $cwd/$dirpath/$1/$1-FCETP.txt
function kjol2ev() {
result=$(bc <<< "scale=10; $1 * 0.0103636")
echo $result
}
ev_lambda=$(kjol2ev $lambda)
echo "ev_lambda : $ev_lambda" &>> $cwd/$dirpath/$1/$1-FCETP.txt
ev_deltag=$(kjol2ev $DeltaG)
echo "ev_deltag : $ev_deltag" &>> $cwd/$dirpath/$1/$1-FCETP.txt
lambda=$(echo "scale=8; $ev_lambda/1.6" | bc)
echo "lambda : $lambda" &>> $cwd/$dirpath/$1/$1-FCETP.txt
DeltaG=$(echo $ev_deltag)
echo "DeltaG : $DeltaG" &>> $cwd/$dirpath/$1/$1-FCETP.txt
python3 $cwd/ET_Cal.py $lambda $DeltaG $V $1
echo "job ended"
}

#all command for putting into the ETCal.sh script no positonal argument 
function append_last_step() {
bash postmd_steps.sh GROMACS_CYPWareData.txt
}
