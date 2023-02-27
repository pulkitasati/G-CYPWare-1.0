export $(xargs <$1)
cwd=$(pwd)
 
##sourcing
source /home/software/ambermpi4py/amber22/amber.sh
module load gromacs-cpu/2022.4
source ./scripts/welcome.sh
source ./scripts/files.sh
source ./scripts/commands.sh

## welcome note
welcome_part1

#run the job
cat $1
mkdir $dirpath
### two form apo and bound 

## complex system all steps
function complex_type() {
bound_state
indexing_bound
bound_all_md_submit
cd $cwd 
slurm_sript $protein-bound_ETcal ETC
append_commands append_last_step $protein-bound_ETcal.sh
sbatch --dependency=afterok:$oxyboundjobid:$reducedboundjobid $protein-bound_ETcal.sh
}


## apo system all steps
function apo_type() {
apo_state
indexing_apo
apo_all_md_submit
cd $cwd 
slurm_sript $protein-apo_ETcal ETC
append_commands append_last_step $protein-apo_ETcal.sh
sbatch --dependency=afterok:$oxyapojobid:$reducedapojobid $protein-apo_ETcal.sh
} 




# Check if the variable is set
if [ -z "$calculation_type" ]; then
    echo "No input value given for calcuation type in input file"
    exit 1
fi

# Perform different actions based on the value of the variable
case $calculation_type in
    apo)
        echo "Performing APO calculation..."
        apo_type
        
        ;;
    bound)
        echo "Performing bound calculation..."
        complex_type
        
        ;;
    both)
        echo "Performing both APO and bound calculations..."
        echo "Performing bound calculation..."
        complex_type
        
        echo "Performing APO calculation..."
        apo_type
        
        ;;
    *)
        echo "Unknown calculation type: $calculation_type"
        exit 1
        ;;
esac