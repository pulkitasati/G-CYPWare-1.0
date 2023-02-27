export $(xargs <$1)
cwd=$(pwd)
##sourcing
source /home/software/ambermpi4py/amber22/amber.sh
module load gromacs-gpu/2022.4
source ./scripts/welcome.sh
source ./scripts/files.sh
source ./scripts/commands.sh
##welcome note part 2
welcome_part2
#run the job
cat $1
# Check if the variable is set
if [ -z "$calculation_type" ]; then
    echo "No input value given for calcuation type in input file"
    exit 1
fi

# Perform different actions based on the value of the variable
case $calculation_type in
    apo)
        echo "Performing APO ET calculation..."
        apo_copy_files_spec_done
        apo_distance_and_energy
        final_calculation apo
        
        ;;
    bound)
        echo "Performing bound ET calculation..."
        bound_copy_files_spec_done
        bound_distance_and_energy
        final_calculation bound
        
        ;;
    both)
        echo "Performing both APO and bound ET calculations..."
        echo "Performing bound ET calculation..."
        bound_copy_files_spec_done
        bound_distance_and_energy
        final_calculation bound
        
        echo "Performing APO ET calculation..."
        apo_copy_files_spec_done
        apo_distance_and_energy
        final_calculation apo
        ;;
    *)
        echo "Unknown calculation type: $calculation_type"
        exit 1
        ;;
esac