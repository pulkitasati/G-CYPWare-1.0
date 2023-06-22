#how to use
# ./docking.sh ligand name only 
#( ligand pdb in current directory)

source /home/software/ambermpi4py/amber22/amber.sh
chmod 777 *
dos2unix *
clear
ligand=$1
pdb4amber -i $ligand.pdb -o compatible-$ligand.pdb 
for i in $(ls compatible-* | grep -v compatible-$ligand.pdb);do echo $i removed;rm $i;done
input_file=compatible-$ligand.pdb
output_file=docked-$ligand.pdb

# Read input file line by line
while IFS= read -r line; do
    # Replace "UNK" with "LIG" using sed
    modified_line=$(echo "$line" | sed 's/UNK/LIG/g')
    # Append modified line to the output file
    echo "$modified_line" >> "$output_file"
done < "$input_file"

echo "Replacement completed. Modified content saved in $output_file."
obabel -i pdb docked-$ligand.pdb -o pdbqt -O docked-$ligand.pdbqt
vina --receptor 1bvy_af.pdbqt --config config-file.txt --ligand docked-$ligand.pdbqt --out docked-$ligand-out.pdbqt --log docked-$ligand-out_log.txt
