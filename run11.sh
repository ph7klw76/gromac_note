#!/bin/bash -l
  
#SBATCH --partition=cpu-epyc-genoa
#SBATCH --job-name=remove_boundary
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:59:59
 
# Load dependencies and libraries
source /app/hpcx/2.17.1/hpcx-init.sh
hpcx_load
 
module load gromacs/2023.2
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo -e "0\nq" | gmx_mpi trjconv -f npt3.gro -s npt3.tpr -pbc mol -o output_whole.gro

# File to store input for make_ndx
echo "2" > commands.txt
for i in $(seq 1 1001); do  # Replace N with the total number of residues
    echo "ri $i" >> commands.txt
done
echo "q" >> commands.txt

# Run gmx_mpi make_ndx
gmx_mpi make_ndx -f npt3.gro -o index.ndx < commands.txt

sleep 2

# Input file and necessary GROMACS files
input_file="nearest_neighbor_distances.txt"
gro_file="output_whole.gro"
index_file="index.ndx"

# Extract the first 120 lines from the input file
temp_file="first_120_lines.txt"
head -n 400 "$input_file" > "$temp_file"

# Loop through each selected line of the temporary file
while read -r line; do
    # Extract the first and second column (residue numbers)
    residue1=$(echo $line | awk '{print $1}')
    residue2=$(echo $line | awk '{print $2}')

    # Generate the output filename based on the residues
    output_file1="${residue1}.pdb"
    output_file2="${residue2}.pdb"

    # Run the GROMACS command to select the residues and write the output
    echo -e "r_${residue1}\nq" | gmx_mpi trjconv -f "$gro_file" -s "$gro_file" -n "$index_file" -o "$output_file1"
    echo -e "r_${residue2}\nq" | gmx_mpi trjconv -f "$gro_file" -s "$gro_file" -n "$index_file" -o "$output_file2"

done < "$temp_file"

# Clean up temporary file
rm "$temp_file"

sleep 2

# Input text file containing the file pairs
input_list="nearest_neighbor_distances.txt"

# Output directory for the combined .xyz files
output_dir="output_xyz_files"
mkdir -p "$output_dir"  # Create the directory if it doesn't exist

# Function to extract atom coordinates from a file
process_pdb() {
    local input_file="$1"
    awk '
    $1 == "ATOM" {
        a = substr($3, 1, 1); # Extract first letter of the 3rd column
        x = $6; y = $7; z = $8; # Extract x, y, z (columns 6, 7, 8)
        printf("%s %s %s %s\n", a, x, y, z)
    }' "$input_file"
}

# Process the first 120 lines of the text file
head -n 400 "$input_list" | while read -r line; do
    # Parse the file names from each line
    xx=$(echo "$line" | awk '{print $1}')
    xxx=$(echo "$line" | awk '{print $2}')
    
    file1="${xx}.pdb"
    file2="${xxx}.pdb"
    
    # Define the output file name
    output_file="${output_dir}/${xx}_${xxx}.xyz"

    # Initialize the output file
    > "$output_file"

    # Extract data from both files
    data1=$(process_pdb "$file1")
    data2=$(process_pdb "$file2")

    # Combine the results
    all_lines="$data1"$'\n'"$data2"

    # Calculate the total number of lines extracted
    total_lines=$(echo "$all_lines" | wc -l)

    # Write to the output file
    {
        echo "$total_lines"
        echo "extracted"
        echo "$all_lines"
    } > "$output_file"

    # Print confirmation for each file pair
    echo "Processed: $file1 and $file2 -> $output_file"
done






