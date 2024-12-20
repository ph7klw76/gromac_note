import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals
import numpy as np
import os

# Function to calculate torsional angles for a molecule
def calculate_torsional_angles(filename, atom_indices_list):
    try:
        # Load the .xyz file
        u = mda.Universe(filename, format="XYZ")
        positions = u.atoms.positions

        # Calculate torsional angles
        torsional_angles = []
        for indices in atom_indices_list:
            dihedral_angle = calc_dihedrals(
                positions[indices[0]],
                positions[indices[1]],
                positions[indices[2]],
                positions[indices[3]]
            )
            angle_in_degrees = np.degrees(dihedral_angle)
            if angle_in_degrees<0:
                angle_in_degrees +=180
            torsional_angles.append(angle_in_degrees)

        return torsional_angles
    except Exception as e:
        print(f"Error processing file {filename}: {e}")
        return None

# Function to save torsional angles and statistics to a file
def save_torsional_data(filename, torsional_data):
    with open(filename, "w") as f:
        for xyz_file, angles in torsional_data.items():
            f.write(f"File: {xyz_file}\n")
            for indices, angle in angles.items():
                f.write(f"  Indices {indices}: {angle:.2f} degrees\n")
            f.write("\n")

        # Calculate and save statistics
        all_angles = [angle for angles in torsional_data.values() for angle in angles.values()]
        if all_angles:
            mean_angle = np.mean(all_angles)
            std_dev_angle = np.std(all_angles)
            f.write(f"Mean of all angles: {mean_angle:.2f} degrees\n")
            f.write(f"Standard deviation of all angles: {std_dev_angle:.2f} degrees\n")

# Main program
if __name__ == "__main__":
    # Define the atom indices for torsional angle calculations
    atom_indices_list = [[26,14,15,24],[91,79,80,89]]
    indices_string = "_".join("_".join(map(str, indices)) for indices in atom_indices_list)
    atom_indices_list = [[value - 1 for value in indices] for indices in atom_indices_list]


    # Get all .xyz files in the current directory
    xyz_files = [f for f in os.listdir() if f.endswith(".xyz")]

    # Dictionary to store torsional data
    torsional_data = {}

    # Loop through each .xyz file
    for xyz_file in xyz_files:
        print(f"Processing file: {xyz_file}")
        torsional_angles = calculate_torsional_angles(xyz_file, atom_indices_list)

        if torsional_angles is not None:
            torsional_data[xyz_file] = {
                tuple(indices): angle for indices, angle in zip(atom_indices_list, torsional_angles)
            }
            print(f"Torsional angles for {xyz_file}:")
            for indices, angle in zip(atom_indices_list, torsional_angles):
                print(f"  Indices {indices}: {angle:.2f} degrees")

        print("------------------------------------")

    # Save data to a file
    output_filename = f'{indices_string}.txt'
    save_torsional_data(output_filename, torsional_data)
    print(f"Torsional data saved to {output_filename}")

