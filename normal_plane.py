import MDAnalysis as mda
import numpy as np
import os

# Function to calculate the normal of a plane defined by three points
def calculate_plane_normal(positions, indices):
    p1, p2, p3 = positions[indices[0]], positions[indices[1]], positions[indices[2]]
    v1 = p2 - p1
    v2 = p3 - p1
    normal = np.cross(v1, v2)
    normal /= np.linalg.norm(normal)  # Normalize the vector
    return normal

# Function to calculate angles between a vector and the coordinate axes
def calculate_angles(normal):
    x_axis = np.array([1, 0, 0])
    y_axis = np.array([0, 1, 0])
    z_axis = np.array([0, 0, 1])

    angle_x = np.degrees(np.arccos(np.dot(normal, x_axis) / np.linalg.norm(normal)))
    angle_y = np.degrees(np.arccos(np.dot(normal, y_axis) / np.linalg.norm(normal)))
    angle_z = np.degrees(np.arccos(np.dot(normal, z_axis) / np.linalg.norm(normal)))

    return angle_x, angle_y, angle_z

# Function to calculate plane normals and their angles for a molecule
def calculate_plane_normals_and_angles(filename, atom_indices_list):
    try:
        # Load the .xyz file
        u = mda.Universe(filename, format="XYZ")
        positions = u.atoms.positions

        # Calculate plane normals and angles
        plane_normals = []
        angles = []
        for indices in atom_indices_list:
            # Calculate the normal of the plane
            normal = calculate_plane_normal(positions, indices)
            plane_normals.append(normal)

            # Calculate angles with respect to x, y, z axes
            angle_x, angle_y, angle_z = calculate_angles(normal)
            angles.append((angle_x, angle_y, angle_z))

            # Print positions and normals for debugging
            print(f"Positions of indices {indices}:", positions[indices])
            print(f"Normal: {normal}")
            print(f"Angles: X: {angle_x:.2f}°, Y: {angle_y:.2f}°, Z: {angle_z:.2f}°")

        return plane_normals, angles
    except Exception as e:
        print(f"Error processing file {filename}: {e}")
        return None, None

# Function to save plane normals and angles to a file
def save_plane_normals_and_angles(filename, plane_normal_data, angles_data):
    with open(filename, "w") as f:
        for xyz_file, (normals, angles) in plane_normal_data.items():
            f.write(f"File: {xyz_file}\n")
            for (indices, normal), angle in zip(normals.items(), angles):
                f.write(f"  Indices {indices}: Normal: {normal}\n")
                f.write(f"  Angles with axes: X: {angle[0]:.2f}°, Y: {angle[1]:.2f}°, Z: {angle[2]:.2f}°\n")
            f.write("\n")

        # Calculate and save statistics
        all_angles = [angle for _, angles in plane_normal_data.values() for angle in angles]
        if all_angles:
            mean_angles = np.mean(all_angles, axis=0)
            std_dev_angles = np.std(all_angles, axis=0)
            f.write(f"Mean Angles: X: {mean_angles[0]:.2f}°, Y: {mean_angles[1]:.2f}°, Z: {mean_angles[2]:.2f}°\n")
            f.write(f"Standard Deviation: X: {std_dev_angles[0]:.2f}°, Y: {std_dev_angles[1]:.2f}°, Z: {std_dev_angles[2]:.2f}°\n")

    print(f"Plane normals and angles saved to {filename}")

# Main program
if __name__ == "__main__":
    # Define the atom indices for plane normal calculations
    atom_indices_list = [[2,20,62],[67,85,127]]
    indices_string = "_".join("_".join(map(str, indices)) for indices in atom_indices_list)
    atom_indices_list = [[value - 1 for value in indices] for indices in atom_indices_list]

    # Get all .xyz files in the current directory
    xyz_files = [f for f in os.listdir() if f.endswith(".xyz")]

    # Dictionary to store plane normals and angles
    plane_normal_data = {}

    # Loop through each .xyz file
    for xyz_file in xyz_files:
        print(f"Processing file: {xyz_file}")
        plane_normals, angles = calculate_plane_normals_and_angles(xyz_file, atom_indices_list)

        if plane_normals is not None:
            plane_normal_data[xyz_file] = (
                {
                    tuple(indices): normal.tolist() for indices, normal in zip(atom_indices_list, plane_normals)
                },
                angles
            )
            print(f"Plane normals and angles for {xyz_file}:")
            for indices, normal, angle in zip(atom_indices_list, plane_normals, angles):
                print(f"  Indices {indices}: Normal: {normal}")
                print(f"  Angles: X: {angle[0]:.2f}°, Y: {angle[1]:.2f}°, Z: {angle[2]:.2f}°")

        print("------------------------------------")

    # Save data to a file
    output_filename = f'{indices_string}_normals_and_angles.txt'
    save_plane_normals_and_angles(output_filename, plane_normal_data, angles)
    print(f"Plane normals and angles saved to {output_filename}")


