# Automated Extraction of Molecular Dynamics Results to Calculation of electronic coupling. 

```python
import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import seaborn as sns

# Function to read and parse a .gro file containing molecular trajectory data
def read_gro_file_with_masses(file_path):
    """
    Reads a .gro file and extracts molecular and atomic information.
    
    Args:
        file_path (str): Path to the .gro file.

    Returns:
        dict: A dictionary where keys are molecule numbers, and values are lists of tuples
              containing atomic coordinates (x, y, z) and atom type.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()[2:]  # Skip the first two header lines

    molecules = {}
    for line in lines:
        if not line.strip():
            continue  # Skip empty lines
        try:
            molecule_number = int(line[0:5])  # Extract molecule number from the first 5 characters
            atom_type = line[10:15].strip()  # Extract atom type
            # Extract x, y, z coordinates, converting them to floats
            x, y, z = map(float, [line[20:28].strip(), line[28:36].strip(), line[36:44].strip()])
            # Append data to the molecule's list or create a new entry
            if molecule_number in molecules:
                molecules[molecule_number].append((x, y, z, atom_type))
            else:
                molecules[molecule_number] = [(x, y, z, atom_type)]
        except ValueError:
            continue  # Skip lines that cannot be parsed
    return molecules

# Function to filter specific carbon atoms based on labels
def extract_specific_carbon_atoms(atoms, carbon_labels):
    """
    Filters atoms based on a list of specific carbon labels.
    
    Args:
        atoms (list): List of atom tuples (x, y, z, atom_type).
        carbon_labels (list): List of carbon labels to filter (e.g., ['C27', 'C14']).

    Returns:
        list: Filtered list of atom tuples matching the carbon labels.
    """
    return [atom for atom in atoms if any(c_label in atom[3] for c_label in carbon_labels)]

# Function to calculate nearest-neighbor distances using KD-Tree
def calculate_kd_tree_nearest_specific_carbon_distances(molecules, carbon_labels):
    """
    Computes the nearest-neighbor distances between specific carbon atoms across molecules.
    
    Args:
        molecules (dict): Dictionary of molecules from `read_gro_file_with_masses`.
        carbon_labels (list): List of specific carbon labels to analyze.

    Returns:
        list: List of nearest distances between specific carbon atoms in different molecules.
    """
    # Extract specific carbon atoms for each molecule
    carbon_atoms = {mol_id: extract_specific_carbon_atoms(atoms, carbon_labels) for mol_id, atoms in molecules.items()}

    # Flatten the data into a list of (molecule ID, atom tuple)
    flat_list = [(mol_id, atom) for mol_id, atoms in carbon_atoms.items() for atom in atoms]
    # Extract atomic coordinates as points
    points = np.array([atom[1][:3] for atom in flat_list])
    # Extract molecule IDs for each point
    mol_ids = np.array([atom[0] for atom in flat_list])

    # Build a KD-Tree for efficient nearest-neighbor search
    tree = cKDTree(points)
    distances, indexes = tree.query(points, k=2)  # Query for two nearest neighbors (self + nearest)

    # Filter and store distances where neighbors are in different molecules
    nearest_distances = []
    for dist, idx, mol_id in zip(distances[:, 1], indexes[:, 1], mol_ids):
        if mol_id != mol_ids[idx]:  # Ensure the neighbor belongs to a different molecule
            nearest_distances.append(dist)

    return nearest_distances

# Function to plot the probability density of distances
def plot_probability_density(distances, title):
    """
    Plots the probability density of distances using a KDE (Kernel Density Estimate).
    
    Args:
        distances (list): List of distances to visualize.
        title (str): Title of the plot.
    """
    plt.figure(figsize=(10, 6))  # Set the figure size
    sns.kdeplot(distances)  # Plot KDE
    plt.title(title)  # Add title
    plt.xlabel('Distance (nm)')  # Label x-axis
    plt.ylabel('Probability Density')  # Label y-axis
    plt.show()  # Display the plot

# Main execution
file_path = 'npt3.gro'  # Path to the .gro file
# User Defined Molecule positions
specific_carbon_labels = ['C27', 'C14', 'C12', 'C2', 'C8', 'C5', 'C15', 'C18', 
                          'C22', 'C29', 'C32', 'C35', 'C40', 'C56', 'C44', 'C53', 'C58', 'C62'] # BEWARE the GauView Number X to PDB numnber
```
## Exhaustive Method [ONLY GET THE NEAREST CONTACT OF ONE MOLECULE TO ANOTHER]

```python
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

# Function to read and parse a .gro file containing molecular trajectory data
def read_gro_file_with_masses(file_path):
    """
    Reads a .gro file and extracts molecular and atomic information.
    
    Args:
        file_path (str): Path to the .gro file.

    Returns:
        dict: A dictionary where keys are molecule numbers, and values are lists of tuples
              containing atomic coordinates (x, y, z) and atom type.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()[2:]  # Skip the first two header lines

    molecules = {}
    for line in lines:
        if not line.strip():
            continue  # Skip empty lines
        try:
            molecule_number = int(line[0:5])  # Extract molecule number from the first 5 characters
            atom_type = line[10:15].strip()  # Extract atom type
            x, y, z = map(float, [line[20:28].strip(), line[28:36].strip(), line[36:44].strip()])  # Extract x, y, z
            if molecule_number in molecules:
                molecules[molecule_number].append((x, y, z, atom_type))
            else:
                molecules[molecule_number] = [(x, y, z, atom_type)]
        except ValueError:
            continue  # Skip lines that cannot be parsed
    return molecules

# Function to filter specific carbon atoms based on labels
def extract_specific_carbon_atoms(atoms, carbon_labels):
    """
    Filters atoms based on a list of specific carbon labels.
    
    Args:
        atoms (list): List of atom tuples (x, y, z, atom_type).
        carbon_labels (list): List of carbon labels to filter (e.g., ['C27', 'C14']).

    Returns:
        list: Filtered list of atom tuples matching the carbon labels.
    """
    return [atom for atom in atoms if any(c_label in atom[3] for c_label in carbon_labels)]

# Function to compute unique nearest-neighbor distances between all specific atoms in different molecules
def calculate_exhaustive_nearest_distances(molecules, carbon_labels):
    """
    Computes unique nearest-neighbor distances between specific atoms in different molecules.
    
    Args:
        molecules (dict): Dictionary of molecules with atomic data.
        carbon_labels (list): List of specific carbon labels to include in the calculation.

    Returns:
        list: List of unique pairs (mol_id1, mol_id2, distance).
    """
    nearest_pairs = []
    molecule_ids = list(molecules.keys())

    for i, mol_id1 in tqdm(enumerate(molecule_ids), total=len(molecule_ids), desc="Calculating distances"):
        atoms1 = extract_specific_carbon_atoms(molecules[mol_id1], carbon_labels)
        for j, mol_id2 in enumerate(molecule_ids):
            if i >= j:
                continue  # Avoid self-pairs and duplicate pairs
            atoms2 = extract_specific_carbon_atoms(molecules[mol_id2], carbon_labels)
            
            # Compute distances between all pairs of specific atoms
            for atom1 in atoms1:
                for atom2 in atoms2:
                    dist = np.linalg.norm(np.array(atom1[:3]) - np.array(atom2[:3]))
                    if dist < 0.5:  # Only include distances less than 0.5
                        nearest_pairs.append((mol_id1, mol_id2, dist))

    # Sort pairs by distance
    nearest_pairs.sort(key=lambda x: x[2])

    # Ensure each molecule appears only once in the final list with the shortest distance
    final_pairs = []
    seen_molecules = set()
    for pair in nearest_pairs:
        mol_id1, mol_id2, dist = pair
        if mol_id1 not in seen_molecules and mol_id2 not in seen_molecules:
            final_pairs.append(pair)
            seen_molecules.add(mol_id1)
            seen_molecules.add(mol_id2)

    return final_pairs

# Function to save nearest-neighbor pairs to a file
def save_unique_nearest_neighbor_pairs(file_name, pairs):
    """
    Save unique molecular pairs with their distances to a file.
    
    Args:
        file_name (str): Path to save the file.
        pairs (list): List of unique pairs with distances.
    """
    with open(file_name, 'w') as file:
        for pair in pairs:
            file.write(f"{pair[0]} {pair[1]} {pair[2]:.6f}\n")

# Function to compute average nearest-neighbor distance
def compute_average_distance(pairs):
    """
    Compute the average distance from the list of unique nearest-neighbor pairs.
    
    Args:
        pairs (list): List of unique pairs with distances.

    Returns:
        float: The average distance.
    """
    distances = [pair[2] for pair in pairs]
    if distances:
        return np.mean(distances)
    return None

# Function to plot the probability density distribution of distances
def plot_probability_density(distances, title):
    """
    Plots the probability density distribution of distances using a KDE (Kernel Density Estimate).
    
    Args:
        distances (list): List of distances to visualize.
        title (str): Title of the plot.
    """
    plt.figure(figsize=(10, 6))
    sns.kdeplot(distances, fill=True)
    plt.title(title)
    plt.xlabel("Distance (nm)")
    plt.ylabel("Probability Density")
    plt.show()

# Main program
if __name__ == "__main__":
    file_path = 'npt3.gro'  # Path to the .gro file
    specific_carbon_labels = ['C60','C55','C53','C40','C38','C41','C46','C33','C28','C19','C24','C25','C3'] # BEWARE the GauView Number X to PDB numnber

    try:
        # Step 1: Read and parse the molecular data
        molecules = read_gro_file_with_masses(file_path)

        # Step 2: Calculate unique nearest-neighbor distances between specific atoms in different molecules
        unique_nearest_pairs = calculate_exhaustive_nearest_distances(molecules, specific_carbon_labels)

        # Step 3: Save unique nearest-neighbor pairs to a file
        output_file_unique = 'nearest_neighbor_distances.txt'
        save_unique_nearest_neighbor_pairs(output_file_unique, unique_nearest_pairs)
        print(f"Nearest neighbor pairs with distances < 0.5 saved to {output_file_unique}")

        # Step 4: Compute the average nearest-neighbor distance
        average_distance = compute_average_distance(unique_nearest_pairs)

        if average_distance is not None:
            print(f"Average nearest-neighbor distance: {average_distance:.6f} nm")
        else:
            print("No valid distances found between molecular centroids.")

        # Step 5: Plot the probability density of distances
        distances = [pair[2] for pair in unique_nearest_pairs]
        plot_probability_density(distances, "Probability Density of Nearest-Neighbor Distances (<0.5 nm)")

    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found. Please ensure the file exists.")
```

# Step 1: Read and parse the molecular data
molecules = read_gro_file_with_masses(file_path)

# Step 2: Calculate nearest-neighbor distances for specific carbon atoms
kd_tree_distances = calculate_kd_tree_nearest_specific_carbon_distances(molecules, specific_carbon_labels)

# Step 3: Visualize the distribution of distances
plot_probability_density(kd_tree_distances, 'Probability Density of KD-Tree Nearest Specific Carbon Distances Between Molecules')


This blog delves into the computational pipeline for analyzing molecular spatial data using Python libraries such as `numpy`, `scipy`, `matplotlib`, and `seaborn`. The script processes molecular trajectory data stored in a `.gro` file to:

1. Identify specific carbon atoms.
2. Compute nearest-neighbor distances using a KD-tree algorithm.
3. Visualize the probability density of these distances to provide insights into the spatial organization of selected molecular components.

---

## 1. Input Data Parsing: Reading the `.gro` File

The `.gro` file format is commonly used in molecular dynamics simulations and contains atomic coordinates. The function `read_gro_file_with_masses` is responsible for parsing the file and extracting atomic positions.

### Key Features:
- **Input**: The function reads the `.gro` file line-by-line, skipping the first two header lines.
- **Data Structure**: Atom data is grouped by molecule number. Each molecule's coordinates and atom types are stored in a dictionary:

```python
molecules[molecule_number] = [(x, y, z, atom_type), ...]
```

# Robustness and Workflow for Carbon Atom Spatial Analysis in Molecular Simulations

## Robustness

- **Error Handling**: 
  - Empty lines and non-numeric values are handled gracefully with error-catching mechanisms.

### Example Output:
A dictionary where:
- **Keys**: Molecule IDs.
- **Values**: Lists of tuples representing atom positions and types.

---

## 1. Filtering Carbon Atoms

The `extract_specific_carbon_atoms` function identifies atoms of interest based on their labels. These labels, such as `C27` or `C14`, are passed as input.

### Implementation:
- **Pattern Matching**: Checks if any of the provided carbon labels (`carbon_labels`) are substrings of the atom type.
- **Output**: A filtered list of carbon atoms for each molecule.

This approach allows flexibility in targeting atoms of interest without modifying the core file-reading functionality.

---

## 2. Spatial Analysis Using KD-Tree

The core of this script involves efficient spatial queries with the KD-tree algorithm, implemented through the `scipy.spatial.cKDTree` class. The function `calculate_kd_tree_nearest_specific_carbon_distances` calculates the shortest distances between specific carbon atoms in different molecules.

### Steps:

#### 1. **Flattening Data**
- Extract the positions of specific carbon atoms into a flat array (`points`).
- Track corresponding molecule IDs for associating distances back to molecules.

#### 2. **KD-Tree Construction**
- Build the KD-tree using the `points` array, enabling $O(\log n)$ query times for spatial searches.

#### 3. **Querying Neighbors**
- For each point, the two nearest neighbors are identified (`k=2`), where the first is the point itself.
- The distance to the nearest neighbor in a different molecule is stored.

#### 4. **Edge Cases**
- Ensures distances are computed **only** between atoms in different molecules using molecule ID comparisons.

---

## 3. Visualization with KDE Plot

The `plot_probability_density` function uses a **Kernel Density Estimate (KDE)** to plot the distribution of nearest-neighbor distances. This is achieved through `seaborn.kdeplot`, providing a smooth representation of data density.

### Key Aspects:
- **Probability Density**: Highlights the most probable spatial separations between the selected carbon atoms.
- **Customizability**: Title, axis labels, and plot size are adjustable to suit various analysis contexts.

---

## Workflow Summary

1. **Input Parsing**:  
   The `.gro` file is parsed into a structured format containing molecules and atoms.

2. **Carbon Atom Filtering**:  
   Specific carbon atoms are identified using flexible label matching.

3. **Distance Calculation**:
   - A **KD-tree** is used for efficient nearest-neighbor queries.
   - Distances between specific carbon atoms in different molecules are computed.

4. **Visualization**:  
   The distances are visualized as a **probability density function**, offering insights into spatial organization.

---

## Application Example

This script can be applied in molecular dynamics studies to analyze:

### 1. **Intermolecular Interactions**
- Identifies spatial patterns in specific atom types.
- Helps infer interaction tendencies or packing arrangements.

### 2. **Structural Organization**
- Understanding how specific atoms are spatially distributed.
- Provides insights into material properties or biological function.

---

## Example Output

If the input `.gro` file contains trajectory data for molecules with targeted carbon atoms, the final KDE plot might reveal:
- **Peaks**: Corresponding to preferred intermolecular distances.
- **Tails**: Indicating occasional long-range separations.

You can save the nearest pair by using the code below:
```python
import numpy as np
from scipy.spatial import cKDTree

# Function to read and parse a .gro file containing molecular trajectory data
def read_gro_file_with_masses(file_path):
    """
    Reads a .gro file and extracts molecular and atomic information.
    
    Args:
        file_path (str): Path to the .gro file.

    Returns:
        dict: A dictionary where keys are molecule numbers, and values are lists of tuples
              containing atomic coordinates (x, y, z) and atom type.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()[2:]  # Skip the first two header lines

    molecules = {}
    for line in lines:
        if not line.strip():
            continue  # Skip empty lines
        try:
            molecule_number = int(line[0:5])  # Extract molecule number from the first 5 characters
            atom_type = line[10:15].strip()  # Extract atom type
            x, y, z = map(float, [line[20:28].strip(), line[28:36].strip(), line[36:44].strip()])  # Extract x, y, z
            if molecule_number in molecules:
                molecules[molecule_number].append((x, y, z, atom_type))
            else:
                molecules[molecule_number] = [(x, y, z, atom_type)]
        except ValueError:
            continue  # Skip lines that cannot be parsed
    return molecules

# Function to filter specific carbon atoms based on labels
def extract_specific_carbon_atoms(atoms, carbon_labels):
    """
    Filters atoms based on a list of specific carbon labels.
    
    Args:
        atoms (list): List of atom tuples (x, y, z, atom_type).
        carbon_labels (list): List of carbon labels to filter (e.g., ['C27', 'C14']).

    Returns:
        list: Filtered list of atom tuples matching the carbon labels.
    """
    return [atom for atom in atoms if any(c_label in atom[3] for c_label in carbon_labels)]

# Function to calculate unique nearest-neighbor distances and pairs
def calculate_kd_tree_nearest_specific_carbon_pairs_unique(molecules, carbon_labels):
    """
    Computes unique nearest-neighbor distances between specific carbon atoms across molecules.
    
    Ensures no duplicate pairs (e.g., (1, 2) and (2, 1)) are recorded.
    """
    carbon_atoms = {mol_id: extract_specific_carbon_atoms(atoms, carbon_labels) for mol_id, atoms in molecules.items()}
    flat_list = [(mol_id, atom) for mol_id, atoms in carbon_atoms.items() for atom in atoms]
    points = np.array([atom[1][:3] for atom in flat_list])
    mol_ids = np.array([atom[0] for atom in flat_list])

    tree = cKDTree(points)
    distances, indexes = tree.query(points, k=2)

    nearest_pairs = set()  # Use a set to ensure unique pairs
    for dist, idx, mol_id, point in zip(distances[:, 1], indexes[:, 1], mol_ids, points):
        neighbor_mol_id = mol_ids[idx]
        if mol_id != neighbor_mol_id:  # Ensure the neighbor belongs to a different molecule
            # Add tuple with sorted molecule IDs to avoid duplicates
            pair = tuple(sorted((mol_id, neighbor_mol_id)))
            nearest_pairs.add((pair[0], pair[1], dist))
    return list(nearest_pairs)

# Function to save nearest-neighbor pairs to a file
def save_unique_nearest_neighbor_pairs(file_name, pairs):
    """
    Save unique molecular pairs with their distances to a file.
    
    Args:
        file_name (str): Path to save the file.
        pairs (list): List of unique pairs with distances.
    """
    with open(file_name, 'w') as file:
        for pair in pairs:
            file.write(f"{pair[0]} {pair[1]} {pair[2]:.6f}\n")

# Main program
if __name__ == "__main__":
    file_path = 'npt3.gro'  # Path to the .gro file
    specific_carbon_labels = ['C27', 'C14', 'C12', 'C2', 'C8', 'C5', 'C15', 'C18',
                              'C22', 'C29', 'C32', 'C35', 'C40', 'C56', 'C44', 'C53', 'C58', 'C62'] # BEWARE the GauView Number X to PDB numnber

    # Step 1: Read and parse the molecular data
    try:
        molecules = read_gro_file_with_masses(file_path)

        # Step 2: Calculate unique nearest-neighbor molecular pairs for specific carbon atoms
        unique_nearest_pairs = calculate_kd_tree_nearest_specific_carbon_pairs_unique(molecules, specific_carbon_labels)

        # Step 3: Save the unique nearest pairs to a file
        output_file_unique = 'nearest_neighbor_unique.txt'
        save_unique_nearest_neighbor_pairs(output_file_unique, unique_nearest_pairs)
        print(f"Nearest neighbor pairs saved to {output_file_unique}")
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found. Please ensure the file exists.")

import matplotlib.pyplot as plt
import seaborn as sns

# Function to plot the probability density distribution of distances
def plot_probability_density(distances, title):
    """
    Plots the probability density distribution of distances using a KDE (Kernel Density Estimate).
    
    Args:
        distances (list): List of distances to visualize.
        title (str): Title of the plot.
    """
    plt.figure(figsize=(10, 6))  # Set figure size
    sns.kdeplot(distances, fill=True)  # Plot KDE with filled area for better visualization
    plt.title(title)  # Add title
    plt.xlabel('Distance (nm)')  # Label x-axis
    plt.ylabel('Probability Density')  # Label y-axis
    plt.show()  # Display the plot

# Extract distances from the unique pairs
distances = [pair[2] for pair in unique_nearest_pairs]

# Plot the probability density distribution
plot_probability_density(distances, 'Probability Density of Nearest Neighbor Distances')
```

## Saving Nearest-Neighbor Pairs
The save_unique_nearest_neighbor_pairs function saves the results to a file:

Input: A list of pairs and a file path.
Output: A .txt file where each line contains:
Two molecule IDs.
The distance between their nearest carbon atoms.
How It Works:

Iterates through the list of pairs.
Writes each pair and their distance to a file with a precision of six decimal places.

Often you want to extract the molecular pair from the gromac file.
In order to do that you need to create an index file using and select the desired index of the residue following the command
where 2 is the residue groupn ri is the index of the residue
```plaintext
echo -e "2\nri 94 |ri 949\nq" |gmx_mpi make_ndx -f npt3.gro -o index.ndx

```

Before able to extract out the pair of molecules it is important to remove boundary condition and rematch the splitted molecule. To run this

```plaintext
echo -e "0\nq" | gmx_mpi trjconv -f npt3.gro -s npt3.tpr -pbc mol -o output_whole.gro
```

the you can extract out the molecular pair by using teh command

```plaintext
echo -e " r_94_r_949\nq" | gmx_mpi trjconv -f output_whole.gro -s output_whole.gro -n index.ndx -o r_94_r_949.pdb
```
However sometimes you want to create all index for each molecules for ease of extracting molecules later

```plaintext
# File to store input for make_ndx
echo "2" > commands.txt
for i in $(seq 1 1001); do  # Replace N with the total number of residues
    echo "ri $i" >> commands.txt
done
echo "q" >> commands.txt
```

It is not advisable to extract all the molecules, so the first 120 lines are extracted.

```plaintext
# Input file and necessary GROMACS files
input_file="nearest_neighbor_centroids.txt"
gro_file="output_whole.gro"
index_file="index.ndx"

# Extract the first 120 lines from the input file
temp_file="first_120_lines.txt"
head -n 120 "$input_file" > "$temp_file"

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
```



However sometimes we are interested in centroid of the molecules, then the program is modified as below

```python
import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import seaborn as sns

# Function to read and parse a .gro file containing molecular trajectory data
def read_gro_file_with_masses(file_path):
    """
    Reads a .gro file and extracts molecular and atomic information.
    
    Args:
        file_path (str): Path to the .gro file.

    Returns:
        dict: A dictionary where keys are molecule numbers, and values are lists of tuples
              containing atomic coordinates (x, y, z) and atom type.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()[2:]  # Skip the first two header lines

    molecules = {}
    for line in lines:
        if not line.strip():
            continue  # Skip empty lines
        try:
            molecule_number = int(line[0:5])  # Extract molecule number from the first 5 characters
            atom_type = line[10:15].strip()  # Extract atom type
            x, y, z = map(float, [line[20:28].strip(), line[28:36].strip(), line[36:44].strip()])  # Extract x, y, z
            if molecule_number in molecules:
                molecules[molecule_number].append((x, y, z, atom_type))
            else:
                molecules[molecule_number] = [(x, y, z, atom_type)]
        except ValueError:
            continue  # Skip lines that cannot be parsed
    return molecules

# Function to filter specific carbon atoms based on labels
def extract_specific_carbon_atoms(atoms, carbon_labels):
    """
    Filters atoms based on a list of specific carbon labels.
    
    Args:
        atoms (list): List of atom tuples (x, y, z, atom_type).
        carbon_labels (list): List of carbon labels to filter (e.g., ['C27', 'C14']).

    Returns:
        list: Filtered list of atom tuples matching the carbon labels.
    """
    return [atom for atom in atoms if any(c_label in atom[3] for c_label in carbon_labels)]

# Function to compute centroids of specific carbon atoms for each molecule
def compute_molecule_centroids(molecules, carbon_labels):
    """
    Compute the centroid (average position) of specified carbon atoms for each molecule.
    
    Args:
        molecules (dict): Dictionary of molecules with atomic data.
        carbon_labels (list): List of specific carbon labels to include in the centroid calculation.

    Returns:
        dict: A dictionary where keys are molecule IDs, and values are centroid positions (x, y, z).
    """
    centroids = {}
    for mol_id, atoms in molecules.items():
        # Filter atoms to include only the specific carbon atoms
        specific_atoms = extract_specific_carbon_atoms(atoms, carbon_labels)
        if specific_atoms:
            # Compute the centroid as the average of all positions
            centroid = np.mean([atom[:3] for atom in specific_atoms], axis=0)
            centroids[mol_id] = centroid
    return centroids

# Function to compute nearest-neighbor distances based on centroids
def calculate_kd_tree_nearest_molecule_distances(centroids):
    """
    Computes unique nearest-neighbor distances between molecular centroids.
    
    Args:
        centroids (dict): Dictionary of molecule IDs and their centroid positions.

    Returns:
        list: List of unique pairs (mol_id1, mol_id2, distance).
    """
    # Extract centroids into an array and track molecule IDs
    mol_ids = np.array(list(centroids.keys()))
    points = np.array(list(centroids.values()))

    # Build KD-tree for centroid positions
    tree = cKDTree(points)
    distances, indexes = tree.query(points, k=2)  # Nearest and second nearest neighbors

    nearest_pairs = set()  # Use a set to avoid duplicate pairs
    for dist, idx, mol_id in zip(distances[:, 1], indexes[:, 1], mol_ids):
        neighbor_mol_id = mol_ids[idx]
        if mol_id != neighbor_mol_id:  # Exclude self-pairs
            # Add pair as a sorted tuple to ensure uniqueness
            pair = tuple(sorted((mol_id, neighbor_mol_id)))
            nearest_pairs.add((pair[0], pair[1], dist))
    return list(nearest_pairs)

# Function to save nearest-neighbor pairs to a file
def save_unique_nearest_neighbor_pairs(file_name, pairs):
    """
    Save unique molecular pairs with their distances to a file.
    
    Args:
        file_name (str): Path to save the file.
        pairs (list): List of unique pairs with distances.
    """
    with open(file_name, 'w') as file:
        for pair in pairs:
            file.write(f"{pair[0]} {pair[1]} {pair[2]:.6f}\n")

# Function to plot the probability density distribution of distances
def plot_probability_density(distances, title):
    """
    Plots the probability density distribution of distances using a KDE (Kernel Density Estimate).
    
    Args:
        distances (list): List of distances to visualize.
        title (str): Title of the plot.
    """
    plt.figure(figsize=(10, 6))  # Set figure size
    sns.kdeplot(distances, fill=True)  # Plot KDE with filled area for better visualization
    plt.title(title)  # Add title
    plt.xlabel('Distance (nm)')  # Label x-axis
    plt.ylabel('Probability Density')  # Label y-axis
    plt.show()  # Display the plot

# Main program
if __name__ == "__main__":
    file_path = 'npt3.gro'  # Path to the .gro file
    specific_carbon_labels = ['C27', 'C14', 'C12', 'C2', 'C8', 'C5', 'C15', 'C18',
                              'C22', 'C29', 'C32', 'C35', 'C40', 'C56', 'C44', 'C53', 'C58', 'C62'] # BEWARE the GauView Number X to PDB numnber

    try:
        # Step 1: Read and parse the molecular data
        molecules = read_gro_file_with_masses(file_path)

        # Step 2: Compute centroids for each molecule
        centroids = compute_molecule_centroids(molecules, specific_carbon_labels)

        # Step 3: Calculate unique nearest-neighbor molecular distances
        unique_nearest_pairs = calculate_kd_tree_nearest_molecule_distances(centroids)

        # Step 4: Save unique nearest-neighbor pairs to a file
        output_file_unique = 'nearest_neighbor_centroids.txt'
        save_unique_nearest_neighbor_pairs(output_file_unique, unique_nearest_pairs)
        print(f"Nearest neighbor pairs saved to {output_file_unique}")

        # Step 5: Extract distances and plot probability density
        distances = [pair[2] for pair in unique_nearest_pairs]
        plot_probability_density(distances, 'Probability Density of Nearest Neighbor Distances (Centroids)')
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found. Please ensure the file exists.")
```


To create .xyz file for viewing

```plaintext
# Input text file containing the file pairs
input_list="nearest_neighbor_centroids.txt"

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
head -n 120 "$input_list" | while read -r line; do
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
```

# You need to filter out those which has the good overall contact using the code below before submitting to Gaussian

```phython
import os
import numpy as np

def extract_coordinates(file_path, indices):
    """
    Extracts coordinates for specified atom indices from an XYZ file.

    Args:
        file_path (str): Path to the XYZ file.
        indices (list): List of atom indices to extract (1-based indexing).

    Returns:
        list: A list of (x, y, z) coordinates as tuples of floats.
    """
    coordinates = []
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            for index in indices:
                # Adjust for 0-based indexing in Python
                line = lines[index + 1].strip().split()
                x, y, z = map(float, line[1:4])
                coordinates.append((x, y, z))
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
    return coordinates

def calculate_average_distance(coords1, coords2):
    """
    Calculates the average distance between two lists of coordinates.

    Args:
        coords1 (list): List of (x, y, z) tuples for the first set of points.
        coords2 (list): List of (x, y, z) tuples for the second set of points.

    Returns:
        float: The average distance between the points in the two lists.
    """
    distances = []
    for coord1 in coords1:
        for coord2 in coords2:
            distance = np.linalg.norm(np.array(coord1) - np.array(coord2))
            distances.append(distance)
    return np.mean(distances)

if __name__ == "__main__":
    # Indices for the first and second molecules
    molecule1_indices = [2, 7, 12, 14, 16, 18, 23, 29, 34, 39, 44, 53, 58, 65, 70, 75, 77, 78, 83, 90, 92, 97]
    molecule2_indices = [104, 109, 114, 116, 118, 120, 125, 131, 136, 141, 146, 155, 160, 167, 172, 177, 179, 180, 185, 192, 194, 199]

    # Output file to store results
    output_file = "average_distances.txt"

    # Open the output file
    with open(output_file, "w") as out_file:
        out_file.write("Filename\tAverage Distance\n")

        # Loop through all .xyz files in the current directory
        for filename in os.listdir("."):
            if filename.endswith(".xyz"):
                try:
                    # Extract coordinates for the specified indices
                    molecule1_coords = extract_coordinates(filename, molecule1_indices)
                    molecule2_coords = extract_coordinates(filename, molecule2_indices)

                    # Calculate the average distance
                    average_distance = calculate_average_distance(molecule1_coords, molecule2_coords)

                    # Write the results to the output file
                    out_file.write(f"{filename}\t{average_distance:.4f}\n")
                    print(f"Processed {filename}: Average Distance = {average_distance:.4f}")
                except Exception as e:
                    print(f"Error processing file {filename}: {e}")

    print(f"Results saved to {output_file}.")
```
# you need to match the filtered molecular pair with the nearest_neighbor_distances.txt and updated that distance
or run

```python
import pandas as pd
# Loading the contents of both files into DataFrames
file_a_path = 'a.txt'
file_b_path = 'b.txt'

# Reading files into DataFrames
df_a = pd.read_csv(file_a_path, sep="\t", header=None, names=["Col1", "Col2", "Col3"])
df_b = pd.read_csv(file_b_path, sep="\t", header=None, names=["Col1", "Col2", "Col3"])

# Find matching rows based on Col1 and Col2
matching_rows = pd.merge(df_a, df_b, on=["Col1", "Col2"])

# Save the matching rows to a new file
matching_rows.to_csv('matching_rows.txt', sep="\t", index=False, header=False)
```



run https://github.com/ph7klw76/gromac_note/blob/main/run11.sh
to create all above 

to create gaussian file to coupling the electronic coupling run the script code below

```plaintext
#!/bin/bash -l
  
#SBATCH --partition=cpu-epyc-genoa
#SBATCH --job-name=extract
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:59:59

 
# Load dependencies and libraries
source /app/hpcx/2.17.1/hpcx-init.sh
hpcx_load
 
module load gromacs/2023.2
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}


#!/bin/bash

# Input text file containing the file pairs 


input_list="nearest_neighbor_centroids.txt"

# Output directory for the generated .gjf files
output_dir="output_gjf_files"
mkdir -p "$output_dir"  # Create the directory if it doesn't exist

# Summary file for processed pairs
summary_file="summary.txt"
> "$summary_file"  # Initialize the summary file

# Function to extract unique first letters from a PDB file
extract_unique_letters() {
    local input_file="$1"
    awk '
    $1 == "ATOM" {
        a = substr($3, 1, 1); # Extract first letter of the 3rd column
        print a
    }' "$input_file" | sort -u
}

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
head -n 120 "$input_list" | while read -r line; do
    # Parse the file names from each line
    xx=$(echo "$line" | awk '{print $1}')
    xxx=$(echo "$line" | awk '{print $2}')
    
    file1="${xx}.pdb"
    file2="${xxx}.pdb"
    
    # Define the output .gjf files
    mkdir -p "$output_dir/${xx}_${xxx}"  # Create the directory if it doesn't exist

    output_pair_file="${output_dir}/${xx}_${xxx}/${xx}_${xxx}.gjf"
    output_file1="${output_dir}/${xx}_${xxx}/${xx}.gjf"
    output_file2="${output_dir}/${xx}_${xxx}/${xxx}.gjf"

    # Check if input files exist
    if [[ ! -f "$file1" || ! -f "$file2" ]]; then
        echo "WARNING: One or both files missing: $file1, $file2" >> "$summary_file"
        continue
    fi

    # Extract data and unique letters
    data1=$(process_pdb "$file1")
    data2=$(process_pdb "$file2")
    unique_letters1=$(extract_unique_letters "$file1" | tr '\n' ' ')
    unique_letters2=$(extract_unique_letters "$file2" | tr '\n' ' ')
    # Combine unique letters from both files
    B=$(echo -e "$unique_letters1\n$unique_letters2" | sort -u | tr '\n' ' ')

    # Write the pair .gjf file
    {
        echo "%nprocshared=16"
        echo "%mem=32GB"
        echo "# gen guess=huckel nosymm pop=nboread"
        echo "# scf=(direct,nosymm)"
        echo ""
        echo "Cb 1"
        echo ""
        echo "0 1"
        echo "$data1"
        echo "$data2"
        echo ""
        echo "$B 0"
        echo "Def2SVP"
        echo "****"
        echo ""
        echo "\$NBO SAO=w53 FAO=W54 \$END"
    } > "$output_pair_file"

    # Write the .gjf file for file1
    {
        echo "%nprocshared=16"
        echo "%mem=32GB"
        echo "# gen nosymm punch(MO)"
        echo "# scf=(direct,nosymm)"
        echo ""
        echo "carbazol"
        echo ""
        echo "0 1"
        echo "$data1"
        echo ""
        echo "$unique_letters1 0"
        echo "Def2SVP"
        echo "****"
        echo ""
    } > "$output_file1"

    # Write the .gjf file for file2
    {
        echo "%nprocshared=16"
        echo "%mem=32GB"
        echo "# gen nosymm punch(MO)"
        echo "# scf=(direct,nosymm)"
        echo ""
        echo "carbazol"
        echo ""
        echo "0 1"
        echo "$data2"
        echo ""
        echo "$unique_letters2 0"
        echo "Def2SVP"
        echo "****"
        echo ""
    } > "$output_file2"

    # Log to the summary file
    echo "Processed: $file1 and $file2" >> "$summary_file"
    echo " - Pair .gjf file: $output_pair_file" >> "$summary_file"
    echo " - Single .gjf file for $file1: $output_file1" >> "$summary_file"
    echo " - Single .gjf file for $file2: $output_file2" >> "$summary_file"

    # Print confirmation
    echo "Processed: $file1 and $file2"
done

# Print completion message
echo "Processing complete. Summary saved in $summary_file."

```

After the file are created you can run the script below:
The will auto submit the jobs from the folder accrodingly

```plaintext
#!/bin/bash -l

#SBATCH --partition=cpu-epyc-genoa
#SBATCH --job-name=submit
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --qos=long
#SBATCH --time=2-23:59:59
#SBATCH --hint=nomultithread

# Input file
INPUT_FILE="nearest_neighbor_centroids.txt"

# Starting line (n) and number of lines to process (m)
START_LINE=1  # Replace with the desired starting line
NUM_LINES=120  # Replace with the number of lines to process

# Read each line in the file
LINE_COUNT=0
while IFS=' ' read -r FIRST SECOND _; do
    # Increment the line counter
    LINE_COUNT=$((LINE_COUNT + 1))
    
    # Process only lines from START_LINE to START_LINE + NUM_LINES - 1
    if [ "$LINE_COUNT" -ge "$START_LINE" ] && [ "$LINE_COUNT" -lt $((START_LINE + NUM_LINES)) ]; then
        # Construct the Gaussian commands
        DIR="${FIRST}_${SECOND}"
        mkdir -p "$DIR"
        sleep 1

        # Create SLURM submission script inside the directory
        SUBMISSION_FILE="${DIR}/submit_gaussian.sh"
        cat > "$SUBMISSION_FILE" << EOF
#!/bin/bash
#SBATCH --partition=cpu-epyc-genoa
#SBATCH --job-name=${FIRST}_${SECOND}
#SBATCH --output=job.out
#SBATCH --error=job.err
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH --ntasks=16
#SBATCH --qos=long
#SBATCH --time=2-23:59:59
#SBATCH --hint=nomultithread

module load gaussian/g09
source \$g09profile

g09 <${FIRST}.gjf> ${FIRST}.log
sleep 3
mv fort.7 fort.7-${FIRST}
g09 <${FIRST}_${SECOND}.gjf> ${FIRST}_${SECOND}.log
g09 <${SECOND}.gjf> ${SECOND}.log
sleep 3
mv fort.7 fort.7-${SECOND}
EOF

        # Make the submission script executable
        chmod +x "$SUBMISSION_FILE"

        # Submit the job from within the directory
        cd "$DIR"
        JOB_ID=$(sbatch "submit_gaussian.sh" | awk '{print $NF}')
        echo "Submitted job with ID: $JOB_ID from directory: $DIR"
        cd ..
    fi

    # Stop processing if we exceed the range
    if [ "$LINE_COUNT" -ge $((START_LINE + NUM_LINES)) ]; then
        break
    fi
done < "$INPUT_FILE"

```
in order to apply the python code at https://github.com/ph7klw76/gaussian_note/blob/main/electronic_coupling.md
we need to create an paramFile.txt

```plaintext
orbitalsA=131 132
orbitalsB=131 132
nBasisFunctsA=676   
nBasisFunctsB=676  
```
BasisFunctions1 represents the total number of basis functions used to calculate Molecule1, which can be obtained from fort.7 by identifying the number of coefficients in the orbital at 1 Alpha.
For orbitalsA, just write down the HOMO index . HOMO-1 will take into play when the energu between HOMO and HOMO-1 is less than 0.3eV
Once the calculation completed you have a given pair will contain data below

![image](https://github.com/user-attachments/assets/3e30bc38-0bbf-4225-84f1-247777017657)

based on the code for single pair at https://github.com/ph7klw76/gaussian_note/blob/main/electronic_coupling.md 
It can be extended to become below for automated calculation as below

```python
import shutil
import os
from pathlib import Path
import ECOUPLING as EC
import glob

# Define the current directory
current_directory = os.getcwd()

# Read the 'nearest_neighbor_centroids.txt' file and process line by line
with open('nearest_neighbor_centroids.txt', 'r') as f1:
    for line in f1:
        # Split line into components and retrieve necessary values
        line_parts = line.split()
        if len(line_parts) < 2:
            continue  # Skip malformed lines (if any)
        
        line0, line1 = line_parts[0], line_parts[1]
        
        # Set up source folder and files to copy
        source_folder = f'{line0}_{line1}'
        file_to_copy = ['FILE.53', 'FILE.54', f'fort.7-{line0}', f'fort.7-{line1}']
        
        # Copy files from the source folder to the current directory
        for file_name in file_to_copy:
            source_file_path = Path(source_folder) / file_name
            try:
                shutil.copy(source_file_path, current_directory)
                print(f"Copied {source_file_path} to {current_directory}")
            except Exception as e:
                print(f"Error copying {source_file_path}: {e}")

        # Create or overwrite 'inFile.in'
        with open('inFile.in', 'w') as f2:
            f2.write('FILE.53\n')
            f2.write('FILE.54\n')
            f2.write(f'fort.7-{line0}\n')
            f2.write(f'fort.7-{line1}\n')

        # Call ECOUPLING method
        try:
            EC.coupling(['a', 'inFile.in', 'paramFile.txt', f'output{line0}_{line1}'])
            print(f"Coupling executed for {line0}_{line1}")
        except Exception as e:
            print(f"Error executing ECOUPLING for {line0}_{line1}: {e}")

# Remove all files starting with 'fort'
files_to_remove = glob.glob('fort*')
for file in files_to_remove:
    try:
        os.remove(file)
        print(f"Removed: {file}")
    except Exception as e:
        print(f"Error removing {file}: {e}")
```



## Conclusion
This Python-based framework showcases the synergy of computational geometry and data visualization in molecular simulations.
By leveraging KD-tree for efficient spatial queries and KDE for insightful visualizations, the script provides a robust tool for analyzing atomic-scale spatial relationships in molecular systems.
