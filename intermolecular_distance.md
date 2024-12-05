# Analyzing Carbon Atom Spatial Relationships in Molecular Simulations Using KD-Tree and Python

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
                          'C22', 'C29', 'C32', 'C35', 'C40', 'C56', 'C44', 'C53', 'C58', 'C62']
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
                              'C22', 'C29', 'C32', 'C35', 'C40', 'C56', 'C44', 'C53', 'C58', 'C62']

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
```plaintext

# Run gmx_mpi make_ndx
gmx_mpi make_ndx -f npt3.gro -o index.ndx < commands.txt


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
                              'C22', 'C29', 'C32', 'C35', 'C40', 'C56', 'C44', 'C53', 'C58', 'C62']

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


## Conclusion
This Python-based framework showcases the synergy of computational geometry and data visualization in molecular simulations.
By leveraging KD-tree for efficient spatial queries and KDE for insightful visualizations, the script provides a robust tool for analyzing atomic-scale spatial relationships in molecular systems.
