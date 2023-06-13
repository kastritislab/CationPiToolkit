import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm

def parse_residues(path, atoms = None, residues = None, exclude_backbone = False, exclude_atoms = None):
    """
    Parses PDB file and filters atoms based on the provided parameters.
    
    Parameters:
        path (str): The path to the PDB file.
        atoms (list, optional): A list of atoms to include in the results. If None, includes all atoms.
        residues (list, optional): A list of residues to include in the results. If None, includes all residues.
        exclude_backbone (bool, optional): If True, exclude backbone atoms.
        exclude_atoms (list, optional): A list of atoms to exclude from the results.
    
    Returns:
        pd.DataFrame: A DataFrame of parsed residues and atoms.
    """
    
    with open(path) as pdb:
        lines = pdb.readlines()
    
    # filter them
    parsed = []
    for line in lines:
        
        if not line.startswith("ATOM"):
            continue

        atom = line[12:16].strip()
        resn = line[17:20].strip()
        resi = int(line[22:26])
        chain = line[20:22].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        
        if exclude_backbone:
            if atom in ['N', 'O', 'C', 'CA']:
                continue
        
        if atoms is not None:
            if atom not in atoms:
                continue
                
        if residues is not None:
            if resn not in residues:
                continue
        
        if exclude_atoms:
            if atom in exclude_atoms:
                continue
            
        
        parsed.append((atom, resn, resi, chain, (x, y, z)))
    
    return pd.DataFrame(parsed, columns = ['atom', 'resn', 'resi', 'chain', 'coords'])

def get_distances(parsed_structure, bait_atoms, prey_atoms, distance_cutoff):
    parsed_structure["coordsx"] = list(map(np.array, parsed_structure.coords))

    # filter bait atoms
    bait_atoms_set = set(bait_atoms)
    bait = parsed_structure[parsed_structure[['resn', 'atom']].apply(tuple, axis=1).isin(bait_atoms_set)]

    # filter the prey atoms
    prey_atoms_dict = {resn: set(atoms) for resn, atoms in prey_atoms}
    mask = parsed_structure.apply(lambda row: row['atom'] in prey_atoms_dict.get(row['resn'], set()), axis=1)
    prey = parsed_structure[mask].copy()

    measured_distances = []

    for row in tqdm(bait.itertuples(), total=len(bait)):
        _, atom1, resn1, resi1, chain1, _, c1 = row

        c2 = np.vstack(prey.coordsx.to_numpy())

        # calculate the distances
        distances = np.sqrt(np.sum((c2 - c1) ** 2, axis=1))
        prey["dists"] = distances

        # make a mask
        mask = distances < distance_cutoff

        # get the residues in distance
        hits = prey[mask].iloc[:, [0, 1, 2, 3, 6]].copy()

        hits["cation_atom"] = atom1
        hits["cation_resn"] = resn1
        hits["cation_resi"] = resi1
        hits["cation_chain"] = chain1

        # rename the result
        rename = {'atom': 'pi_atom', 'resn': 'pi_resn', 'resi': 'pi_resi', 'chain': 'pi_chain'}
        hits.rename(columns=rename, inplace=True)

        # order them
        order = ['cation_atom', 'cation_resn', 'cation_resi', 'cation_chain', 'pi_atom', 'pi_resn', 'pi_resi',
                 'pi_chain', 'dists', ]

        measured_distances.append(hits[order])

    return pd.concat(measured_distances)

def catpi_finder(path, residues=['LYS', 'ARG', 'PHE', 'TRP', 'TYR'], 
                 exclude_backbone=True, exclude_atoms=["CB", "NH1", "NH2", "NE1", "NE2", "OH"], 
                 bait_atoms=[('ARG', 'CZ'), ('LYS', 'NZ')], 
                 prey_atoms=[('TRP', ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']),
                             ('PHE', ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']),
                             ('TYR', ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'])], 
                 chains=None, min_interactions=5, mean_threshold=5, std_threshold=0.75):
    """
    Executes the default workflow for identifying cation-pi interactions in a protein structure.

    Parameters:
        path (str): Path to the PDB file.
        residues (list of str, optional): Residues to include in the search for interactions. Defaults to ['LYS', 'ARG', 'PHE', 'TRP', 'TYR'].
        exclude_backbone (bool, optional): Whether to exclude backbone atoms. Defaults to True.
        exclude_atoms (list of str, optional): Specific atoms to exclude from the search. Defaults to ["CB", "NH1", "NH2", "NE1", "NE2", "OH"].
        bait_atoms (list of tuples, optional): Atoms to consider as 'bait' in the interaction. Defaults to [('ARG', 'CZ'), ('LYS', 'NZ')].
        prey_atoms (list of tuples, optional): Atoms to consider as 'prey' in the interaction. Each tuple contains a residue name and a list of atom names. Defaults to [('TRP', ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']), ('PHE', ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']), ('TYR', ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'])].
        chains (list of str, optional): Specific chains to consider for interactions. If None, all chains are considered. Defaults to None.
        min_interactions (int, optional): Minimum number of interactions required for consideration. Defaults to 5.
        mean_threshold (float, optional): Mean distance threshold for considering an interaction. Defaults to 5.
        std_threshold (float, optional): Standard deviation threshold for considering an interaction. Defaults to 0.75.

    Returns:
        pandas.DataFrame: A DataFrame containing the potential interactions identified.
    """

    
    # parse the pdb structure
    parsed_structure = parse_residues(path, residues = residues, exclude_backbone = exclude_backbone, exclude_atoms = exclude_atoms)
    
    # get the distance
    distances = get_distances(parsed_structure, bait_atoms, prey_atoms, distance_cutoff = mean_threshold * 1.4)
    
    # filter the distances according the input variables
    if chains is not None:
        df = distances.loc[(distances.cation_chain.isin(chains)) | (distances.pi_chain.isin(chains)) ]
    else:
        df = distances
        
    # aggregated_distances
    aggregated_distances = df.groupby(['cation_resn', 'cation_resi', 'cation_chain', 'pi_resn', 'pi_resi', 'pi_chain']).agg({'dists': ['mean', 'std'], 'cation_resn': 'size'})
    aggregated_distances.columns = ['mean_dist', 'std_dist', 'n']
    aggregated_distances = aggregated_distances.reset_index()


    # filter accoring the input variable
    selected_distances = aggregated_distances.loc[(aggregated_distances.mean_dist < mean_threshold) & (aggregated_distances.std_dist < std_threshold) & (aggregated_distances.n >= min_interactions)]
    
    # and sort the dataframe
    selected_distances = selected_distances.sort_values(by=['cation_chain', 'cation_resi']).reset_index(drop = True)

    
    # save the interactions
    ofile = path[:-4] + "_catpi.csv"
    selected_distances.to_csv(ofile, index = False)
    
    return selected_distances

def main():
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description='Process PDB files to find cation-pi interactions.')

    # Add arguments to the parser
    parser.add_argument('path', type=str, help='Path to the PDB file.')
    parser.add_argument('--residues', nargs='+', default=['LYS', 'ARG', 'PHE', 'TRP', 'TYR'], help='List of residues to consider for interactions.')
    parser.add_argument('--exclude_backbone', action='store_true', help='If set, exclude backbone atoms.')
    parser.add_argument('--exclude_atoms', nargs='+', default=["CB", "NH1", "NH2", "NE1", "NE2", "OH"], help='List of atoms to exclude from the results.')
    parser.add_argument('--chains', nargs='*', default=None, help='List of chains to consider for interactions.')
    parser.add_argument('--min_interactions', type=int, default=6, help='Minimum number of interactions for consideration.')
    parser.add_argument('--mean_threshold', type=float, default=5, help='Mean distance threshold.')
    parser.add_argument('--std_threshold', type=float, default=0.75, help='Standard deviation threshold.')
    parser.add_argument('--bait_atoms', nargs='+', default=[('ARG', 'CZ'), ('LYS', 'NZ')], help='List of bait atoms for the interactions. Each bait atom is represented as a tuple of residue name and atom name.')
    parser.add_argument('--prey_atoms', nargs='+', default=[('TRP', ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']),
                                                            ('PHE', ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']),
                                                            ('TYR', ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'])], 
                       help='List of prey atoms for the interactions. Each prey atom is represented as a tuple of residue name and a list of atom names.')

    # Parse the arguments
    args = parser.parse_args()

    # Use the arguments in the default_workflow function
    interactions = catpi_finder(
        args.path, 
        residues=args.residues,
        exclude_backbone=args.exclude_backbone,
        exclude_atoms=args.exclude_atoms,
        bait_atoms=args.bait_atoms,
        prey_atoms=args.prey_atoms,
        chains=args.chains,
        min_interactions=args.min_interactions,
        mean_threshold=args.mean_threshold,
        std_threshold=args.std_threshold
    )

    

# If the script is run directly (python CationPiToolkit.py), call main()
if __name__ == "__main__":
    main()  