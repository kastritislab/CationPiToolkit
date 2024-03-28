# CationPiToolkit

This toolkit analyses a protein structure, given as a PDB file, for potential cation pi interactions. It is recommended to visually inspect the identified sites of interaction.

## Installation

Clone the repository using git:

```
git clone https://github.com/kastritislab/CationPiToolkit
```

Navigate to the downloaded folder:

```
cd CationPiToolkit
```

and install the toolkit:

```
pip install .
```

This will install the `CationPiToolkit` command to your path.

## Usage

To use the toolkit, run the `CationPiToolkit` command followed by the path to your pdb file, and any other optional arguments:

```
CationPiToolkit "/path/to/file.pdb" 
```

This will analyze the given pdb file for potential cation pi interactions and save the result as a csv file.

Optional arguments:

* `--residues`: A list of residues to consider for interactions.
  * Default: `['LYS', 'ARG', 'PHE', 'TRP', 'TYR']`
  * Usage: `--residues LYS CYS PHE`

* `--exclude_backbone`: If set, exclude backbone atoms.
  * Default: `False`
  * Usage: `--exclude_backbone` (without value sets this to True)

* `--exclude_atoms`: A list of atoms to exclude from the results.
  * Default: `["CB", "NH1", "NH2", "NE1", "NE2", "OH"]`
  * Usage: `--exclude_atoms CA CB NH1 NH2`

* `--bait_atoms`: A list of bait atoms for the interactions. Each bait atom is represented as a tuple of residue name and atom name.
  * Default: `[('ARG', 'CZ'), ('LYS', 'NZ')]`
  * Usage: `--bait_atoms CYS CA --bait_atoms SER OH`

* `--prey_atoms`: A list of prey atoms for the interactions. Each prey atom is represented as a tuple of residue name and a list of atom names.
  * Default: `[('TRP', ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']), ('PHE', ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']), ('TYR', ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'])]`
  * Usage: `--prey_atoms GLY CA N O C --prey_atoms CYS CB SD`

* `--chains`: A list of chains to consider for interactions.
  * Default: `None`
  * Usage: `--chains A a B b`

* `--min_interactions`: Minimum number of interactions for consideration.
  * Default: `6`
  * Usage: `--min_interactions 10`

* `--mean_threshold`: Mean distance threshold.
  * Default: `5`
  * Usage: `--mean_threshold 4.5`

* `--std_threshold`: Standard deviation threshold.
  * Default: `0.75`
  * Usage: `--std_threshold 1.0`

These are intended as reasonable default values, but you may need to adjust them for your specific use case.


## Acknoledgment

The writing of this software was assited by ChatGPT.

## Citation

If this software is used in your research, please cite the following paper: 
`to be published`

