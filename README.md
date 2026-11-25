# Template-based Crystal Structure Substitution Generator

A pymatgen-based tool to generate candidate crystal structures by substituting elements on Wyckoff sites or original element groups from predefined structural templates.

## Overview

This code performs template-matching and site substitution for crystal structure prediction:

- Load template structures stored in a `template.csv` file (POSCAR string + metadata)
- Two substitution modes:
  - `wyckoff`: group atoms by symmetry-equivalent Wyckoff sites
  - `element`: group atoms by the original element symbols in the template
- Match total atom count, number of unique sites / element types with the target composition
- Enumerate valid assignments of site groups to match the target stoichiometry exactly
- Rescale the cell volume using `DLSVolumePredictor`
- Export all generated candidate structures as VASP POSCAR (`*.vasp`) files

## Dependencies

```
pip install pymatgen numpy pandas
```

Python ≥ 3.8 recommended.

## File Layout

```
.
├── template
│   └── template.csv       # Template library: atom_num, unique_site, element_num, structure_str (poscar)
├── results                 # Default output folder for generated structures
└── substitute.py          # Main substitution class and helper functions
```

## Template CSV Format

| Column | Description |
| --- | --- |
| `atom_num` | Total number of atoms in the template structure |
| `unique_site` | Number of symmetry-unique Wyckoff sites |
| `element_num` | Number of distinct elements in the template |
| `structure_str` | Full POSCAR string of the template |

## Quick Usage

```
from substitute import Substitute

sub = Substitute(max_atom=30)

# Generate candidates by element grouping
sub.get_templates_from_composition(
    target_formula="UB2",
    mode="element",
    atom_num_match=False,
    set_num_match=True,
    save_path="./results"
)

# Generate candidates by Wyckoff site grouping
sub.get_templates_from_composition(
    target_formula="UB2",
    mode="wyckoff",
    atom_num_match=False,
    set_num_match=True,
    save_path="./results"
)
```

### Parameters

- `target_formula`: Chemical formula string, e.g. `"UB2"`, `"MgO"`
- `mode`: `"wyckoff"` or `"element"`
- `atom_num_match`: `True` requires exact total atom match; `False` allows integer multiples
- `set_num_match`: `True` requires exact match of unique-site / element count; `False` allows template count ≥ target
- `save_path`: Root directory to write output VASP files


## Citation

If you use this code in your research, please cite pymatgen and the corresponding DLS volume predictor work.

## License

MIT