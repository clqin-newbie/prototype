from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.structure_prediction.volume_predictor import DLSVolumePredictor
from pathlib import Path
import pandas as pd 
import numpy as np 
class Substitute:
    """
    A class to perform structure substitution based on template structures.
    It matches target compositions with pre-stored templates and generates substituted crystal structures.
    """
    def __init__(self, max_atom=10):
        """
        Initialize the Substitute class.

        Args:
            max_atom (int): Maximum total number of atoms allowed for loaded template entries.
        """
        self.template = pd.read_csv("./template/template.csv")
        self.template = self.template[self.template['atom_num'] <= max_atom]

    def _template_match(self, target_composition, mode, set_num_match, atom_num_match):
        """
        Match suitable template rows according to the target composition and matching rules.

        Args:
            target_composition (Composition): Target chemical composition.
            mode (str): Matching mode, either 'wyckoff' or 'element'.
            set_num_match (bool): If True, require exact match of unique-site / element count;
                                  if False, require template count >= target count.
            atom_num_match (bool): If True, require exact total atom count match;
                                   if False, allow template atom count to be an integer multiple.

        Returns:
            pd.DataFrame: Filtered template dataframe.
        """
        # Exact match of total atom count
        if atom_num_match:
            template = self.template[self.template['atom_num'] == target_composition.num_atoms]
        else:
            # Template atom count is an integer multiple of target atom count
            template = self.template[self.template['atom_num'] % target_composition.num_atoms == 0]

        if mode == 'wyckoff':
            if set_num_match:
                template = template[template['unique_site'] == len(target_composition.as_dict())]
            else:
                template = template[template['unique_site'] >= len(target_composition.as_dict())]
        else:
            if set_num_match:
                template = template[template['element_num'] == len(target_composition.as_dict())]
            else:
                template = template[template['element_num'] >= len(target_composition.as_dict())]

        template = template.reset_index(drop=True)
        return template

    def get_templates_from_composition(self, target_formula, mode='element', atom_num_match=False, set_num_match=True, save_path='./results'):
        """
        Generate substituted structures from a given chemical formula and save them as VASP POSCAR files.

        Args:
            target_formula (str): Chemical formula of the target phase, e.g. "UB2".
            mode (str): Matching strategy: 'element' (match by original elements in template)
                        or 'wyckoff' (match by symmetry-equivalent Wyckoff sites).
            atom_num_match (bool): Whether to enforce exact total atom number match.
            set_num_match (bool): Whether to enforce exact match of unique-site / element count.
            save_path (str): Root directory for output VASP files.
        """
        target_composition = Composition(target_formula)

        save_path = Path(save_path).joinpath(target_formula)
        save_path.mkdir(exist_ok=True, parents=True)
        template = self._template_match(target_composition, mode=mode, set_num_match=set_num_match, atom_num_match=atom_num_match)

        num = 0
        for idx in range(template.shape[0]):
            template_sturcture = Structure.from_str(template['structure_str'][idx], fmt='poscar')
            if mode == 'wyckoff':
                s_list = self._get_substitue_struc_wyckoff(target_composition, template_sturcture)
            else:
                s_list = self._get_substitue_struc_element(target_composition, template_sturcture)
            for struc in s_list:
                struc.to_file(save_path.joinpath(f"{num+1}.vasp"), fmt='poscar')
                num += 1

    def _get_substitue_struc_wyckoff(self, target_composition, template_sturcture):
        """
        Generate substituted structures using symmetry-equivalent Wyckoff site groups.

        Args:
            target_composition (Composition): Target chemical composition.
            template_sturcture (Structure): Template crystal structure.

        Returns:
            list[Structure]: List of substituted candidate structures.
        """
        template_composition = template_sturcture.composition
        template_atom_num = template_composition.num_atoms
        target_atom_num = target_composition.num_atoms
        sga = SpacegroupAnalyzer(template_sturcture)
        sym_data = sga.get_symmetry_dataset()
        eq_atoms = sym_data["equivalent_atoms"]
        unique_value = np.unique(eq_atoms)
        indices_num = []
        indices_list = []
        for i in unique_value:
            indices = tuple([int(j) for j in np.where(eq_atoms == i)[0]])
            indices_list.append(indices)
            indices_num.append(len(indices))
        struc_list = self._sub(target_composition, template_atom_num, target_atom_num, indices_num, indices_list, template_sturcture)
        return struc_list

    def _sub(self, target_composition, template_atom_num, target_atom_num, indices_num, indices_list, template_sturcture):
        """
        Core substitution routine: assign target elements to grouped site indices and rescale volume.

        Args:
            target_composition (Composition): Target composition.
            template_atom_num (int): Total atoms in template.
            target_atom_num (int): Total atoms in target formula.
            indices_num (list[int]): Atom count for each site group.
            indices_list (list[tuple]): Atom index tuples for each site group.
            template_sturcture (Structure): Template structure object.

        Returns:
            list[Structure]: Valid substituted and scaled structures.
        """
        target_composition *= template_atom_num / target_atom_num
        target_num = list(target_composition.as_dict().values())
        target_elem = list(target_composition.as_dict().keys())
        data_list = self._can_form_exact(indices_num, indices_list, target_num=target_num, ele_num=len(target_num))

        struc_list = []
        for indices in data_list:
            target_structure = template_sturcture.copy()
            for index, elem in zip(indices, target_elem):
                for i in index:
                    target_structure[i] = elem
            struc_list.append(scale_structure(target_structure.get_sorted_structure()))
        return struc_list

    def _get_substitue_struc_element(self, target_composition, template_sturcture):
        """
        Generate substituted structures by grouping atoms according to their original element symbols in template.

        Args:
            target_composition (Composition): Target chemical composition.
            template_sturcture (Structure): Template crystal structure.

        Returns:
            list[Structure]: List of substituted candidate structures.
        """
        template_composition = template_sturcture.composition
        template_atom_num = template_composition.num_atoms
        target_atom_num = target_composition.num_atoms
        indices_num = []
        indices_list = []
        for e in template_sturcture.composition.elements:
            indices = template_sturcture.indices_from_symbol(e.symbol)
            indices_list.append(indices)
            indices_num.append(len(indices))
        struc_list = self._sub(target_composition, template_atom_num, target_atom_num, indices_num, indices_list, template_sturcture)
        return struc_list

    def _can_form_exact(self, indices_num, indices_list, target_num, ele_num):
        """
        Find all ordered partitions of site groups such that summed atom counts exactly match target stoichiometry.
        All source groups must be consumed without duplication.

        Args:
            indices_num (list[int]): Atom count of each source site group.
            indices_list (list[tuple]): Atom index tuples of each source site group.
            target_num (list[int]): Required atom counts for each target element.
            ele_num (int): Number of distinct target elements.

        Returns:
            list[list[tuple]]: Valid index-group assignments matching target stoichiometry.
        """
        data_list = []
        idx = list(range(len(indices_num)))
        comb_list = partitions_ordered_nonempty(idx, ele_num)
        for comb in comb_list:
            new_indices_num = []
            new_indices_list = []
            for data in comb:
                num = 0
                indices = ()
                for i in data:
                    num += indices_num[i]
                    indices += indices_list[i]
                new_indices_list.append(indices)
                new_indices_num.append(num)
            if new_indices_num == target_num:
                data_list.append(new_indices_list)
        return data_list


def partitions_ordered_nonempty(lst, m):
    """
    Partition the input list into m ordered, non-empty groups.
    Group order matters; every group must contain at least one element.

    Args:
        lst (list): Input index list to partition.
        m (int): Number of groups required.

    Returns:
        list[list[list]]: All valid ordered partition schemes.
    """
    from itertools import product

    n = len(lst)
    if n < m:
        return []
    result = []
    for assign in product(range(m), repeat=n):
        if set(assign) == set(range(m)):
            groups = [[] for _ in range(m)]
            for val, g in zip(lst, assign):
                groups[g].append(val)
            result.append(groups)
    return result


def get_num_element_dict(composition):
    """
    Build a dictionary mapping stoichiometric count to element symbols.

    Args:
        composition (Composition): pymatgen Composition object.

    Returns:
        dict[int, list[str]]: {stoichiometric_number: [list_of_elements]}
    """
    data_dict = {}
    for k, v in composition.as_dict().items():
        v = int(v)
        if data_dict.get(v):
            data_dict[v].append(k)
        else:
            data_dict[v] = [k]
    return data_dict


def scale_structure(sturcture):
    """
    Rescale the volume of a candidate structure using the DLS volume predictor.
    Fall back to the original structure if prediction fails.

    Args:
        sturcture (Structure): Input crystal structure.

    Returns:
        Structure: Volume-scaled structure or original structure on exception.
    """
    predictor = DLSVolumePredictor(cutoff=5.0, min_scaling=0.5, max_scaling=1.5)
    try:
        struct_pred = predictor.get_predicted_structure(sturcture)
        return struct_pred
    except Exception:
        return sturcture
