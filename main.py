from sub import Substitute
import pandas as pd

sub = Substitute()
sub.get_templates_from_composition('BN2O', mode='wyckoff', atom_num_match=False, set_num_match=True)

