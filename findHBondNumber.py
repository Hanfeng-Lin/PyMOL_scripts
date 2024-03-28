from pymol import cmd, stored
from get_raw_distances import get_raw_distances


def findHBondNumber(success_list: list):
    for ternary_object in success_list:
        cmd.dist("HBond", ternary_object + ' and chain A', ternary_object + ' and chain P', mode=2)
        H_bond_info = get_raw_distances("HBond")
        if len(H_bond_info) == 0:
            print(ternary_object + ' has 0 protein-protein H-bond, deleting ' + ternary_object + '...')
            cmd.delete(ternary_object)
            cmd.delete("HBond")
        else:
            print(ternary_object + ' has ' + str(len(H_bond_info)) + ' protein-protein H-bond.')
            cmd.delete("HBond")
    print("H-bond filtering completed.\n")