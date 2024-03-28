from pymol import cmd, stored
from protacAlign_new import protacAlign_new, getAlignAtoms
from checkOverlapColor import checkOverlapColor
import re, os
from ternaryReduction_new import ternaryReduction_new
from findHBondNumber import findHBondNumber
from protac_bin import protac_bin, ternary_bin
from tqdm import tqdm
import time

'''
To use this script in PyMOL, make sure there are only PROTAC conformations in the workspace.
PyMOL> run autoProtac.py
PyMOL> autoProtac targetProteinObject, E3LigaseObject

Credit: Hanfeng Lin
Jin Wang Lab, Baylor College of Medicine
'''


@cmd.extend
def autoProtac(targetProtein, E3ligase, ternary_reduction=False):
    assert targetProtein in cmd.get_object_list(), \
        "no target protein detected in the workspace, check your spelling"
    assert E3ligase in cmd.get_object_list(), \
        "no E3 ligase protein detected in the workspace, check your spelling"

    start = time.perf_counter()
    print("Removing side chain...")
    cmd.remove("sidechain", quiet=0)
    print("Removing hydrogen...")
    cmd.remove("hydrogen", quiet=0)

    # Unifying chain identification and res_ID
    print("Renaming chain IDs and small molecules...")
    cmd.alter(targetProtein + " and organic", 'chain="X"')
    cmd.alter("organic", 'resi="1"')
    cmd.alter("*", 'segi=""')
    cmd.alter(targetProtein + " and polymer", 'chain="P"')
    cmd.alter(E3ligase + " and organic", 'chain="Y"')
    cmd.alter(E3ligase + " and polymer", 'chain="A"')

    # Get PROTAC object names
    print("Reading PROTAC information...")
    protac_objects = cmd.get_object_list()
    protac_objects.remove(targetProtein)
    protac_objects.remove(E3ligase)
    protac_name = re.match(r'^[A-Za-z0-9]*(?=_)',
                           protac_objects[0]).group()  # Rename protac res_name according to object name
    print(protac_name)
    assert len(protac_name) <= 4, \
        "Please ensure your conformation library pdb file is within 4 characters. e.g. 4172.pdb"
    print("Reformatting protac chain ID and name...")
    cmd.alter("organic and not (chain X or chain Y)", 'resn="' + protac_name + '"')
    cmd.alter("organic and not (chain X or chain Y)", 'chain="Z"')

    print('protac list: ' + str(protac_objects[0]) + " ~ " + str(protac_objects[-1]))

    align_atoms = getAlignAtoms()
    warhead_align_atoms = align_atoms[0]
    e3ligand_align_atoms = align_atoms[1]
    linker_warhead_align_atoms = align_atoms[2]
    linker_e3ligand_align_atoms = align_atoms[3]

    # Create conformation bins. The list contains filenames.pdb
    protac_bin_list = protac_bin(bin_size=300)

    # Open each bin for alignment
    current_ternary_count = 0

    for protac_bins in protac_bin_list:
        cmd.load(protac_bins)
        print("Loading " + protac_bins + " ...")

        # Do alignment for each end and remove crashed structures
        protac_objects = cmd.get_object_list()
        protac_objects.remove(targetProtein)
        protac_objects.remove(E3ligase)
        regex = re.compile(protac_name + r'_.*$')
        protac_objects = [i for i in protac_objects if regex.match(i)]
        for protac in protac_objects:
            protacAlign_new(targetProtein, E3ligase, protac, warhead_align_atoms, e3ligand_align_atoms,
                            linker_warhead_align_atoms, linker_e3ligand_align_atoms)
            for i in checkOverlapColor(targetProtein, E3ligase, protac):
                current_ternary_count += 1
                print("Productive Ternary Conformations detected so far: " + str(current_ternary_count))
        os.remove(protac_bins)
        ternary_bin(
            protac_bins[0:-4])  # save the ternary complexes detected in this protac bin. [0:-4] removes extension name

    # Recall all the successful ternary complexes from bins
    ternary_bin_list = []
    for file in os.listdir():
        if re.match(r"ternary_protac_bin_.*", file):
            ternary_bin_list.append(re.match(r"ternary_protac_bin_.*", file).group())
    for ternary_bins in ternary_bin_list:
        cmd.load(ternary_bins, format='pse', partial=1)  # partial=1: merge with current session
        print("Loading " + ternary_bins + " ...")
        os.remove(ternary_bins)

    success_list = cmd.get_object_list()
    success_list.remove(targetProtein)
    success_list.remove(E3ligase)
    print("Successful ternary complexes: " + str(len(success_list)) + "\n")
    align_finish = time.perf_counter()

    # Delete ternary complex if no hydrogen bond between proteins
    # findHBondNumber(success_list)

    success_list = cmd.get_object_list()
    success_list.remove(targetProtein)
    success_list.remove(E3ligase)

    if ternary_reduction:
        ternaryReduction_new(success_list)

    success_list = cmd.get_object_list()
    success_list.remove(targetProtein)
    success_list.remove(E3ligase)
    if success_list:
        cmd.extra_fit("chain A", success_list[0])

    cmd.order("*", "yes")
    reduction_finish = time.perf_counter()

    print("Non-redundant Ternary complexes: " + str(len(cmd.get_object_list()[2:])) + "\n")
    print(
        f"Alignment finished in {(align_finish - start) / 3600:0.3f} hours / {(align_finish - start) / 60:0.2f} minutes")
    print(
        f"Ternary reduction finished in {(reduction_finish - align_finish) / 3600:0.3f} hours / {(reduction_finish - align_finish) / 60:0.2f} minutes\n")
    print(
        f"Overall computation hours: {(reduction_finish - start) / 3600:0.3f} hours / {(reduction_finish - start) / 60:0.2f} minutes\n")
    print("Thanks for using autoProtac!")


cmd.auto_arg[0]['autoProtac'] = [cmd.object_sc, 'object', ',']
cmd.auto_arg[1]['autoProtac'] = [cmd.object_sc, 'object', '']
