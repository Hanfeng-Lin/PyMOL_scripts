"""
Align warhead/E3_ligand align_atoms with linker align_atoms. If RMSD < 0.4, keep this linker conformer.
An additional minimization step is needed to deal with protein-linker clashes.

Usage: Put all the Rosetta docked protein with their ligands into the workspace. Also add the omega linker conformations
Make sure chain X and Y are warhead and E3 ligand.
"""
import os

from lxml import etree
from pymol import cmd, stored
from protac_bin import protac_bin, good_linker_bin
import re


# Retrieve warhead and E3 ligand aligning atom names. Linker resname is based on input.
def getAlignAtoms(linker):
    protac_xml = etree.parse("PROTAC.xml")
    stored.ligands = []  # This list will store three ligands: warhead, E3 ligand
    cmd.iterate("first chain X", "stored.ligands.append(resn)")  # returns warhead resn
    cmd.iterate("first chain Y", "stored.ligands.append(resn)")  # returns E3 ligand resn
    warhead = stored.ligands[0]
    e3ligand = stored.ligands[1]
    warhead_Align_atoms = protac_xml.xpath('//ligand[@title="' + warhead + '"]/align_atoms/atom/text()')
    e3ligand_Align_atoms = protac_xml.xpath('//ligand[@title="' + e3ligand + '"]/align_atoms/atom/text()')
    linker_warhead_Align_atoms = protac_xml.xpath('//linker[@title="' + linker + '"]/warhead_align_atoms/atom/text()')
    linker_e3ligand_Align_atoms = protac_xml.xpath('//linker[@title="' + linker + '"]/E3_align_atoms/atom/text()')
    assert warhead_Align_atoms, \
        "Unrecognized warhead. Please add the ligand in the target protein structure to the PROTAC.xml database"
    assert e3ligand_Align_atoms, \
        "Unrecognized E3 ligand. Please add the ligand in the E3 ligase structure to the PROTAC.xml database"
    assert linker_e3ligand_Align_atoms, \
        "Unrecognized linker. Please add this linker to the PROTAC.xml database. Or, did you forget to set the <E3_align_atoms> tag for this linker in the PROTAC.xml?"
    assert linker_warhead_Align_atoms, \
        "Unrecognized linker. Please add this linker to the PROTAC.xml database. Or, did you forget to set the <warhead_align_atoms> tag for this linker in the PROTAC.xml?"
    assert len(warhead_Align_atoms) == len(linker_warhead_Align_atoms)
    assert len(e3ligand_Align_atoms) == len(linker_e3ligand_Align_atoms)

    return warhead_Align_atoms, e3ligand_Align_atoms, linker_warhead_Align_atoms, linker_e3ligand_Align_atoms


@cmd.extend
def linker_builder(docked_object, rmsd_cutoff: float = 0.4):
    cmd.remove("hydrogen", quiet=0)
    # get PROTAC name (linker name) from the workspace
    first_protac_object = None
    for workspace_object in cmd.get_object_list():
        if cmd.select(workspace_object + " and polymer") == 0:
            first_protac_object = workspace_object
            break
    protac_name = re.match(r'^[A-Za-z0-9]*(?=_)', first_protac_object).group()
    print("PROTAC name is: " + protac_name)

    # Rename omega conformations
    print("Renaming conformations...")
    cmd.alter("organic and not (chain X or chain Y)", 'resn="' + protac_name + '"')
    cmd.alter("organic and not (chain X or chain Y)", 'chain="Z"')

    print("Inquiring align atoms from PROTAC.xml...")
    align_atoms = getAlignAtoms(protac_name)

    warhead_align_atoms = align_atoms[0]
    e3ligand_align_atoms = align_atoms[1]
    linker_warhead_align_atoms = align_atoms[2]
    linker_e3ligand_align_atoms = align_atoms[3]

    protac_bin_list = protac_bin(bin_size=200)
    good_linker_bin_list = []
    for protac_bins in protac_bin_list:
        cmd.load(protac_bins)
        print("Loading " + protac_bins + " ...")
        # Get current protac objects in the workspace
        protac_objects = []
        for workspace_object in cmd.get_object_list():
            if cmd.select(workspace_object + " and polymer") == 0:
                protac_objects.append(workspace_object)

        for protac in protac_objects:
            select_ligand_object_str = []
            select_linker_object_str = []

            # Linker alignment for both warhead and E3 ligand
            counter = 0
            for atom in warhead_align_atoms:
                counter += 1
                cmd.select("select_target_" + str(counter), docked_object + " and chain X and name " + atom)
                select_ligand_object_str.append("select_target_" + str(counter))
            for atom in e3ligand_align_atoms:
                counter += 1
                cmd.select("select_E3ligand_" + str(counter), docked_object + " and chain Y and name " + atom)
                select_ligand_object_str.append("select_E3ligand_" + str(counter))

            counter = 0
            for atom in linker_warhead_align_atoms:
                counter += 1
                cmd.select("select_warhead_" + str(counter), protac + " and chain Z and name " + atom)
                select_linker_object_str.append("select_warhead_" + str(counter))
            for atom in linker_e3ligand_align_atoms:
                counter += 1
                cmd.select("select_E3_" + str(counter), protac + " and chain Z and name " + atom)
                select_linker_object_str.append("select_E3_" + str(counter))

            pair_fit_list = []  # making the cmd.pair_fit *args
            for i in range(max(len(select_ligand_object_str), len(select_linker_object_str))):
                if select_linker_object_str:
                    pair_fit_list.append(select_linker_object_str.pop(0))
                if select_ligand_object_str:
                    pair_fit_list.append(select_ligand_object_str.pop(0))

            print("-----------------------------------")
            print(protac)
            rmsd = cmd.pair_fit(*pair_fit_list)
            if rmsd > float(rmsd_cutoff):
                cmd.delete(protac)
            else:
                print(protac + " fits well in " + docked_object)
                cmd.create(protac + "_"+f'{rmsd:.3f}', protac)
                cmd.delete(protac)
                continue
        # Temporarily save good linker so far
        good_linker_bin_list.append(good_linker_bin(protac_bins))

    for good_linker_pse in good_linker_bin_list:
        print("Loading " + good_linker_pse + "...")
        cmd.load(good_linker_pse, format='pse', partial=1)  # partial=1: merge with current session
        os.remove(good_linker_pse)
    cmd.delete("select*")
    protac_objects = []
    for workspace_object in cmd.get_object_list():
        if cmd.select(workspace_object + " and polymer") == 0:
            protac_objects.append(workspace_object)
    print("Good linkers for " + docked_object + ": " + str(len(protac_objects)))
    min_item = min(protac_objects, key=lambda x: float(x.split('_')[-1]))
    print("Best linker with the smallest RMSD: "+ str(min_item))
    # Clean up files
    for protac_bins in protac_bin_list:
        os.remove(protac_bins)


cmd.auto_arg[0]['linker_builder'] = [cmd.object_sc, 'object', '']
