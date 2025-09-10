"""
Align warhead/E3_ligand align_atoms with linker align_atoms. If RMSD < 0.4, keep this linker conformer.
An additional minimization step is needed to deal with protein-linker clashes.

Usage: Put all the Rosetta docked structures with their ligands into the workspace. Also load the omega linker conformations.
Make sure chain X and Y are warhead and E3 ligand.

Then run the following pymol command: 

(dock all ternary objects in one session with default cutoff=0.4):
> linker_builder

(dock all ternary objects in one session):
> linker_builder rmsd_cutoff=0.4

(dock a specific ternary object in one session):
> linker_builder obj1 rmsd_cutoff=0.4
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
def linker_builder(docked_objects: str = "", rmsd_cutoff: float = 0.4):
    """
    Build linkers for one or multiple docked ternary objects.
    docked_objects: space-separated list of PyMOL object names. If empty, auto-detect ternary objects.
    """

    cmd.remove("hydrogen", quiet=0)

    # If user doesn't specify, find all ternary objects
    if not docked_objects:
        docked_objects = [
            obj for obj in cmd.get_object_list()
            if cmd.select(obj + " and chain A and polymer") > 0
            and cmd.select(obj + " and chain P and polymer") > 0
        ]
    else:
        docked_objects = docked_objects.split()

    print("Processing ternary objects: " + ", ".join(docked_objects))

    # Get PROTAC name (linker name) from the workspace
    first_protac_object = None
    for workspace_object in cmd.get_object_list():
        if cmd.select(workspace_object + " and polymer") == 0:
            first_protac_object = workspace_object
            break
    protac_name = re.match(r'^[A-Za-z0-9]*(?=_)', first_protac_object).group()
    print("PROTAC name is: " + protac_name)

    # Rename omega conformations
    cmd.alter("organic and not (chain X or chain Y)", f'resn="{protac_name}"')
    cmd.alter("organic and not (chain X or chain Y)", 'chain="Z"')

    # Inquire align atoms from PROTAC.xml
    align_atoms = getAlignAtoms(protac_name)
    warhead_align_atoms, e3ligand_align_atoms, linker_warhead_align_atoms, linker_e3ligand_align_atoms = align_atoms

    print("warhead_align_atoms:", warhead_align_atoms)
    print("e3ligand_align_atoms:", e3ligand_align_atoms)
    print("linker_warhead_align_atoms:", linker_warhead_align_atoms)
    print("linker_e3ligand_align_atoms:", linker_e3ligand_align_atoms)

    # Generate protac bins once for all ternaries
    protac_bin_list = protac_bin(bin_size=200)

    # --- Save each ternary object into its own pse bin ---
    ternary_bins = []
    for ternary in docked_objects:
        pse_name = f"{ternary}.pse"
        cmd.save(pse_name, ternary)
        ternary_bins.append(pse_name)

    # Track successful merged ternaries
    success_pse_list = []

    # Process each ternary separately
    for ternary_bin in ternary_bins:
        ternary_name = os.path.splitext(os.path.basename(ternary_bin))[0]
        print("\n=== Working on ternary object:", ternary_name, "===")

        # Load just this ternary
        cmd.load(ternary_bin)
        docked_object = ternary_name

        good_linker_bin_list = []

        # Reload protac conformers for this ternary
        for protac_bins in protac_bin_list:
            cmd.load(protac_bins)
            print("Loading " + protac_bins + " ...")

            protac_objects = [
                obj for obj in cmd.get_object_list()
                if cmd.select(obj + " and polymer") == 0
            ]

            for protac in protac_objects:
                select_ligand_object_str = []
                select_linker_object_str = []

                counter = 0
                for atom in warhead_align_atoms:
                    counter += 1
                    cmd.select(f"select_target_{counter}", docked_object + " and chain X and name " + atom)
                    select_ligand_object_str.append(f"select_target_{counter}")
                for atom in e3ligand_align_atoms:
                    counter += 1
                    cmd.select(f"select_E3ligand_{counter}", docked_object + " and chain Y and name " + atom)
                    select_ligand_object_str.append(f"select_E3ligand_{counter}")

                counter = 0
                for atom in linker_warhead_align_atoms:
                    counter += 1
                    cmd.select(f"select_warhead_{counter}", protac + " and chain Z and name " + atom)
                    select_linker_object_str.append(f"select_warhead_{counter}")
                for atom in linker_e3ligand_align_atoms:
                    counter += 1
                    cmd.select(f"select_E3_{counter}", protac + " and chain Z and name " + atom)
                    select_linker_object_str.append(f"select_E3_{counter}")

                # Pair fit linker to ternary
                pair_fit_list = []
                for i in range(max(len(select_ligand_object_str), len(select_linker_object_str))):
                    if select_linker_object_str:
                        pair_fit_list.append(select_linker_object_str.pop(0))
                    if select_ligand_object_str:
                        pair_fit_list.append(select_ligand_object_str.pop(0))

                rmsd = cmd.pair_fit(*pair_fit_list)
                if rmsd > float(rmsd_cutoff):
                    cmd.delete(protac)
                else:
                    print(f"{protac} fits well in {docked_object} (RMSD={rmsd:.3f})")
                    cmd.create(protac + "_" + f"{rmsd:.3f}", protac)
                    cmd.delete(protac)

            good_linker_bin_list.append(good_linker_bin(protac_bins))

        # Reload all good linkers for this ternary
        for good_linker_pse in good_linker_bin_list:
            print("Loading " + good_linker_pse + "...")
            cmd.load(good_linker_pse, format='pse', partial=1)
            os.remove(good_linker_pse)

        cmd.delete("select*")
        protac_objects = [
            obj for obj in cmd.get_object_list()
            if cmd.select(obj + " and polymer") == 0
        ]

        if not protac_objects:
            print("No good linkers found for " + docked_object)
        else:
            # Select best linker (lowest RMSD)
            min_item = min(protac_objects, key=lambda x: float(x.split('_')[-1]))
            print("Best linker for " + docked_object + ": " + str(min_item))

            # Merge best linker with ternary
            merged_name = docked_object + "_with_linker_rmsd_" + min_item.split('_')[-1]
            cmd.create(merged_name, docked_object + " or " + min_item)
            print("Merged best linker into " + merged_name)

            # Save successful ternary to pse
            merged_pse = merged_name + ".pse"
            cmd.save(merged_pse, merged_name)
            success_pse_list.append(merged_pse)
            print("Saved successful ternary as " + merged_pse)

        # Clean up this ternary before moving on
        cmd.delete(docked_object)

    # Cleanup protac bins and ternary bins
    for protac_bins in protac_bin_list:
        os.remove(protac_bins)
    for ternary_bin in ternary_bins:
        os.remove(ternary_bin)

    # Reload successful ternaries at the end
    print("\n=== Reloading successful ternaries into session ===")
    for merged_pse in success_pse_list:
        cmd.load(merged_pse, partial=1)
        print("Loaded " + merged_pse)
        os.remove(merged_pse)



cmd.auto_arg[0]['linker_builder'] = [cmd.object_sc, 'object', '']
