from lxml import etree
from pymol import cmd, stored


def overlapping_atoms(ligand_name):
    # Read xml from file
    protac_xml = etree.parse("PROTAC.xml")
    # print(etree.tostring(PROTACxml))

    # xpath '//' search everything, returns with an Element object (which is a list)
    # with a '/text()' at the end of xpath can generate a list of node text
    overlapping_atom_list = protac_xml.xpath('//ligand[@title="' + ligand_name + '"]/overlapping_atoms/atom/text()')
    # ligand_type = protac_xml.xpath('//ligand[@title="' + ligand_name + '"]')[0].getparent().getparent().tag

    return overlapping_atom_list


def useless_linker_atoms(ligand_name):
    # Read xml from file
    protac_xml = etree.parse("PROTAC.xml")
    useless_atom_list = protac_xml.xpath(
        '//linker[@title="' + ligand_name + '"]/useless_atoms/atom/text()')  # How to make this a list!!!

    return useless_atom_list


@cmd.extend
def merge_PROTAC():
    cmd.remove("hydrogen", quiet=0)
    for workspace_object in cmd.get_object_list():
        stored.chains = set()
        cmd.iterate(workspace_object + " and not metal", "stored.chains.add(chain)")
        if len(stored.chains) == 5:
            ternary_object = workspace_object

            # The first part will delete overlapping atoms

            # get ligand residue name from pymol
            stored.ligands = []  # This list will store three ligands: warhead, E3 ligand, and linker.
            cmd.iterate("first chain X", "stored.ligands.append(resn)")  # returns resn
            cmd.iterate("first chain Y", "stored.ligands.append(resn)")
            cmd.iterate("first chain Z", "stored.ligands.append(resn)")
            warhead = stored.ligands[0]
            e3_ligand = stored.ligands[1]
            linker = stored.ligands[2]
            if overlapping_atoms(warhead):
                print("deleting " + warhead + " " + str(overlapping_atoms(warhead)))
                for atom in overlapping_atoms(warhead):
                    select_atom_string = ternary_object + " and chain X and name " + atom
                    cmd.remove(select_atom_string)

            if overlapping_atoms(e3_ligand):
                print("deleting " + e3_ligand + " " + str(overlapping_atoms(e3_ligand)))
                for atom in overlapping_atoms(e3_ligand):
                    select_atom_string = ternary_object + " and chain Y and name " + atom
                    cmd.remove(select_atom_string)

            linker_useless_atoms: list = useless_linker_atoms(linker)
            if linker_useless_atoms:
                print("deleting " + linker + " " + str(useless_linker_atoms(linker)))
                for atom in linker_useless_atoms:
                    select_atom_string = ternary_object + " and chain Z and name " + atom
                    cmd.remove(select_atom_string)

            # The second part will bond bridging points

            protac_xml = etree.parse("PROTAC.xml")

            warhead_bridge_point: list = protac_xml.xpath('//ligand[@title="' + warhead + '"]/bridge_point/text()')
            e3ligand_bridge_point: list = protac_xml.xpath('//ligand[@title="' + e3_ligand + '"]/bridge_point/text()')
            linker_bridge_point_warhead: list = protac_xml.xpath(
                '//linker[@title="' + linker + '"]/bridge_point_warhead/text()')
            linker_bridge_point_e3ligand: list = protac_xml.xpath(
                '//linker[@title="' + linker + '"]/bridge_point_E3ligand/text()')

            for atom1 in warhead_bridge_point:
                for atom2 in linker_bridge_point_warhead:
                    select_atom1_string = ternary_object + " and chain X and name " + atom1
                    select_atom2_string = ternary_object + " and chain Z and name " + atom2
                    print("Bridging \n" + select_atom1_string + "\n" + select_atom2_string)
                    cmd.bond(select_atom1_string, select_atom2_string)

            for atom1 in e3ligand_bridge_point:
                for atom2 in linker_bridge_point_e3ligand:
                    select_atom1_string = ternary_object + " and chain Y and name " + atom1
                    select_atom2_string = ternary_object + " and chain Z and name " + atom2
                    print("Bridging \n" + select_atom1_string + "\n" + select_atom2_string)
                    cmd.bond(select_atom1_string, select_atom2_string)

            # The third part will merge three items into one

            cmd.alter(ternary_object + " and resn " + warhead, 'resi="1"')
            cmd.alter(ternary_object + " and resn " + e3_ligand, 'resi="1"')
            cmd.alter(ternary_object + " and resn " + linker, 'resi="1"')
            cmd.alter(ternary_object + " and resn " + warhead, 'segi=""')
            cmd.alter(ternary_object + " and resn " + e3_ligand, 'segi=""')
            cmd.alter(ternary_object + " and resn " + linker, 'segi=""')
            cmd.alter(ternary_object + " and resn " + warhead, 'chain="Z"')
            cmd.alter(ternary_object + " and resn " + e3_ligand, 'chain="Z"')
            cmd.alter(ternary_object + " and resn " + warhead, 'resn="' + linker + '"')
            cmd.alter(ternary_object + " and resn " + e3_ligand, 'resn="' + linker + '"')

            cmd.rename()

        else:
            print(workspace_object + " is not a ternary complex, omitted.")
