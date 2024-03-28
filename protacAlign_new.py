from lxml import etree
from pymol import cmd, stored


def getAlignAtoms():
    protac_xml = etree.parse("PROTAC.xml")
    stored.ligands = []  # This list will store three ligands: warhead, E3 ligand, and linker.
    cmd.iterate("first chain X", "stored.ligands.append(resn)")  # returns resn
    cmd.iterate("first chain Y", "stored.ligands.append(resn)")
    cmd.iterate("first chain Z", "stored.ligands.append(resn)")
    warhead = stored.ligands[0]
    e3ligand = stored.ligands[1]
    linker = stored.ligands[2]
    print(f'warhead name {warhead} \n'
          f'e3ligand name {e3ligand} \n'
          f'linker name {linker}')
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


def protacAlign_new(targetProtein, E3ligase, protac, warhead_align_atoms, e3ligand_align_atoms,
                    linker_warhead_align_atoms, linker_e3ligand_align_atoms: list):

    select_target_object_str = []
    select_linker_warhead_object_str = []
    select_e3ligand_object_str = []
    select_e3_object_str = []

    # warhead alignment

    counter = 0
    for atom in warhead_align_atoms:
        counter += 1
        cmd.select("select_target_" + str(counter), targetProtein + " and chain X and name " + atom)
        select_target_object_str.append("select_target_" + str(counter))

    counter = 0
    for atom in linker_warhead_align_atoms:
        counter += 1
        cmd.select("select_warhead_" + str(counter), protac + " and chain Z and name " + atom)
        select_linker_warhead_object_str.append("select_warhead_" + str(counter))

    pair_fit_list = []  # making the cmd.pair_fit *args
    for i in range(max(len(select_target_object_str), len(select_linker_warhead_object_str))):
        if select_target_object_str:
            pair_fit_list.append(select_target_object_str.pop(0))
        if select_linker_warhead_object_str:
            pair_fit_list.append(select_linker_warhead_object_str.pop(0))

    cmd.pair_fit(*pair_fit_list)

    # E3 alignment

    counter = 0
    for atom in e3ligand_align_atoms:
        counter += 1
        cmd.select("select_E3ligand_" + str(counter), E3ligase + " and chain Y and name " + atom)
        select_e3ligand_object_str.append("select_E3ligand_" + str(counter))

    counter = 0
    for atom in linker_e3ligand_align_atoms:
        counter += 1
        cmd.select("select_E3_" + str(counter), protac + " and chain Z and name " + atom)
        select_e3_object_str.append("select_E3_" + str(counter))

    pair_fit_list = []  # making the cmd.pair_fit *args
    for i in range(max(len(select_e3ligand_object_str), len(select_e3_object_str))):
        if select_e3ligand_object_str:
            pair_fit_list.append(select_e3ligand_object_str.pop(0))
        if select_e3_object_str:
            pair_fit_list.append(select_e3_object_str.pop(0))

    cmd.pair_fit(*pair_fit_list)
    cmd.delete("select*")


