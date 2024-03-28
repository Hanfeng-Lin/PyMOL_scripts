import os
import re
from pymol import cmd

"""
Save every 300 protac conformations to a bin
"""


def protac_bin(bin_size=300):
    # Identify protac objects
    protac_objects = []
    first_protac_object = None
    for workspace_object in cmd.get_object_list():
        if cmd.select(workspace_object + " and polymer") == 0:
            first_protac_object = workspace_object
            break
    protac_name = re.match(r'.*(?=_)', first_protac_object).group()
    for workspace_object in cmd.get_object_list():
        if protac_name in workspace_object:
            protac_objects.append(workspace_object)

    bin_total = len(protac_objects) // bin_size
    if len(protac_objects) % bin_size != 0:
        bin_total = bin_total + 1
    print("Conformation number: " + str(len(protac_objects)) +
          "\nA total of " + str(bin_total) + " bins will be generated.")
    for bin_number in range(bin_total):
        cmd.select("none")
        bin_count = 0
        for protac in protac_objects:
            cmd.select("sele or " + protac)
            bin_count += 1
            if bin_count == bin_size:
                break
        print("Saving protac_bin_%02d.pdb ..." % bin_number)
        cmd.multisave("protac_bin_%02d.pdb" % bin_number, "sele", format="pdb")
        for object_name in cmd.get_names('objects', 0, 'sele'):
            cmd.delete(object_name)
        # Remove the saved protac from the list
        protac_objects = protac_objects[bin_size:]

    protac_bin_list = []
    for file in os.listdir():
        if re.match(r'protac_bin_.*\.pdb$', file):
            protac_bin_list.append(re.match(r'protac_bin_.*\.pdb$', file).group())

    return protac_bin_list


def ternary_bin(protac_bins):
    ternary_objects = cmd.get_object_list()
    regex = re.compile(r'.*_\d+ter$')
    ternary_objects = [i for i in ternary_objects if regex.match(i)]
    print("ternary objects:" + str(ternary_objects))
    print("Saving ternary_" + protac_bins + ".pse ...")
    cmd.save("ternary_" + protac_bins + ".pse", format="pse")
    for object_name in ternary_objects:
        cmd.delete(object_name)


def good_linker_bin(protac_bins):
    protac_objects = []
    for workspace_object in cmd.get_object_list():
        if cmd.select(workspace_object + " and polymer") == 0:
            protac_objects.append(workspace_object)
    cmd.select("none")
    for protac in protac_objects:
        cmd.select("sele or " + protac)
    print("Saving good_" + protac_bins[0:-4] + ".pse" + " ...")
    cmd.save("good_" + protac_bins[0:-4] + ".pse", format="pse")
    # clean up after saving
    for object_name in cmd.get_names('objects', 0, 'sele'):
        print("Clean up after saving...")
        cmd.delete(object_name)

    return "good_" + protac_bins[0:-4] + ".pse"
