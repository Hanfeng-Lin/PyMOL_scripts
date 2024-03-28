from pymol import cmd, stored


# backbone-only, no-hydrogen protein should be over 3.0 A.

def checkOverlapColor(sel1, sel2, protac):  # targetProtein, E3ligase, protac
    cutoff_distance = 2.6

    cutoffColorDictionary = {round(cutoff_distance + 1.2, 1): "0xff496c",
                             round(cutoff_distance + 1.0, 1): "0xffa343",
                             round(cutoff_distance + 0.8, 1): "0xffcf48",
                             round(cutoff_distance + 0.6, 1): "0xb2ec5d",
                             round(cutoff_distance + 0.4, 1): "0x45cea2",
                             round(cutoff_distance + 0.2, 1): "0x1974d2",
                             round(cutoff_distance, 1): "0x7851a9"}  # the minimal allowed distance between any two atoms
    for cutoff, color in cutoffColorDictionary.items():
        overlap = cmd.select(
            sel1 + " and polymer within " + str(cutoff) + " of " + sel2 + " and polymer")  # return overlap atom number
        if overlap == 0:
            protac_overlap = cmd.select(
                protac + " within " + str(cutoff) + " of ((" + sel1 + " or " + sel2 + ") and polymer)")
            if protac_overlap == 0:
                print("\n!!!!!!" + protac + ": no clashes detected within " + str(cutoff) + " Angstrom.")
                yield protac
                cmd.color(color, sel1 + " and elem C")
                cmd.create(str(cutoff) + "_" + protac + "ter", sel1 + " or " + sel2 + " or " + protac)
                break
            else:
                print(protac + ": detected " + str(overlap) + " backbone atoms overlapping within " + str(
                    cutoff) + " Angstrom, but detected " + str(
                    protac_overlap) + " atoms linker overlapping.")
        else:
            print(protac + ": detected " + str(overlap) + " backbone atoms overlapping within " + str(
                cutoff) + " Angstrom)")
    cmd.delete(protac)
