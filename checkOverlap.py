from pymol import cmd, stored


def checkOverlap(sel1, sel2, protac):
    cutoff = 2.0  # the minimal allowed distance between any two atoms
    overlap = cmd.select(sel1 + " and polymer within " + str(cutoff) + " of " +
                         sel2 + " and polymer")  # return overlap atom number
    if overlap == 0:
        protacOverlap = cmd.select(protac + " within " + str(cutoff) + " of ((" +
                                   sel1 + " or " + sel2 + ") and polymer)")
        if protacOverlap == 0:
            print("\n!!!!!!" + protac + ": no clashes detected within " + str(cutoff) + " Angstrom.")
            yield protac
            cmd.create("ternary_" + protac, sel1 + " or " + sel2 + " or " + protac)
            cmd.delete(protac)
        else:
            print(protac + ": detected " + str(overlap) + " backbone atoms overlapping, but detected " + str(
                protacOverlap) + " atoms linker overlapping.")
            cmd.delete(protac)
    else:
        print(protac + ": detected " + str(overlap) + " backbone atoms overlapping. (Cutoff = " +
              str(cutoff) + " Angstrom)")
        cmd.delete(protac)


