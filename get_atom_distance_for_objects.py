from pymol import cmd, stored
from get_raw_distances import get_raw_distances
from protac_bin import protac_bin
import matplotlib.pyplot as plt
import os, math
"""
This scripts helps export distance objects in pymol into a csv file in a state-wise manner.
Usage: get_trajectory_distance("measure01", "measure02", ...)
csv file will be found in the pymol current working directory
"""


def get_atom_distance_for_objects(atom1, atom2):
    cmd.delete("dist_*")
    distance_list = []

    # Create conformation bins. The list contains filenames.pdb
    protac_bin_list = protac_bin(bin_size=300)

    # Open each bin for alignment
    current_ternary_count = 0

    for protac_bins in protac_bin_list:
        cmd.load(protac_bins)
        print("Loading " + protac_bins + " ...")
        for object in cmd.get_object_list():
            cmd.distance('dist_'+object, object+' and name '+atom1, object+' and name '+atom2)
            dist = get_raw_distances('dist_'+object)[0][2]
            distance_list.append(float(f'{dist:.2f}'))
            print("Calculating distance for "+object)
        os.remove(protac_bins)
        cmd.delete('*')

    distance_list.sort()
    print(distance_list)

    # Plot distance distribution

    interval_width = 0.1

    # Determine the minimum and maximum values in the data
    min_value = min(distance_list)
    max_value = max(distance_list)
    min_value_rounded = math.floor(min_value / interval_width) * interval_width
    max_value_rounded = math.ceil(max_value / interval_width) * interval_width

    # Create the intervals
    intervals = []
    for i in range(int((max_value_rounded - min_value_rounded) / interval_width) + 1):
        intervals.append(
            (round(min_value_rounded + i * interval_width, 1), round(min_value_rounded + (i + 1) * interval_width, 1)))

    # Count the frequency of numbers falling within each interval
    frequency_count = [sum(1 for num in distance_list if start <= num < end) for start, end in intervals]

    # Extract the interval labels for plotting (using the minimum value of each interval)
    interval_labels = [f'{start:.1f}' for start, _ in intervals]

    # Plotting
    plt.bar(range(len(intervals)), frequency_count, align='edge', width=1)
    plt.xlabel('Intervals')
    plt.ylabel('Frequency')
    plt.title('Frequency Plot with Intervals')

    # Set x-ticks positions and labels
    plt.xticks(range(len(intervals) + 1), interval_labels + [f'{intervals[-1][1]:.1f}'])

    plt.tight_layout()
    plt.show()


cmd.extend('get_atom_distance_for_objects', get_atom_distance_for_objects)
