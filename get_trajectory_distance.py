from pymol import cmd, stored
from get_raw_distances import get_raw_distances
import pandas as pd

"""
This scripts helps export distance objects in pymol into a csv file in a state-wise manner.
Usage: get_trajectory_distance("measure01", "measure02", ...)
csv file will be found in the pymol current working directory
"""


def get_trajectory_distance(*distance_objects):
    total_state = cmd.count_states()
    state_list = range(1, total_state + 1)
    csv_dict = {'state': state_list}

    for dist_object in distance_objects:
        distance_list = []

        for state in range(total_state):
            try:
                current_distance = get_raw_distances(dist_object, state=state + 1)[0][2]
                # print distance from tuple(obj1,obj2,dist)
            except:
                current_distance = "NA"
            print("reading state " + str(state + 1) + " in " + dist_object + ": " + str(current_distance) + " A.")
            distance_list.append(current_distance)

        csv_dict[dist_object] = distance_list

    df = pd.DataFrame(csv_dict)
    df.to_csv('distance.csv', index=False)


cmd.extend('get_trajectory_distance', get_trajectory_distance)
