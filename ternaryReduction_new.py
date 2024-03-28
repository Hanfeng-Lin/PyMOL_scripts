from pymol import cmd, stored

"""
This new ternaryReduction method will initially try to classify a new incoming conformation into a cluster with 
similar rmsd members by searching against the first member of each cluster. If similar cluster doesn't exist, 
a new cluster will be created. Then rmsd search will be conducted starting from the cluster it belongs to, 
and moving along to other similar clusters. This increases the chance to find any similar conformation at an early time 
point, thus reducing the system complexity from pairwise search O(n(n-1)) 
"""


def ternaryReduction_new(success_list: list, rmsd_cutoff=0.5, cluster_cutoff=5, confident_cutoff=8):
    delete_set = set()

    # Clustering every single conformation according to their RMSD core
    cluster_list = [[success_list[0]]]  # A nested list indicating the conformation in different clusters
    for conf_i in success_list[1:len(success_list)]:
        print("\nChecking " + conf_i + "...")

        # Store the cluster searching rmsd score in a dict containing {cluster id : rmsd score}
        RMSD_dict = {}
        cluster_counter = -1
        for cluster in cluster_list:  # check rmsd against the first conformation in each cluster, and write into dict
            cluster_counter += 1
            rmsd = cmd.align(conf_i + " and (chain A or chain P)", cluster[0] +
                             " and (chain A or chain P)", cycles=0)[0]
            print(conf_i + " vs cluster_" + str(cluster_counter) + ": RMSD = " + str(rmsd))
            RMSD_dict[cluster_counter] = rmsd  # add rmsd score to the dict

        # Generate a list only containing cluster_id sorted by rmsd
        RMSD_sorted_list = sorted(RMSD_dict.items(), key=lambda x: x[1])  # a list containing (cluster_id, rmsd)
        cluster_id_sorted = []  # This list describes, for an object, the rmsd closest clusterID to the further clusterID
        for item in RMSD_sorted_list:
            cluster_id_sorted.append(item[0])

        # align conf_i against everything, starting from the lowest RMSD cluster members
        for cluster_id in cluster_id_sorted:
            for conf_j in cluster_list[cluster_id]:
                if conf_j in delete_set:
                    print(conf_j + " has already been in delete_list")
                    continue
                rmsd = cmd.align(conf_i + " and (chain A or chain P)", conf_j + " and (chain A or chain P)", cycles=0)[0]
                if rmsd >= confident_cutoff:
                    print(conf_i + " is different from " + conf_j + " in cluster_" + str(cluster_id) + ". RMSD = " + str(rmsd))
                    print(conf_i + " is believed to be different from all other conformations with a rmsd > " + str(confident_cutoff))
                    break
                elif rmsd >= rmsd_cutoff:
                    print(conf_i + " is different from " + conf_j + " in cluster_" + str(cluster_id) + ". RMSD = " + str(rmsd))
                else:
                    print(conf_i + " is similar with " + conf_j + " in cluster_" + str(cluster_id) + ". RMSD = " + str(rmsd) + ". \nRemoving " + conf_i)
                    delete_set.add(conf_i)
                    break
            else:  # only execute when it's no break in the inner loop
                continue
            break

        # Survived structures will be explicitly aligned with E3 ligase
        if conf_i not in delete_set:
            cmd.align(str(conf_i) + " and chain A", cluster_list[0][0] + " and chain A")
            # Add conf_i to the lowest rmsd cluster
            if any(rmsd <= cluster_cutoff for rmsd in RMSD_dict.values()):
                print("Appending " + conf_i + " to the cluster_" + str(cluster_id_sorted[0]))
                cluster_list[cluster_id_sorted[0]].append(conf_i)
            # If rmsd for all cluster > cluster_cutoff, start a new cluster
            if all(rmsd > cluster_cutoff for rmsd in RMSD_dict.values()):
                print("No cluster suitable for " + conf_i + ". Generating cluster_" + str(len(cluster_list)) + "...")
                cluster_list.append([conf_i])

    # delete redundant ternary complexes
    print("Objects to be deleted: " + str(delete_set))
    for ternary in delete_set:
        cmd.delete(ternary)
    print(str(len(delete_set)) + " objects has been deleted \n")

    # print cluster information
    cluster_number = -1
    for cluster in cluster_list:
        cluster_number += 1
        print("cluster_" + "%02d" % cluster_number + ": " + ' '.join(str(x) for x in cluster))
        for x in cluster:
            cmd.set_name(x, "clu" + "%02d" % cluster_number + "_" + x)
