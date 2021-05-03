def load_hmmer_profiles(hmmer_file):
    profiles = {}
    current_profile = ''
    with open(hmmer_file, 'r') as f:
        for line in f:
            current_profile += line
            if line.startswith('ACC   '):
                profile_name = line.split()[-1].strip('\n').split('.')[0]
            if line.startswith('//'):
                profiles[profile_name] = current_profile
                current_profile = ''
    return profiles


def get_all_marker_set_markers(markers_sets_file):
    all_markers = set()
    sets_dict = {}
    with open(markers_sets_file, 'r') as f:
        for line in f:
            line = line.strip('\n \t').split('\t')
            marker_sets = line[-1].replace(']), set([', ';')
            for char in '[](\'() \"':
                marker_sets = marker_sets.replace(char, '').replace('set', '')
            marker_sets = [[marker.split('.')[0] for marker in marker_set.split(',')] for marker_set in
                           marker_sets.split(';') if not marker_set.split(',') == ['']]
            sets_dict[line[2]] = marker_sets
            marker_sets = [marker for marker_set in marker_sets for marker in marker_set]
            all_markers.update(marker_sets)
    return all_markers, sets_dict


def tigrfam2pfam_dict(tigrfam2pfam_file):
    tigrfam2pfam_dict = {}
    with open(tigrfam2pfam_file, 'r') as f:
        for line in f:
            line = line.strip('\n \t').split('\t')
            tigrfam, pfam = line[1], line[0].replace('pfam', 'PF')
            if not tigrfam2pfam_dict.get(tigrfam):
                tigrfam2pfam_dict[tigrfam] = [pfam]
            else:
                tigrfam2pfam_dict[tigrfam].append(pfam)
            if not tigrfam2pfam_dict.get(pfam):
                tigrfam2pfam_dict[pfam] = [tigrfam]
            else:
                tigrfam2pfam_dict[pfam].append(tigrfam)
    return tigrfam2pfam_dict


def remove_unused_checkm_hmm_profiles(hmmer_file, markers_sets_file, tigrfam2pfam_file, out_file):
    unused = []
    profiles = load_hmmer_profiles(hmmer_file)
    profiles_accs = list(profiles.keys())
    all_markers, all_markers_dict = get_all_marker_set_markers(markers_sets_file)
    tf2pf = tigrfam2pfam_dict(tigrfam2pfam_file)
    with open(out_file, 'w') as of:
        for profile in profiles_accs:
            if profile in all_markers or any(alt_acc in all_markers for alt_acc in tf2pf.get(profile, [])):
                # if any(alt_acc in all_markers for alt_acc in tf2pf.get(profile, [])) and not profile in all_markers:
                    # print(profile, tf2pf.get(profile, []))
                of.write(profiles[profile])
            else:
                # print('not found', profile, tf2pf.get(profile, 'No alt found.'))
                unused.append([profile, tf2pf.get(profile, '')])
    return unused

hmmer_file = '/Users/oskar.hickl/Downloads/checkm_data_2015_01_16/hmms/checkm.hmm'
markers_sets_file = '/Users/oskar.hickl/Downloads/checkm_data_2015_01_16/taxon_marker_sets.tsv'
tigrfam2pfam_file = '/Users/oskar.hickl/Downloads/checkm_data_2015_01_16/pfam/tigrfam2pfam.tsv'
out_file = '/Users/oskar.hickl/Downloads/checkm_data_2015_01_16/hmms/checkm_binny_test.hmm'
old_binny_file = '/Users/oskar.hickl/Downloads/checkm_data_2015_01_16/hmms/checkm_binny.hmm'


unused_profiles = remove_unused_checkm_hmm_profiles(hmmer_file, markers_sets_file, tigrfam2pfam_file, out_file)

tf2pf = tigrfam2pfam_dict(tigrfam2pfam_file)


markers_set, markers_dict = get_all_marker_set_markers(markers_sets_file)
old_binny_hmms = set(list(load_hmmer_profiles(old_binny_file).keys()))
new_binny_hmms = set(list(load_hmmer_profiles(out_file).keys()))

not_in_old = new_binny_hmms - old_binny_hmms
not_in_new = old_binny_hmms - new_binny_hmms

example_domtblout = '/Users/oskar.hickl/Downloads/checkm_data_2015_01_16/hmms/prokka.faa.checkm.domtblout.v4.hmmscan'


def parse_domtblout(domtblout_file, i_eval_threshold=1e-10):
    # Gets non-overlapping hits with i-eval < than specified (or 1e-10 by default) in i-eval ascending order
    domtblout_data = {}
    with open(domtblout_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                 line = line.strip('\n \t').split()
                 gene, acc, i_eval, start, end = line[0], line[4].split('.')[0], float(line[12]), int(line[15]), int(line[16])
                 if i_eval <= i_eval_threshold:
                    if not domtblout_data.get(gene):
                        domtblout_data[gene] = [[acc, i_eval, start, end]]
                    else:
                        if (not any(acc == other_acc[0] for other_acc in domtblout_data[gene])
                                and not any(alt_acc == other_acc[0] for other_acc in domtblout_data[gene] for alt_acc in tf2pf.get(acc, []))):
                            # print('test1', [other_acc[0] for other_acc in domtblout_data[gene]])
                            # print('test2', acc, [alt_acc for alt_acc in tf2pf.get(acc, [])])
                            domtblout_data[gene].append([acc, i_eval, start, end])
                        elif any(acc == other_acc[0] for other_acc in domtblout_data[gene]):
                            for ind, i in enumerate(domtblout_data[gene]):
                                if acc == i[0]:
                                    if i_eval < i[1]:
                                        # print('Better hit for same ID {0} ({1}): Old i-eval(ID: {2}) = {3}, new i-eval(ID: {4}) = {5}'.format(
                                        #     acc, tf2pf.get(acc, []), i[0],  i[1], acc, i_eval))
                                        domtblout_data[gene][ind] = [acc, i_eval, start, end]
                        elif any(alt_acc == other_acc[0] for other_acc in domtblout_data[gene] for alt_acc in tf2pf.get(acc, [])):
                            for ind, i in enumerate(domtblout_data[gene]):
                                if any(alt_acc == i[0] for alt_acc in tf2pf.get(acc, [])):
                                    if i_eval < i[1]:
                                        # print('Better hit for same ID {0} ({1}): Old i-eval(ID: {2}) = {3}, new i-eval(ID: {4}) = {5}'.format(
                                        #     acc, tf2pf.get(acc, []), i[0], i[1], acc, i_eval))
                                        domtblout_data[gene][ind] = [acc, i_eval, start, end]
                        # domtblout_data[gene].append([acc, i_eval, start, end])
    for k in domtblout_data:
        domtblout_data[k].sort(key=lambda x: x[1], reverse=False)

    for k, v in domtblout_data.items():
        print(k, v)
        new_v = [v[0]]
        for ind, i in enumerate(v):
            if ind > 0:
                if not any(pos in range(v[index][2], v[index][3]) for pos in range(v[ind][2], v[ind][3]) for index in range(ind)):
                    new_v.append(i)
        domtblout_data[k] = new_v
        # if len(new_v) > 1:
        #     print(k, v)
        #     print(k,new_v)
    return domtblout_data


gene_markers = parse_domtblout(example_domtblout)
