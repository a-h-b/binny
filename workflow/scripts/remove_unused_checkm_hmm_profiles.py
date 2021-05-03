#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 6 08:43:36 2021

@author: oskar.hickl
"""

import sys
from pathlib import Path


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
    with open(markers_sets_file, 'r') as f:
        for line in f:
            line = line.strip('\n \t').split('\t')
            marker_sets = line[-1].replace(']), set([', ';')
            for char in '[](\'() \"':
                marker_sets = marker_sets.replace(char, '').replace('set', '')
            marker_sets = [[marker.split('.')[0] for marker in marker_set.split(',')] for marker_set in
                           marker_sets.split(';') if not marker_set.split(',') == ['']]
            marker_sets = [marker for marker_set in marker_sets for marker in marker_set]
            all_markers.update(marker_sets)
    return all_markers


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


def remove_unused_checkm_hmm_profiles(hmmer_file, markers_sets_file, tigrfam2pfam_file, out_path):
    profiles = load_hmmer_profiles(hmmer_file)
    profiles_accs = list(profiles.keys())
    all_markers = get_all_marker_set_markers(markers_sets_file)
    tf2pf = tigrfam2pfam_dict(tigrfam2pfam_file)
    out_path_pfam = Path(out_path + '/checkm_pf/checkm_filtered_pf.hmm')
    out_path_tigrfam = Path(out_path + '/checkm_tf/checkm_filtered_tf.hmm')
    out_path_pfam.parent.mkdir(parents=True, exist_ok=True)
    out_path_tigrfam.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path_pfam, 'w') as of_pfam:
        with open(out_path_tigrfam, 'w') as of_tigrfam:
            for profile in profiles_accs:
                if profile in all_markers or any(alt_acc in all_markers for alt_acc in tf2pf.get(profile, [])):
                    # if any(alt_acc in all_markers for alt_acc in tf2pf.get(profile, [])) and not profile in all_markers:
                        # print(profile, tf2pf.get(profile, []))
                    if profile.startswith('PF'):
                        of_pfam.write(profiles[profile])
                    elif profile.startswith('TIGR'):
                        of_tigrfam.write(profiles[profile])
                # else:
                #     print('not found', profile, tf2pf.get(profile, 'No alt found.'))


#hmmer_file = sys.argv[1]
#markers_sets_file = sys.argv[2]
#tigrfam2pfam_file = sys.argv[3]
#out_file = sys.argv[4]


#remove_unused_checkm_hmm_profiles(hmmer_file, markers_sets_file, tigrfam2pfam_file, out_file)
