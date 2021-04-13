# embedr = '/Users/oskar.hickl/local_tools/EMBEDR'
#
# sys.path.append(embedr)
# from embedr import EMBEDR

########################################################################################################################
# checkm_hmmer_search2prokka_gff(prokka_checkm_marker_hmm_out, raw_annot)
# annot_df_checkm, annot_dict_checkm = gff2ess_gene_df(annot_file_checkm, target_attribute='checkm_marker', get_dict=True)
# test_bin = [marker.split('.')[0] for marker in annot_df_checkm[annot_df_checkm['contig'] == 'S5C9352']['essential'].tolist()[0].split(',')]
# best_bin_marker_set = chose_checkm_marker_set(test_bin, taxon_marker_sets, tigrfam2pfam_data)
########################################################################################################################


# # Standard scaler
# scaler = StandardScaler().fit(X)
# X_scaled = scaler.transform(X)

# def replaceZeroes(data):
#   min_nonzero = np.min(data[np.nonzero(data)])
#   data[data == 0] = min_nonzero * 0.0001
#   return data


########################################################################################################################
# mys = 'GGGATATCGATCGATCGATTATATATGACGTCGCGCGATTAGTCGCGCCCGTGATCCGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCC'
#
# def gc_skew(seq):
#     g = seq.count('G')
#     c = seq.count('C')
#     return (g - c) / (g + c)
#
# gc_skew_list = []
# counter = 0
# n_window = 1
# last_window = len(mys) // 10 + (len(mys) % 10 > 0)
# current_window = ''
# for nt in mys:
#     if counter < 10:
#         current_window += nt
#         counter += 1
#     elif counter == 10:
#         gc_skew_list.append(gc_skew(current_window))
#         current_window = nt
#         counter = 1
#         n_window += 1
#     elif n_window == last_window:
#         gc_skew_list.append(gc_skew(current_window))
########################################################################################################################


# print('Using n_neighbors of {0}.'.format(perp))
# X_embedded = umap.UMAP(n_neighbors=1000, min_dist=0.01, n_components=n_dim, random_state=0, densmap=True, dens_lambda=10,
#                        verbose=True , metric='cosine').fit_transform(X_scaled)  #
# data_Y = [X_embedded]


# use_densemap = True
# umap_param_dict = {'n_neighbors': perp, 'min_dist': 0.0, 'densmap': use_densemap,
#                    'dens_lambda': 0.1, 'metric': 'cosine'}  # , 'n_components': n_dim, 'random_state': 0
# tsne_param_dict = {'n_iter': 750, 'exaggeration': 12, 'early_exag_iter': 250, 'exag_mom': 0.5, 'momentum': 0.8,
#                    'metric': 'cosine', 'initialization': 'pca',
#                    'openTSNE_params': {'learning_rate': 'auto'}}
# embedr_obj = EMBEDR(perplexity=perp,
#                     random_state=0,
#                     dimred_alg='umap',
#                     n_components=n_dim,
#                     dimred_params=umap_param_dict,
#                     # project_dir='/Users/oskar.hickl/local_tools/EMBEDR/Embeddings',
#                     n_data_embed=1,
#                     n_null_embed=5,
#                     n_jobs=threads,
#                     cache_results=False,
#                     verbose=0)
# data_Y = embedr_obj.fit_transform(X_scaled)


# bad_points = 0
# embed_thresh = 0.05
# for i in embedr_obj.pVals:
#     if i > embed_thresh:
#         bad_points += 1
# print('Percentage of low quality embeddings above p-val threshold of {0}: {1}.'.format(embed_thresh, bad_points / len(
#     embedr_obj.pVals)))
#
# print('cont: {0}, kmer: {1}, n_dim: {2}, dense: {3}, embed: {4}, n_nei: {5}, %_bad: {6}.'.format(min_cont_size,
#                                                                                                ','.join([str(i) for i in kmer_sizes]),
#                                                                                                n_dim,
#                                                                                                use_densemap,
#                                                                                                embedding_tries,
#                                                                                                perp,
#                                                                                                round(bad_points / len(embedr_obj.pVals), 3)))

# if n_dim == 2:
#     fig, ax = plt.subplots(1, 1)
#     ax = embedr_obj.plot(ax=ax)
#     ax.title.set_text('cont:{0},kmer:{1},dim:{2},dense:{3},embed:{4},n_nei:{5},bad:{6}.'.format(min_cont_size,
#                                                                                                ','.join([str(i) for i in kmer_sizes]),
#                                                                                                n_dim,
#                                                                                                use_densemap,
#                                                                                                embedding_tries,
#                                                                                                perp,
#                                                                                                round(bad_points / len(embedr_obj.pVals), 3)))
# elif n_dim == 3:
#     fig = plt.figure()
#     ax = Axes3D(fig)
#     ax = embedr_obj.plot(ax=ax)
# plt.show()


# elif n_dim == 3:
#     fig = plt.figure()
#     ax = Axes3D(fig)
#     ax.embedr_obj.plot()
#     plt.show()

# print(len(embedr_obj.pVals))

# good_points_list_indices = [index for index, i in enumerate(embedr_obj.pVals) if i <= 1.00]  # embed_thresh
# ###
# # bad_points_list_indices = [index for index, i in enumerate(embedr_obj.pVals) if i > 0.1]
# # round_X_bad_contigs = [X_contigs[index] for index in bad_points_list_indices]
# ###
# if not good_points_list_indices:
#     break
# round_X_contigs_good = [round_X_contigs[index] for index in good_points_list_indices]
# round_X_good = [round_X[index] for index in good_points_list_indices]
#
# round_X_bad_contigs = [round_X_contigs[index] for index, i in enumerate(round_X_contigs) if index not in good_points_list_indices]
# round_X_bad = [round_X[index] for index, i in enumerate(round_X) if index not in good_points_list_indices]
# X_embedded_bad = [data_Y[0][index] for index, i in enumerate(data_Y[0]) if index not in good_points_list_indices]
#
# dim_range = [i + 1 for i in range(n_dim)]
# bad_coord_df = pd.DataFrame(data=X_embedded_bad, index=None, columns=['dim' + str(i) for i in dim_range])
# bad_coord_df['contig'] = round_X_bad_contigs
# # Reorder
# bad_coord_df = bad_coord_df[['contig'] + ['dim' + str(i) for i in dim_range]]
# bad_coords_file = 'intermediary/bad_contig_coordinates_dynamic_marker_sets_v04.tsv'
# bad_coord_df.to_csv(bad_coords_file, sep='\t', index=False, header=False)
# bad_contig_data_df = load_and_merge_cont_data(annot_file, mg_depth_file, bad_coords_file, dims=n_dim)
#
# if len(round_X_bad_contigs) != len(round_X_bad):
#     print('fuuuuuuuuck')
#     raise Exception

# X_scaled = clr(round_X)
#
# embedr_obj = EMBEDR(random_state=0,
#                     dimred_alg='umap',
#                     n_components=n_dim,
#                     dimred_params=umap_param_dict,
#                     # project_dir='/Users/oskar.hickl/local_tools/EMBEDR/Embeddings',
#                     n_data_embed=3,
#                     n_null_embed=5,
#                     n_jobs=threads,
#                     cache_results=False,
#                     verbose=1)
# data_Y = embedr_obj.fit_transform(X_scaled)
#
# ax1 = embedr_obj.plot()
# plt.show()
#
# print(len(embedr_obj.pVals))
# bad_points = 0
# for i in embedr_obj.pVals:
#     if i >= 0.05:
#         bad_points += 1
# print(bad_points / len(embedr_obj.pVals))
#
# good_points_list_indices = [index for index, i in enumerate(embedr_obj.pVals) if i <= 0.05]
# if not good_points_list_indices:
#     break
# round_X_contigs = [round_X_contigs[index] for index in good_points_list_indices]
# round_X = [round_X[index] for index in good_points_list_indices]

# if len(round_X_contigs_good ) != len(round_X_good ):
#     print('fuuuuuuuuck')
#     raise Exception

# X_embedded = [data_Y[0][index] for index in good_points_list_indices]


