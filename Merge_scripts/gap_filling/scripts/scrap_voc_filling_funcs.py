def mask_overlap(*vars_list):
    mask = vars_list[0].notna()
    for s in vars_list[1:]:
        mask &= s.notna()

    vars_overlap = [s.loc[mask] for s in vars_list]
def fitting_metrics(x_overlap, y_overlap):

    slope, intercept = np.polyfit(x_overlap, y_overlap, 1)
    rmse = np.sqrt(np.mean((x_overlap - y_overlap)**2))
    
    #R squared calculation
    r = np.corrcoef(x_overlap, y_overlap)[0, 1]
    r2 = r**2

    return slope, intercept, rmse, r, r2
def fitting_metrics_revised(x_overlap, y_overlap):

    slope, intercept = np.polyfit(x_overlap, y_overlap, 1)
    rmse = np.sqrt(np.mean((x_overlap - y_overlap)**2))
    
    #R squared calculation
    r = np.corrcoef(x_overlap, y_overlap)[0, 1]
    r2 = r**2

    return slope, intercept, rmse, r, r2
def odr_improves_revised(metrics_init_input, metrics_corr_input, y_overlap):
def process_one_voc_stage1_revised(var_species_name, duplicate_species_name, var_species_data, duplicate_species_data):
    cache_file_stage1 = f'{CACHE_DIR}/{var_species_name}_stage1.pkl'
    
    # ---- load if already processed ----
    # if os.path.exists(cache_file_stage1):
    #     with open(cache_file_stage1, 'rb') as f:
    #         return pickle.load(f)

        # If any of the results need to be modified, comment out the return pickle.load(f) line and replace with code below with alterations. 
        # This one changes a key previously named 'metrics_ml_udaq_corrected' into 'metrics_ml_udaq_corr'
        #     results = pickle.load(f)
        # print('results: ', results)
        # if 'metrics_ml_udaq_corrected' in results:
        #     results['metrics_ml_udaq_corr'] = results.pop('metrics_ml_udaq_corrected')
        # with open(cache_file, 'wb') as f:
        #     pickle.dump(results, f)

    # -----------------------------
    # CASE: All NaNs in var
    # -----------------------------
    if var_species_data.isna().all():

        print('Processing ', var_species_name, ': All NaNs in var')

        stage1_results = {
            'x_species_name': duplicate_species_name,
            'y_species_name': var_species_name,
            'case': 'All NaNs in var',
            'points_considered_in_odr': np.nan,
            'odr_eq_adj': np.nan,
            'metrics_init': (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan),
            'metrics_corr': (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan),
            'did_slope_improve': np.nan, 
            'slope_distance_from_1_init': np.nan,
            'slope_distance_from_1_corr': np.nan,
            'slope_distance_improvement_val': np.nan,
            'intercept_error_init': np.nan,
            'intercept_error_corr': np.nan,
            'did_rmse_improve': np.nan, 
            'rmse_percent_improvement': np.nan,
            'rmse_normalized_init': np.nan,
            'rmse_normalized_corr': np.nan,
            'r2_at_least_half': np.nan,
            'total_slope_intercept_rmse_score_init': np.nan, 
            'total_slope_intercept_rmse_score_corr': np.nan,
            'should_correction_be_applied_from_score_eval': np.nan,
            'var_data_init_full': var_species_data,
            'var_data_init_overlap': pd.Series(np.nan, index=df_all_measured_species.index),
            'dup_data_init_full': duplicate_species_data,
            'dup_data_init_overlap': pd.Series(np.nan, index=df_all_measured_species.index),
            'dup_data_corr_full': pd.Series(np.nan, index=df_all_measured_species.index),
            'dup_data_corr_overlap': pd.Series(np.nan, index=df_all_measured_species.index)
        }
    elif isinstance(duplicate_species_name, list):
        print('Processing ', var_species_name, ': Has duplicate list')

        #output from ODR fitting: return {'slope': odr_intercept, 'intercept': odr_slope,
        # 'rmse': rmse, 'norm_rmse': norm_rmse, 'r2': r2,'overlap_points_counted': overlap_points_count, 'score':score_estimate}
        run_odr_outputs = [run_odr_revised(x_data_odr = data, y_data_odr= var_species_data, mask_type = None) for data in duplicate_species_data]
        for i, out in enumerate(run_odr_outputs):
            print(i, len(out))

        # #Metrics/Evaluation of ODR fitting       
        # #next_outputs = [next_function(out, other_input) for xvar_overlap, yvar_overlap in (xvar_overlap, yvar_overlap, xvar_init_data, yvar_init_data, correction_eq, xvals_corrected, xcorrected_overlap, points_considered_in_odr)]
        # metrics_init = [fitting_metrics(xvar_overlap, yvar_overlap) for xvar_overlap, yvar_overlap, xvar_init_data, yvar_init_data, correction_eq, xvals_corrected, xcorrected_overlap, points_considered_in_odr in run_odr_outputs]
        # metrics_corr  = [fitting_metrics(xcorrected_overlap, yvar_overlap) for xvar_overlap, yvar_overlap, xvar_init_data, yvar_init_data, correction_eq, xvals_corrected, xcorrected_overlap, points_considered_in_odr in run_odr_outputs]
        
        # # #Quality control decision based off of if improvements are made by applying correction to UDAQ data
        # improvement_stats = [odr_improves(metrics_init_input = metrics_init_tuple, metrics_corr_input = metrics_corr_tuple, y_overlap = var_overlap) for metrics_init_tuple, metrics_corr_tuple, (_, var_overlap, *_) in zip(metrics_init, metrics_corr, run_odr_outputs)]
       
        # stage1_results = {
        #     'x_species_name': duplicate_species_name,
        #     'y_species_name': var_species_name,
        #     'case': 'Has duplicate list, has duplicate data not all NaNs, has vardata not all NaNs',
        #     'points_considered_in_odr': [out[7] for out in run_odr_outputs],
        #     'odr_eq_adj': [out[4] for out in run_odr_outputs], #equation used to correct the UDAQ data
        #     'metrics_init': metrics_init,
        #     'metrics_corr': metrics_corr,
        #     'did_slope_improve': [out[0] for out in improvement_stats], 
        #     'slope_distance_from_1_init': [out[1] for out in improvement_stats], 
        #     'slope_distance_from_1_corr': [out[2] for out in improvement_stats], 
        #     'slope_distance_improvement_val': [out[3] for out in improvement_stats], 
        #     'intercept_error_init': [out[4] for out in improvement_stats], 
        #     'intercept_error_corr': [out[5] for out in improvement_stats], 
        #     'did_rmse_improve': [out[6] for out in improvement_stats], 
        #     'rmse_percent_improvement': [out[7] for out in improvement_stats], 
        #     'rmse_normalized_init': [out[8] for out in improvement_stats], 
        #     'rmse_normalized_corr': [out[9] for out in improvement_stats], 
        #     'r2_at_least_half': [out[10] for out in improvement_stats], 
        #     'total_slope_intercept_rmse_score_init': [out[11] for out in improvement_stats], 
        #     'total_slope_intercept_rmse_score_corr': [out[12] for out in improvement_stats], 
        #     'should_correction_be_applied_from_score_eval': [out[13] for out in improvement_stats], 
        #     'var_data_init_full': var_species_data,
        #     'var_data_init_overlap': [out[1] for out in run_odr_outputs],
        #     'dup_data_init_full': [out[2] for out in run_odr_outputs],
        #     'dup_data_init_overlap': [out[0] for out in run_odr_outputs],
        #     'dup_data_corr_full': [out[5] for out in run_odr_outputs],
        #     'dup_data_corr_overlap': [out[6] for out in run_odr_outputs]
        # }
        
    # -----------------------------
    # CASE: All NaNs in Duplicate
    # -----------------------------
    elif duplicate_species_data.isna().all():

        print('Processing ', var_species_name, ': All NaNs in duplicate or no duplicate')

        stage1_results = {
            'x_species_name': duplicate_species_name,
            'y_species_name': var_species_name,
            'case': 'All NaNs in duplicate or no duplicate',
            'points_considered_in_odr': np.nan,
            'odr_eq_adj': np.nan,
            'metrics_init': (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan),
            'metrics_corr': (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan),
            'did_slope_improve': np.nan, 
            'slope_distance_from_1_init': np.nan,
            'slope_distance_from_1_corr': np.nan,
            'slope_distance_improvement_val': np.nan,
            'intercept_error_init': np.nan,
            'intercept_error_corr': np.nan,
            'did_rmse_improve': np.nan, 
            'rmse_percent_improvement': np.nan,
            'rmse_normalized_init': np.nan,
            'rmse_normalized_corr': np.nan,
            'r2_at_least_half': np.nan,
            'total_slope_intercept_rmse_score_init': np.nan, 
            'total_slope_intercept_rmse_score_corr': np.nan,
            'should_correction_be_applied_from_score_eval': np.nan,
            'var_data_init_full': var_species_data,
            'var_data_init_overlap': pd.Series(np.nan, index=df_all_measured_species.index),
            'dup_data_init_full': duplicate_species_data,
            'dup_data_init_overlap': pd.Series(np.nan, index=df_all_measured_species.index),
            'dup_data_corr_full': pd.Series(np.nan, index=df_all_measured_species.index),
            'dup_data_corr_overlap': pd.Series(np.nan, index=df_all_measured_species.index)
        }

    else:
        print('Processing ', var_species_name, ': Has duplicate, has duplicate data not all NaNs, has vardata not all NaNs')
        
        #ODR fitting on UDAQ data
        #output from ODR fitting: return {'slope': odr_intercept, 'intercept': odr_slope,
        # 'rmse': rmse, 'norm_rmse': norm_rmse, 'r2': r2,'overlap_points_counted': overlap_points_count, 'score':score_estimate}
        run_odr_outputs = run_odr_revised(x_data_odr = duplicate_species_data, 
                                        y_data_odr = var_species_data, 
                                        mask_type = None),
        

        #Metrics/Evaluation of ODR fitting
        (vars_overlap, vars_original) = mask_overlap(duplicate_species_data, var_species_data)
        
        metrics_init = fitting_metrics(vars_overlap[0], vars_overlap[1])

        metrics_corr  = fitting_metrics(xcorrected_overlap, vars_overlap[1])

        #Quality control decision based off of if improvements are made by applying correction to UDAQ data
        (did_slope_improve, slope_distance_i, 
        slope_distance_c, slope_distance_improvement_val, 
        intercept_err_i, intercept_err_c,
        did_rmse_improve, rmse_percent_improvement,
        rmse_norm_i, rmse_norm_c,
        r2_at_least_half, score_i, score_c,
        apply_correction_based_off_score) = odr_improves(metrics_init, metrics_corr, vars_overlap[1])

        stage1_results = {
            'x_species_name': duplicate_species_name,
            'y_species_name': var_species_name,
            'case': 'Has duplicate, has duplicate data not all NaNs, has vardata not all NaNs',
            'points_considered_in_odr': points_considered_in_odr,
            'odr_eq_adj': correction_eq, #equation used to correct the UDAQ data
            'metrics_init': metrics_init,
            'metrics_corr': metrics_corr,
            'did_slope_improve': did_slope_improve, 
            'slope_distance_from_1_init': slope_distance_i,
            'slope_distance_from_1_corr': slope_distance_c,
            'slope_distance_improvement_val': slope_distance_improvement_val,
            'intercept_error_init': intercept_err_i,
            'intercept_error_corr': intercept_err_c,
            'did_rmse_improve': did_rmse_improve, 
            'rmse_percent_improvement': rmse_percent_improvement,
            'rmse_normalized_init': rmse_norm_i, 
            'rmse_normalized_corr': rmse_norm_c,
            'r2_at_least_half': r2_at_least_half, 
            'total_slope_intercept_rmse_score_init': score_i, 
            'total_slope_intercept_rmse_score_corr': score_c,
            'should_correction_be_applied_from_score_eval': apply_correction_based_off_score,
            'var_data_init_full': var_species_data,
            'var_data_init_overlap': vars_overlap[1],
            'dup_data_init_full': duplicate_species_data,
            'dup_data_init_overlap': vars_overlap[0],
            'dup_data_corr_full': xvals_corrected,
            'dup_data_corr_overlap': xcorrected_overlap
        }
        print('Resolved ' + var_species_name + ' Stage1, one VOC')
    
    # save so ODR never recomputes
    with open(cache_file_stage1, 'wb') as f:
        pickle.dump(stage1_results, f)

    return stage1_results
def process_all_vocs_stage1_revised():
    stage1_results_all = []
    
    for row_idx in range(0, len(df_duplicates_and_tracers.index)):
        spec_name = df_duplicates_and_tracers['Varname'].iloc[row_idx]
        print('Processing Stage 1 for: ', spec_name)
        dup_name = df_duplicates_and_tracers['Duplicate_name'].iloc[row_idx]

        print('type(dup_name): ', type(dup_name))
        if pd.isna(dup_name):
            dup_species_data = pd.Series(np.nan, index=df_all_measured_species.index)
        elif pd.notna(dup_name) and ';' in dup_name:
            dup_name_list = dup_name.split('; ')
            print('Has more than one dup:', dup_name_list)
            dup_name = dup_name_list
            dup_species_data = [df_all_measured_species[dup_name[0]], df_all_measured_species[dup_name[0]]]
            print('dup_species_data type should be a list: ', type(dup_species_data))

        else:
            print('Processing Stage 1 with ', dup_name)
            dup_species_data = df_all_measured_species[dup_name]

        stage1_results = process_one_voc_stage1_revised(var_species_name = spec_name, 
                                                duplicate_species_name = dup_name, 
                                                var_species_data = df_all_measured_species[spec_name], 
                                                duplicate_species_data = dup_species_data)
        
        stage1_results_all.append(stage1_results)
    return stage1_results_all
def build_summary_stage1_revised(stage1_results_all):
def tracer_overlap(x_species_series, y_species_series):
    overlap_numerator = (~y_species_series.isna().values & ~x_species_series.isna()).sum()
    overlap_comparison = (~y_species_series.isna().values & ~x_species_series.isna()).sum() / (~y_species_series.isna()).sum()
    return overlap_numerator, overlap_comparison
###############################################################
# def stage2_tracer_scoring_draft():
#     corr_matrices = {}
#     tracer_dfs = {}
#     good_tracers_names = {}
#     tracer_filtered_data = {}
#     #Get overlap between tracers and the target species. 
#     for target, tracer in species_tracer_dict.items():
#         print('target: ', target)
#         print('tracer: ', len(tracer))
#         tracer_overlap_list = []
#         correlation_between_tracer_and_target_list = []
#         for tracer_length in range(0, len(tracer)):
#             #print(tracer[tracer_length])

#             overlap_numerator, overlap_comparison = tracer_overlap(x_species_series = df_all_measured_species[tracer[tracer_length]], 
#                                                 y_species_series = df_all_measured_species[target])

#             tracer_overlap_list.append(overlap_comparison)
#             #Get the correlation between the target species and the tracer

#             run_odr_outputs = run_odr(xdata_odr = df_all_measured_species[tracer[tracer_length]],
#                                        y_data_odr = df_all_measured_species[target])
#             #outputs from ODR are:
#             #run_odr_outputs = xvar_overlap, yvar_overlap, xvar_init_data, yvar_init_data, correction_eq, xvals_corrected, xcorrected_overlap, points_considered_in_odr
#                     #Metrics/Evaluation of ODR fitting
#             (vars_overlap, vars_original) = mask_overlap(df_all_measured_species[tracer[tracer_length]], df_all_measured_species[target])
#             metrics_init = fitting_metrics(xvar_overlap, yvar_overlap)


#             if overlap_numerator == 0:
#                 r_val = np.nan
#             else:
#                 (vars_overlap, vars_original) = mask_overlap(df_all_measured_species[tracer[tracer_length]], df_all_measured_species[target])
#                 metrics_init = fitting_metrics(vars_overlap[0], vars_overlap[1])
#                 #I only set the metrics_init to output r squared so calculate r value
#                 r_val = np.sqrt(metrics_init[3])
#             correlation_between_tracer_and_target_list.append(r_val)
#         df_tracers = pd.DataFrame({'Tracer': tracer, 'Overlap': tracer_overlap_list, 'Corr with target': correlation_between_tracer_and_target_list})
#         #print(df_tracers)

#         #Drop any tracers where overlap is less than 30%
#         df_tracers = df_tracers[df_tracers['Overlap'] >= 0.3]
#         df_tracers = df_tracers[df_tracers['Corr with target'] >= 0.2]        
#         # print(df_tracers)
#         # print(df_tracers['Tracer'])
#         print('df_tracers: ', df_tracers)
#         good_tracers = df_tracers.loc[df_tracers["Overlap"] >= 0.3, "Tracer"]
#         print('good_tracers:', good_tracers)
#         filtered_data = df_all_measured_species[df_tracers['Tracer']]
#         print('filtered_data:', filtered_data)

#         tracer_dfs[target] = df_tracers
#         good_tracers_names[target] = good_tracers
#         tracer_filtered_data[target] = filtered_data

#         corr_matrix = filtered_data.corr()
#         corr_matrices[target] = corr_matrix

#     cache_file_stage2_corr_matrix = f'{CACHE_DIR}/corr_matrices.pkl'
#     cache_file_stage2_tracer_dfs = f'{CACHE_DIR}/tracer_dfs.pkl'
#     cache_file_stage2_good_tracers_names = f'{CACHE_DIR}/good_tracers_names.pkl'
#     cache_file_stage2_tracer_filtered_data = f'{CACHE_DIR}/tracer_filtered_data.pkl'

#     with open(cache_file_stage2_corr_matrix, "wb") as f:
#         pickle.dump(corr_matrices, f)

#     with open(cache_file_stage2_tracer_dfs, "wb") as f:
#         pickle.dump(tracer_dfs, f)

#     with open(cache_file_stage2_good_tracers_names, "wb") as f:
#         pickle.dump(good_tracers_names, f)

#     with open(cache_file_stage2_tracer_filtered_data, "wb") as f:
#         pickle.dump(tracer_filtered_data, f)
    
#     #return corr_matrices
# def stage_2_correlation_clustering():
    cache_file_stage2_tracer_filtered_data = f'{CACHE_DIR}/tracer_filtered_data.pkl'
    
    with open(cache_file_stage2_tracer_filtered_data, "rb") as f:
        tracer_filtered_data = pickle.load(f)

    tracer_names = tracer_filtered_data["Ethane_WAS"].columns

    cache_file_stage2_tracer_dfs = f'{CACHE_DIR}/tracer_dfs.pkl'
    with open(cache_file_stage2_tracer_dfs, "rb") as f:
        tracer_dfs = pickle.load(f)
    print('tracer_dfs["Ethane_WAS"]:', tracer_dfs['Ethane_WAS'])
    tracer_df_for_species = tracer_dfs['Ethane_WAS']

    cache_file_stage2_corr_matrix = f'{CACHE_DIR}/corr_matrices.pkl'
    with open(cache_file_stage2_corr_matrix, "rb") as f:
        corr_matrices = pickle.load(f)

    # Use it immediately
    print(corr_matrices["Ethane_WAS"])

    # corr_matrix from earlier
    corr = corr_matrices["Ethane_WAS"].corr()

    # convert to distance matrix
    distance = 1 - np.abs(corr)

    # convert to condensed form for clustering
    dist_condensed = squareform(distance)

    # hierarchical clustering
    Z = linkage(dist_condensed, method='average')
    #set a distance threshold 
    threshold = 0.2
    clusters = fcluster(Z, threshold, criterion='distance')
    print(clusters)

    cluster_df = pd.DataFrame({"Tracer": tracer_names, "Cluster": clusters})

    print(cluster_df.sort_values("Cluster"))

    selected_tracers = []

    for cluster in cluster_df["Cluster"].unique():
        
        tracers = cluster_df.loc[
            cluster_df["Cluster"] == cluster, "Tracer"
        ]
        
        subset = tracer_df_for_species[tracer_df_for_species['Tracer'].isin(tracers)]
        
        best = subset.sort_values("Corr with target", ascending=False).iloc[0]
        
        selected_tracers.append(best["Tracer"])

    print(selected_tracers)

    #Compute VIF
    vif_threshold = 5

    
    # X_clean = X.dropna()

    # vif = pd.Series([variance_inflation_factor(X_clean.values, i) for i in range(X_clean.shape[1])],
    #                 index=X_clean.columns)

    # while max(VIF) > vif_threshold:
    # remove tracer with largest VIF
    # recompute VIF

    #Remove remaining collinear predictors
    #Run regression

            'x_species_name': duplicate_species_name,
            'y_species_name': var_species_name,
            'case': 'Has duplicate, has duplicate data not all NaNs, has vardata not all NaNs',
            'points_considered_in_odr': run_odr_outputs['points_considered_in_odr'],
            'odr_eq_adj': run_odr_outputs['correction_eq'], #equation used to correct the UDAQ data
            'odr_slope': run_odr_outputs['odr_slope'],
            'odr_intercept': run_odr_outputs['odr_intercept'],
            'metrics_init_rmse': run_odr_outputs['rmse_init'],
            'metrics_init_norm_rmse': run_odr_outputs['norm_rmse_init'],
            'metrics_init_r2': run_odr_outputs['r2_init'],
            'metrics_corr_rmse': run_odr_outputs['rmse_corr'],
            'metrics_corr_norm_rmse': run_odr_outputs['norm_rmse_corr'],
            'metrics_corr_r2': run_odr_outputs['r2_corr'],
            'did_rmse_improve': run_odr_outputs['did_rmse_improve'],
            'rmse_percent_improvement': run_odr_outputs['rmse_percent_improvement'],
            'did_r2_improve': run_odr_outputs['did_r2_improve'],
            'use_correction_eval': run_odr_outputs['use_correction_eval'],
            'var_data_init_full': var_species_data,
            'var_data_init_overlap': run_odr_outputs['yvar_overlap'],
            'dup_data_init_full': duplicate_species_data,
            'dup_data_init_overlap': run_odr_outputs['xvar_overlap'],
            'dup_data_corr_full': run_odr_outputs['xvals_corrected'],
            'dup_data_corr_overlap': run_odr_outputs['xvals_corrected_overlap']


varname': r['y_species_name'],
            'Duplicate_name': r['x_species_name'],
            'points_considered_in_odr': r['points_considered_in_odr'],
            'case': r['case'],
            'odr_eq_adj': r['odr_eq_adj'],
            'odr_slope': r['odr_slope'],
            'odr_intercept': r['odr_intercept'],
            'metrics_init_rmse': r['metrics_init_rmse'],
            'metrics_init_norm_rmse': r['metrics_init_norm_rmse'],
            'metrics_init_r2': r['metrics_init_r2'],
            'metrics_corr_rmse': r['metrics_corr_rmse'],
            'metrics_corr_norm_rmse':r['metrics_corr_norm_rmse'],
            'metrics_corr_r2': r['metrics_corr_r2'],
            'did_rmse_improve': r['did_rmse_improve'], 
            'rmse_percent_improvement': r['rmse_percent_improvement'],
            'use_correction_eval': r['use_correction_eval']

'x_species_name': duplicate_species_name,
            'y_species_name': var_species_name,
            'case': 'Has duplicate list, has duplicate data not all NaNs, has vardata not all NaNs',
            'points_considered_in_odr': [odr_results['points_considered_in_odr'] for odr_results in run_odr_outputs],
            'odr_eq_adj': [odr_results['correction_eq'] for odr_results in run_odr_outputs], #equation used to correct the UDAQ data
            'odr_slope': [odr_results['odr_slope'] for odr_results in run_odr_outputs],
            'odr_intercept': [odr_results['odr_intercept'] for odr_results in run_odr_outputs],
            'metrics_init_rmse': [odr_results['rmse_init'] for odr_results in run_odr_outputs],
            'metrics_init_norm_rmse': [odr_results['norm_rmse_init'] for odr_results in run_odr_outputs],
            'metrics_init_r2': [odr_results['r2_init'] for odr_results in run_odr_outputs],
            'metrics_corr_rmse': [odr_results['rmse_corr'] for odr_results in run_odr_outputs],
            'metrics_corr_norm_rmse': [odr_results['norm_rmse_corr'] for odr_results in run_odr_outputs],
            'metrics_corr_r2': [odr_results['r2_corr'] for odr_results in run_odr_outputs],
            'did_rmse_improve': [odr_results['did_rmse_improve'] for odr_results in run_odr_outputs],
            'rmse_percent_improvement': [odr_results['rmse_percent_improvement'] for odr_results in run_odr_outputs],
            'did_r2_improve': [odr_results['did_r2_improve'] for odr_results in run_odr_outputs],
            'use_correction_eval': [odr_results['use_correction_eval'] for odr_results in run_odr_outputs],
            'var_data_init_full': var_species_data,
            'var_data_init_overlap': [odr_results['yvar_overlap'] for odr_results in run_odr_outputs],
            'dup_data_init_full': duplicate_species_data,
            'dup_data_init_overlap': [odr_results['xvar_overlap'] for odr_results in run_odr_outputs],
            'dup_data_corr_full': [odr_results['xvals_corrected'] for odr_results in run_odr_outputs],
            'dup_data_corr_overlap': [odr_results['xvals_corrected_overlap'] for odr_results in run_odr_outputs]

def scatterplots_for_comparing_init_and_corr_odr_fits(odr_fit_results, stage_type, flag_type):
    if flag_type == 'Duplicate list':
        for dup in range(0, 1):
            x_voc_name = odr_fit_results['x_species_name'][dup]
            y_voc_name = odr_fit_results['y_species_name'][dup]

            x_init_overlap = odr_fit_results['dup_data_init_overlap'][dup]
            x_corr_overlap = odr_fit_results['dup_data_corr_overlap'][dup]
            y_overlap = odr_fit_results['var_data_init_overlap'][dup]

            rmse_init = odr_fit_results['metrics_init_rmse'][dup]
            rmse_corr = odr_fit_results['metrics_corr_rmse'][dup]

            r2_init = odr_fit_results['metrics_init_r2'][dup]
            r2_corr = odr_fit_results['metrics_corr_r2'][dup]

            #scatterplot of ML vs UDAQ initial
            fig, ax = plt.subplots(1, 2, figsize=(10,10), tight_layout=True)
            ax[0].scatter(x_init_overlap, y_overlap, s=10, alpha=0.5)

            #To draw regression line, we need a continuous line. Since some species are in ppt, we need
            #to select an appropriate scale to the step
            if np.nanmax(x_init_overlap) < 0.1:
                step = 0.001
            else:
                step = 0.1

            xrange_init = np.arange(0,np.nanmax(x_init_overlap), step)

            ax[0].plot(xrange_init, (slope_init * xrange_init + intercept_init))
            ax[0].set_title('Initial')
            if dup == 0:
                label0 = odr_fit_results['x_species_name'][0]
            elif dup == 1:
                label0 = odr_fit_results['x_species_name'][1]
        
            ax[0].set_xlabel(label0 + ' (ppb)')
            ax[0].set_ylabel(str(y_voc_name) + ' (ppb)')

            ax[0].text(0.05, 0.92, "R$^2$= " + str(round(r2_init, 3)), transform=ax[0].transAxes) 
            ax[0].text(0.05, 0.90, "RMSE:  " + str(round(rmse_init, 3)), transform=ax[0].transAxes)

            ax[1].scatter(x_corr_overlap, y_overlap, s=10, alpha=0.5)

            if np.nanmax(x_corr_overlap) < 0.1:
                step = 0.001
            else:
                step = 0.1

            xrange_corr = np.arange(0,np.nanmax(x_corr_overlap),step)
            
            ax[1].plot(x_corr_overlap, (slope_corr * xrange_corr + intercept_corr))
            ax[1].set_title('Corrected')
            
            ax[1].set_xlabel(label0 + ' (ppb)')
            ax[1].set_ylabel(y_voc_name + ' (ppb)')

            ax[1].text(0.05, 0.96, "Slope = " + str(round(slope_corr, 3)), transform=ax[1].transAxes)
            ax[1].text(0.05, 0.94, "Intercept = " + str(round(intercept_corr, 3)), transform=ax[1].transAxes)
            ax[1].text(0.05, 0.92, "R$^2$= " + str(round(r2_corr, 3)), transform=ax[1].transAxes) 
            ax[1].text(0.05, 0.90, "RMSE:  " + str(round(rmse_corr, 3)), transform=ax[1].transAxes)
            
            plt.savefig(stage_data_dirs[stage_type] + 'plots/init_corr_scatterplot_comparison_' + y_voc_name + '_' + label0 + '.png', dpi = 150)
            plt.show()
            plt.close()
    else:    
        x_voc_name = odr_fit_results['x_species_name']
        y_voc_name = odr_fit_results['y_species_name']

        x_init_overlap = odr_fit_results['dup_data_init_overlap']
        x_corr_overlap = odr_fit_results['dup_data_corr_overlap']
        y_overlap = odr_fit_results['var_data_init_overlap']

        rmse_init = odr_fit_results['metrics_init_rmse']
        rmse_corr = odr_fit_results['metrics_corr_rmse']

        r2_init = odr_fit_results['metrics_init_r2']
        r2_corr = odr_fit_results['metrics_corr_r2']

        #scatterplot of ML vs UDAQ initial
        print('Plotting array in shape: ', x_init_overlap.shape, y_overlap.shape)
        fig, ax = plt.subplots(1, 2, figsize=(10,10), tight_layout=True)
        ax[0].scatter(x_init_overlap, y_overlap, s=10, alpha=0.5)

        #To draw regression line, we need a continuous line. Since some species are in ppt, we need
        #to select an appropriate scale to the step
        if np.nanmax(x_corr_overlap) < 0.1:
            step = 0.001
        elif x_voc_name == 'H2O_CRDS' or y_voc_name == 'H2O_CRDS':
            step = 1000
        else:
            step = 0.1

        xrange_init = np.arange(0,np.nanmax(x_init_overlap), step)

        ax[0].plot(xrange_init, (slope_init * xrange_init + intercept_init))
        ax[0].set_title('Initial')
        ax[0].set_xlabel(str(x_voc_name) + ' (ppb)')
        ax[0].set_ylabel(str(y_voc_name) + ' (ppb)')

        ax[0].text(0.05, 0.92, "R$^2$= " + str(round(r2_init, 3)), transform=ax[0].transAxes) 
        ax[0].text(0.05, 0.90, "RMSE:  " + str(round(rmse_init, 3)), transform=ax[0].transAxes)


        ax[1].scatter(x_corr_overlap, y_overlap, s=10, alpha=0.5)

        if np.nanmax(x_corr_overlap) < 0.1:
            step = 0.001
        elif x_voc_name == 'H2O_CRDS' or y_voc_name == 'H2O_CRDS':
            step = 1000
        else:
            step = 0.1

        xrange_corr = np.arange(0,np.nanmax(x_corr_overlap),step)
        
        ax[1].plot(xrange_corr, (slope_corr * xrange_corr + intercept_corr))
        ax[1].set_title('Corrected')
        ax[1].set_xlabel(str(x_voc_name) + ' (ppb)')
        ax[1].set_ylabel(str(y_voc_name) + ' (ppb)')

        ax[1].text(0.05, 0.92, "R$^2$= " + str(round(r2_corr, 3)), transform=ax[1].transAxes) 
        ax[1].text(0.05, 0.90, "RMSE:  " + str(round(rmse_corr, 3)), transform=ax[1].transAxes)

        plt.savefig(stage_data_dirs[stage_type] + 'plots/init_corr_scatterplot_comparison_' + str(y_voc_name) + '_' +str(x_voc_name) + '.png', dpi = 150)
        plt.show()
        plt.close()