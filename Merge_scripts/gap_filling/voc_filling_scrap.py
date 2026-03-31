def scatterplot_co_and_voc_result_input(df_co_and_voc_results, idx):
    #Enter results from plotting scatterplot relationship between VOC and CO
    #Looking for: Clear positive trend, not a vertical cloud, no obvious curvature or thresholds
    #Vertical cloud means: VOC varies independently of CO, CO is not controlling or constraining VOC variability, CO is a weak predictor
    outputs = ['Trend', 'Vertical_cloud', 'Curvature', 'Thresholds']
    if df_co_and_voc_results.loc[idx, outputs].notna().all():
        return df_co_and_voc_results  # already done

    trend_input = input("Is the scatterplot trend positive? Press 't' for Positive, 'f' for Negative, 'nt' for No Trend, 'u' for Unsure.")
    if trend_input.lower()=='t':
        trend_eval = 'Positive'
    elif trend_input.lower() =='f':
        trend_eval = 'Negative'
    elif trend_input.lower()=='nt':
        trend_eval = 'No Trend'
    elif trend_input.lower()=='u':
        trend_eval = 'Unsure'
    else:
        print('try again')

    vertical_cloud_input = input("Is there a vertical cloud (column at one CO value)? Press 't' for True, 'f' for False.")
    if vertical_cloud_input.lower()=='t':
        vertical_cloud_eval = True
    elif vertical_cloud_input.lower() =='f':
        vertical_cloud_eval = False
    else:
        print('try again')

    curvature_input = input("Is there any curvature? Press 't' for yes, 'f' for no.")
    if curvature_input.lower()=='t':
        curvature_eval = True
    elif curvature_input.lower() =='f':
        curvature_eval = False
    else:
        print('try again')

    thresholds_input = input("Are there any thresholds (Flat or noisy VOC at low CO, then sudden linear trend after CO exceeds some value)? Press 't' for yes, 'f' for no.")
    if thresholds_input.lower()=='t':
        thresholds_eval = True
    elif thresholds_input.lower() =='f':
        thresholds_eval = False
    else:
        print('try again')

    df_co_and_voc_results.loc[idx, 'Trend'] = trend_eval
    df_co_and_voc_results.loc[idx, 'Vertical_cloud'] = vertical_cloud_eval
    df_co_and_voc_results.loc[idx, 'Curvature'] = curvature_eval
    df_co_and_voc_results.loc[idx, 'Thresholds'] = thresholds_eval

    return df_co_and_voc_results
def scatterplot_log_co_and_log_voc(voc_nonans, co_nonans, df_co_and_voc_results, idx, vocname_udaq):
    #log-log scatter plot
    #Good sign: Linear-ish cloud, Constant spread across range
    #Bad sign: Two regimes, Strong curvature, Fan-shaped scatter
    mask_pos = (voc_nonans > 0) & (co_nonans > 0)
    voc_p = voc_nonans[mask_pos]
    co_p  = co_nonans[mask_pos]

    plt.figure(figsize=(5,5))
    plt.scatter(np.log10(co_p), np.log10(voc_p),
                s=10, alpha=0.3)
    plt.xlabel('log CO')
    plt.ylabel('log ' + vocname_udaq + ' (ppb)')
    plt.tight_layout()
    #plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/co_tracer_plots/scatterplot_logco_vs_logvoc/log_co_with_log_'+ str(vocname_udaq) + '_scatterplot.png', dpi =300)
    plt.show()

    #Quantify correlation with Pearson and Spearman
    #where >= 0.6 Spearman is likely promising
    r_p, _ = pearsonr(co_p, voc_p)
    r_s, _ = spearmanr(co_p, voc_p)

    print(f'Pearson r {vocname_udaq} = {r_p:.2f}')
    print(f'Spearman r {vocname_udaq}= {r_s:.2f}')

    outputs = ['Pearson_correlation_coefficient_logscatter', 'Spearman_correlation_coefficient_logscatter']
    if df_co_and_voc_results.loc[idx, outputs].notna().all():
        return df_co_and_voc_results  # already done

    df_co_and_voc_results.loc[idx, 'Pearson_correlation_coefficient_logscatter'] = r_p
    df_co_and_voc_results.loc[idx, 'Spearman_correlation_coefficient_logscatter'] = r_s

    return voc_p, co_p, df_co_and_voc_results
def scatterplot_log_co_and_log_voc_result_input(df_co_and_voc_results, idx):
    #Enter results from plotting scatterplot relationship between log VOC and log CO
    #Good sign: Linear-ish cloud, Constant spread across range
    #Bad sign: Two regimes, Strong curvature, Fan-shaped scatter
    outputs = ['Linear_cloud_logscatter', 'Constant_spread_logscatter', 'Multiple_regimes_logscatter', 'Strong_curvature_logscatter', 'Fan_scatter_logscatter']
    if df_co_and_voc_results.loc[idx, outputs].notna().all():
        return df_co_and_voc_results  # already done
    
    linear_cloud_input = input("Is the scatter linear-ish? Press 't' for Linear, 'f' for nonlinear, 'nt' for No Trend, 'u' for Unsure.")
    if linear_cloud_input.lower()=='t':
        linear_cloud_eval = 'Linear'
    elif linear_cloud_input.lower() =='f':
        linear_cloud_eval = 'Nonlinear'
    elif linear_cloud_input.lower()=='nt':
        linear_cloud_eval = 'No Trend'
    elif linear_cloud_input.lower()=='u':
        linear_cloud_eval = 'Unsure'
    else:
        print('try again')

    constant_spread_input = input("Is there a constant spread across the range? Press 't' for Yes, 'f' for No.")
    if constant_spread_input.lower()=='t':
        constant_spread_eval = True
    elif constant_spread_input.lower() =='f':
        constant_spread_eval = False
    else:
        print('try again')

    multiple_regimes_input = input("Are there multiple regimes? Press 't' for One, 'f' for Multiple.")
    if multiple_regimes_input.lower()=='t':
        multiple_regimes_eval = 'One'
    elif multiple_regimes_input.lower() =='f':
        multiple_regimes_eval = 'Multiple'
    else:
        print('try again')

    strong_curvature_input = input("Is there strong curvature? Press 't' for yes, 'f' for no.")
    if strong_curvature_input.lower()=='t':
        strong_curvature_eval = True
    elif strong_curvature_input.lower() =='f':
        strong_curvature_eval = False
    else:
        print('try again')

    fanscatter_input = input("Is there a fan-shaped scatter? Press 't' for yes, 'f' for no.")
    if fanscatter_input.lower()=='t':
        fanscatter_eval = True
    elif fanscatter_input.lower() =='f':
        fanscatter_eval = False
    else:
        print('try again')

    df_co_and_voc_results.loc[idx, 'Linear_cloud_logscatter'] = linear_cloud_eval
    df_co_and_voc_results.loc[idx, 'Constant_spread_logscatter'] = constant_spread_eval
    df_co_and_voc_results.loc[idx, 'Multiple_regimes_logscatter'] = multiple_regimes_eval
    df_co_and_voc_results.loc[idx, 'Strong_curvature_logscatter'] = strong_curvature_eval
    df_co_and_voc_results.loc[idx, 'Fan_scatter_logscatter'] = fanscatter_eval

    return df_co_and_voc_results
def fit_log_model(voc_p, co_p, df_co_and_voc_results, idx):
    #Fit the log-log model, inspect residuals
    #Look for:
    # Slope near constant across subsets
    # Reasonable R2 (0.4–0.8 is common, not bad)
    # No obvious leverage points dominating

    x_logco = sm.add_constant(np.log10(co_p))
    y_logvoc = np.log10(voc_p)

    model_logfit = sm.OLS(y_logvoc, x_logco).fit()
    print(model_logfit.summary())

    print(f"Slope = {model_logfit.params[1]:.3f}," f"Intercept = {model_logfit.params[0]:.3f}," f"R$^2$ = {model_logfit.rsquared:.3f}")
    model_logfit_slope = model_logfit.params[1]
    model_logfit_intercept = model_logfit.params[0]
    model_logfit_rsquared = model_logfit.rsquared

    outputs = ['Slope_logfit', 'Intercept_logfit', 'Rsquared_logfit']
    if df_co_and_voc_results.loc[idx, outputs].notna().all():
        return df_co_and_voc_results  # already done

    df_co_and_voc_results.loc[idx, 'Slope_logfit'] = model_logfit_slope
    df_co_and_voc_results.loc[idx, 'Intercept_logfit'] = model_logfit_intercept
    df_co_and_voc_results.loc[idx, 'Rsquared_logfit'] = model_logfit_rsquared

    return model_logfit, df_co_and_voc_results
def residuals_vs_co(model_logfit, vocname_udaq):
    #Residuals vs CO (checks extrapolation risk)
    #Bad sign:
    # Residuals systematically increase at low or high CO
    # Suggests regime change or missing covariate

    resid_logfit = model_logfit.resid

    plt.figure(figsize=(5,4))
    plt.scatter(np.log10(co_p), resid_logfit, s=10, alpha=0.3)
    plt.axhline(0, color='k')
    plt.xlabel('log CO')
    plt.ylabel('Residual (log '+ str(vocname_udaq) +')')
    plt.tight_layout()
    #plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/co_tracer_plots/scatterplot_logco_vs_residual/log_co_with_residual_log_'+ str(vocname_udaq) + '_scatterplot.png', dpi =300)
    plt.show()

    return resid_logfit
def residuals_vs_co_result_input(df_co_and_voc_results, idx):
    outputs = ['Residual_systemic_increase_at_low_or_high_co']
    if df_co_and_voc_results.loc[idx, outputs].notna().all():
        return df_co_and_voc_results  # already done
    
    residuals_systemic_increase_input = input("Is there a systemic increase in residuals at low or high CO? Press 'l' for fan shape (widening) at low CO, 'h' for fan shape (widening) at high CO, 'n' for Neither.")
    if residuals_systemic_increase_input.lower()=='l':
        residuals_systemic_increase_eval = 'Fan at low CO'
    elif residuals_systemic_increase_input.lower() =='h':
        residuals_systemic_increase_eval = 'Fan at high CO'
    elif residuals_systemic_increase_input.lower() =='n':
        residuals_systemic_increase_eval = 'No'
    else:
        print('try again')

    df_co_and_voc_results.loc[idx, 'Residuals_systemic_increase_at_low_or_high_co'] = residuals_systemic_increase_eval

    return df_co_and_voc_results

def residuals_vs_hour_of_day(voc_p, resid_logfit, vocname_udaq):
    #Residuals vs hour of day (checks diurnal confounding)
    #If structure remains:
    # CO alone isn’t enough
    # Consider adding NOx or stratifying by time/season

    #Calculate residuals
    residuals = resid_logfit
    hours = voc_p.index.hour  # get hour of day

    #Make a DataFrame for grouping
    resid_df = pd.DataFrame({'hour': hours, 'residual': residuals})

    #Calculate mean and std per hour
    hourly_mean = resid_df.groupby('hour')['residual'].mean()
    hourly_std  = resid_df.groupby('hour')['residual'].std()

    #Flag hours with mean residual exceeding threshold (e.g., ±0.1 log units)
    threshold = 0.1
    flagged_hours = hourly_mean[hourly_mean.abs() > threshold].index.tolist()
    print("Hours with systematic residual bias:", flagged_hours)

    plt.figure(figsize=(10,5))
    plt.errorbar(hourly_mean.index, hourly_mean, yerr=hourly_std, fmt='o', capsize=4, label='Mean ± Std')
    plt.axhline(0, color='k', linestyle='--')
    plt.xlabel('Hour of Day')
    plt.ylabel('Residuals (log' + str(vocname_udaq) + ')')
    plt.title('Residuals vs Hour of Day')
    plt.xticks(range(0,24))
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/co_tracer_plots/scatterplot_residuals_vs_hour_of_day/residuals_log_'+ str(vocname_udaq) + 'vs_hour_of_day_scatterplot.png', dpi =300)
    plt.show()

    # hours = voc_p.index.hour

    # plt.figure(figsize=(6,4))
    # plt.scatter(hours, resid_logfit, s=10, alpha=0.3)
    # plt.axhline(0, color='k')
    # plt.xlabel('Hour of day')
    # plt.ylabel('Residual (log '+ str(vocname_udaq) +')')
    # #plt.title('Residuals vs hour of day')
    # plt.tight_layout()
    # #plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/co_tracer_plots/scatterplot_residuals_vs_hour_of_day/residuals_log_'+ str(vocname_udaq) + 'vs_hour_of_day_scatterplot.png', dpi =300)
    # plt.show()

def residuals_vs_hour_of_day_result_input(df_co_and_voc_results, idx):
    outputs = ['Structure_remaining']
    if df_co_and_voc_results.loc[idx, outputs].notna().all():
        return df_co_and_voc_results  # already done
    
    structure_remaining_input = input("Does structure remain? Press 't' for Yes, 'f' for No.")
    if structure_remaining_input.lower()=='t':
        structure_remaining_eval = 'Linear'
    elif structure_remaining_input.lower() =='f':
        structure_remaining_eval = 'Nonlinear'
    else:
        print('try again')
    
    df_co_and_voc_results.loc[idx, 'Structure_remaining'] = structure_remaining_eval

    return df_co_and_voc_results

def stability_check(voc_p, co_p, df_co_and_voc_results, idx):
    #Stability check 
    mid = voc_p.index[len(voc_p)//2]

    early = voc_p.index < mid
    late  = voc_p.index >= mid

    def fit_loglog(v, c):
        X = sm.add_constant(np.log10(c))
        y = np.log10(v)
        res = sm.OLS(y, X).fit()
        slope = res.params[1]      # slope
        intercept = res.params[0]  # intercept
        slope_se = res.bse[1]      # standard error of slope
        r2 = res.rsquared
        return slope, intercept, slope_se, r2
    
    outputs = ['Slope_early', 'Intercept_early', 'Slope_standard_error_early', 'R2_early', 'Slope_late', 'Intercept_late',  'Slope_standard_error_late',  'R2_late']
    
    if df_co_and_voc_results.loc[idx, outputs].notna().all():
        return df_co_and_voc_results  # already done

    slope_early, int_early, slope_se_early, r2_early = fit_loglog(voc_p[early], co_p[early])
    slope_late,  int_late,  slope_se_late,  r2_late  = fit_loglog(voc_p[late],  co_p[late])

   
    df_co_and_voc_results.loc[idx, 'Slope_early'] = slope_early
    df_co_and_voc_results.loc[idx, 'Intercept_early'] = int_early
    df_co_and_voc_results.loc[idx, 'Slope_standard_error_early'] = slope_se_early
    df_co_and_voc_results.loc[idx, 'R2_early'] = r2_early
    df_co_and_voc_results.loc[idx, 'Slope_late'] = slope_late
    df_co_and_voc_results.loc[idx, 'Intercept_late'] = int_late
    df_co_and_voc_results.loc[idx, 'Slope_standard_error_late'] = slope_se_late
    df_co_and_voc_results.loc[idx, 'R2_late'] = r2_late 

    return df_co_and_voc_results

for col in df_ml_udaq_initial_and_filled_ml.columns:
    if col.startswith('Filled_') and col.endswith('_ML_with_UDAQ_Adjusted'):
        filled_vocname = col[len('Filled_'):-len('_ML_with_UDAQ_Adjusted')]
        print('filled_voc_name: ', filled_vocname)
        filled_colname = f'Filled_{filled_vocname}_ML_with_UDAQ_Adjusted'

        if filled_vocname not in reverse_mapping:
            continue
        
        ml_name_init = reverse_mapping[filled_vocname]
        init_colname = f'NOAA_{ml_name_init}_Initial'
        print('name_init: ', ml_name_init)

        if init_colname not in df_ml_udaq_initial_and_filled_ml.columns:
            continue
        
        if df_ml_udaq_initial_and_filled_ml[filled_colname].isna().all():
            print(filled_vocname, ' has no measurements.')

        else:
#             plot_carbon_monoxide_tracer_timeseries(
#                 df_index = df_ml_udaq_initial_and_filled_ml.index, 
#                 filled_col = df_ml_udaq_initial_and_filled_ml[filled_colname], 
#                 vocname_udaq = filled_vocname
#             )

            (voc_nonans, co_nonans) = scatterplot_co_and_voc(
                filled_col = df_ml_udaq_initial_and_filled_ml[filled_colname],
                vocname_udaq = filled_vocname
            )
            
            co_tracer_voc_regression(
                filled_col = df_ml_udaq_initial_and_filled_ml[filled_colname],
                vocname_udaq = filled_vocname)
            
#             (voc_p, co_p, pearson_r_val, spearman_r_val) = scatterplot_log_co_and_log_voc(
#                 voc_nonans, 
#                 co_nonans, 
#                 vocname_udaq = filled_vocname
#             )

#             (model_logfit, model_logfit_slope, model_logfit_rsquared) = fit_log_model(voc_p, co_p)
    
#             resid_logfit = residuals_vs_co(
#                 model_logfit,
#                 vocname_udaq = filled_vocname
#             )

#             residuals_vs_hour_of_day(
#                 voc_p, 
#                 resid_logfit, 
#                 vocname_udaq = filled_vocname
#             )
        
        #Decision rule (simple and honest)
        # I’d be comfortable using CO-based gap filling if:
        # Spearman ≥ 0.6
        # Log–log relationship looks linear
        # Residuals are flat vs hour
        # Slopes stable across time
        # CV error < natural variability
        # If not → don’t force it.

#FOR APPLYING VOC TRACER
# calculated_fields = ['Pearson_correlation_coefficient_logscatter', 'Spearman_correlation_coefficient_logscatter', 'Slope_logfit', 'Intercept_logfit', 'Rsquared_logfit', 'Slope_early', 'Intercept_early', 'Slope_standard_error_early', 'R2_early', 'Slope_late',  'Intercept_late',  'Slope_standard_error_late',  'R2_late']
# input_required_fields = ['Trend', 'Vertical_cloud', 'Curvature', 'Thresholds', 'Residual_systemic_increase_at_low_or_high_co', 'Structure_remaining', 'Linear_cloud_logscatter', 'Constant_spread_logscatter', 'Multiple_regimes_logscatter', 'Strong_curvature_logscatter', 'Fan_scatter_logscatter']
# all_fields = ['ML_VOC_Name'] + ['UDAQ_VOC_Name'] + calculated_fields + input_required_fields
# tracer_csv_savepath = dirpath + '/Merge_scripts/calibration_adjustments/co_tracer_info.csv'

# if os.path.exists(tracer_csv_savepath):
#     df_co_and_voc_results = pd.read_csv(tracer_csv_savepath)
# else:
#     df_co_and_voc_results = pd.DataFrame(columns=all_fields)
    
# #############################

# def get_row_index(df_co_and_voc_results, vocname_ml, vocname_udaq):
#     matches = df_co_and_voc_results.index[df_co_and_voc_results['ML_VOC_Name'] == vocname_ml]
#     if len(matches) > 0:
#         return matches[0]
    
#     # create new row
#     row = {k: pd.NA for k in all_fields}
#     row['ML_VOC_Name'] = vocname_ml
#     row['UDAQ_VOC_Name'] = vocname_udaq

#     df_co_and_voc_results.loc[len(df_co_and_voc_results)] = row
#     return len(df_co_and_voc_results) - 1

# for col in df_ml_udaq_initial_and_filled_ml.columns:
#     if col.startswith('Filled_') and col.endswith('_ML_with_UDAQ_Adjusted'):
#         filled_vocname = col[len('Filled_'):-len('_ML_with_UDAQ_Adjusted')]
#         print('filled_voc_name: ', filled_vocname)
#         filled_colname = f'Filled_{filled_vocname}_ML_with_UDAQ_Adjusted'

#         if filled_vocname not in reverse_mapping:
#             continue
        
#         ml_name_init = reverse_mapping[filled_vocname]
#         init_colname = f'NOAA_{ml_name_init}_Initial'
#         print('name_init: ', ml_name_init)

#         if init_colname not in df_ml_udaq_initial_and_filled_ml.columns:
#             continue

#         idx = get_row_index(df_co_and_voc_results, vocname_ml = ml_name_init, vocname_udaq = filled_vocname)
        
#         if df_ml_udaq_initial_and_filled_ml[filled_colname].isna().all():
#             print(filled_vocname, ' has no measurements.')
        
#             # explicitly mark all auto + annotation fields as NaN
#             for col in calculated_fields + input_required_fields:
#                 df_co_and_voc_results.loc[idx, col] = np.nan

#             df_co_and_voc_results.to_csv(tracer_csv_savepath, index=False)
#             continue

#         else:
#             plot_carbon_monoxide_tracer_timeseries(
#                 df_index = df_ml_udaq_initial_and_filled_ml.index, 
#                 filled_col = df_ml_udaq_initial_and_filled_ml[filled_colname], 
#                 vocname_udaq = filled_vocname)

#             (voc_nonans, co_nonans) = scatterplot_co_and_voc(
#                 filled_col = df_ml_udaq_initial_and_filled_ml[filled_colname],
#                 vocname_udaq = filled_vocname)

#             df_co_and_voc_results = scatterplot_co_and_voc_result_input(df_co_and_voc_results, idx)

#             (voc_p, co_p, df_co_and_voc_results) = scatterplot_log_co_and_log_voc(
#                 voc_nonans, 
#                 co_nonans, 
#                 df_co_and_voc_results, 
#                 idx,
#                 vocname_udaq = filled_vocname)

#             df_co_and_voc_results = scatterplot_log_co_and_log_voc_result_input(df_co_and_voc_results, idx)

#             (model_logfit, df_co_and_voc_results) = fit_log_model(voc_p, co_p, df_co_and_voc_results, idx)

#             resid_logfit = residuals_vs_co(
#                 model_logfit,
#                 vocname_udaq = filled_vocname)
            
#             df_co_and_voc_results = residuals_vs_co_result_input(df_co_and_voc_results, idx)

#             residuals_vs_hour_of_day(
#                 voc_p, 
#                 resid_logfit, 
#                 vocname_udaq = filled_vocname)
            
#             df_co_and_voc_results = residuals_vs_hour_of_day_result_input(df_co_and_voc_results, idx)

#             df_co_and_voc_results = stability_check(voc_p, co_p, df_co_and_voc_results, idx)
            
#             #save VOC data
#             df_co_and_voc_results = df_co_and_voc_results.to_csv(tracer_csv_savepath, index=False)

# new_voc_interpolated_2hr_dict = {}

# for col in df_ml_udaq_initial_and_filled_ml.columns:
#     if col.startswith('Filled_') and col.endswith('_ML_with_UDAQ_Adjusted'):
#         filled_vocname = col[len('Filled_'):-len('_ML_with_UDAQ_Adjusted')]
#         print('filled_voc_name: ', filled_vocname)
#         filled_colname = f'Filled_{filled_vocname}_ML_with_UDAQ_Adjusted'

#         if filled_vocname not in reverse_mapping:
#             continue
        
#         ml_name_init = reverse_mapping[filled_vocname]
#         init_colname = f'NOAA_{ml_name_init}_Initial'
#         print('name_init: ', ml_name_init)

#         if init_colname not in df_ml_udaq_initial_and_filled_ml.columns:
#             continue
        
#         if df_ml_udaq_initial_and_filled_ml[filled_colname].isna().all():
#             print(filled_vocname, ' has no measurements.')
#             co_and_voc_results.append({
#                 'ML_VOC_Name': ml_name_init,
#                 'UDAQ_VOC_Name': filled_vocname,
#                 'Trend': np.nan,
#                 'Vertical_cloud': np.nan,
#                 'Curvature': np.nan,
#                 'Thresholds': np.nan,
#                 'Pearson_correlation_coefficient_logscatter': np.nan,
#                 'Spearman_correlation_coefficient_logscatter': np.nan,
#                 'Linear_cloud_logscatter': np.nan, 
#                 'Constant_spread_logscatter': np.nan, 
#                 'Multiple_regimes_logscatter': np.nan, 
#                 'Strong_curvature_logscatter': np.nan, 
#                 'Fan_scatter_logscatter': np.nan,
#                 'Slope_logfit': np.nan, 
#                 'Rsquared_logfit': np.nan
#             })

#         else:
#             plot_carbon_monoxide_tracer_timeseries(
#                 df_index = df_ml_udaq_initial_and_filled_ml.index, 
#                 filled_col = df_ml_udaq_initial_and_filled_ml[filled_colname], 
#                 vocname_udaq = filled_vocname
#             )

#             (voc_nonans, co_nonans) = scatterplot_co_and_voc(
#                 filled_col = df_ml_udaq_initial_and_filled_ml[filled_colname],
#                 vocname_udaq = filled_vocname
#                 )

#             (trend_eval, vertical_cloud_eval, curvature_eval, thresholds_eval) = scatterplot_co_and_voc_result_input()

#             co_and_voc_results.append({
#                 'ML_VOC_Name': ml_name_init,
#                 'UDAQ_VOC_Name': filled_vocname,
#                 'Trend': trend_eval,
#                 'Vertical_cloud': vertical_cloud_eval,
#                 'Curvature': curvature_eval,
#                 'Thresholds': thresholds_eval
#             })

#             (voc_p, co_p, pearson_r_val, spearman_r_val) = scatterplot_log_co_and_log_voc(
#                 voc_nonans, 
#                 co_nonans, 
#                 vocname_udaq = filled_vocname
#             )
#             co_and_voc_results.append({
#                 'Pearson_correlation_coefficient_logscatter': pearson_r_val, 
#                 'Spearman_correlation_coefficient_logscatter': spearman_r_val
#             })

#             (linear_cloud_eval, constant_spread_eval, multiple_regimes_eval, strong_curvature_eval, fanscatter_eval) = scatterplot_log_co_and_log_voc_result_input()
#             co_and_voc_results.append({
#                 'Linear_cloud_logscatter': linear_cloud_eval, 
#                 'Constant_spread_logscatter': constant_spread_eval, 
#                 'Multiple_regimes_logscatter': multiple_regimes_eval, 
#                 'Strong_curvature_logscatter': strong_curvature_eval, 
#                 'Fan_scatter_logscatter': fanscatter_eval
#             })

#             (model_logfit, model_logfit_slope, model_logfit_rsquared) = fit_log_model(voc_p, co_p)
#             co_and_voc_results.append({
#                 'Slope_logfit': model_logfit_slope, 
#                 'Rsquared_logfit': model_logfit_rsquared
#             })

#             resid_logfit = residuals_vs_co(
#                 model_logfit,
#                 vocname_udaq = filled_vocname)

#             residuals_vs_hour_of_day(
#                 voc_p, 
#                 resid_logfit, 
#                 vocname_udaq = filled_vocname)
            
#             #save VOC data
            
#             df_co_and_voc_results = pd.DataFrame(co_and_voc_results).to_csv(tracer_savepath, index=False)
            
            
#string_inputs: 
#Trend: 'Positive', 'Negative', 'No Trend', 'Unsure'
#Linear_cloud_logscatter: 'Linear', 'Nonlinear', 'No Trend', 'Unsure'
#Multiple_regimes_logscatter: 'One', 'Multiple'







#         # plot_merged_data_time_series(
#         #     df_index = df_ml_udaq_initial_and_filled_ml.index,
#         #     init_col = df_ml_udaq_initial_and_filled_ml[init_colname],
#         #     filled_col = df_ml_udaq_initial_and_filled_ml[filled_colname],
#         #     vocname = filled_voc_name
#         # )

#         # diurnal_cycle_post_merge_check(
#         #     df_index = df_ml_udaq_initial_and_filled_ml.index,
#         #     init_col = df_ml_udaq_initial_and_filled_ml[init_colname],
#         #     filled_col = df_ml_udaq_initial_and_filled_ml[filled_colname],
#         #     vocname = filled_voc_name
#         # )

#         ds_filled_interp_final = apply_interpolation_to_small_gaps_2hrs(
#             df_index = df_ml_udaq_initial_and_filled_ml.index,
#             init_col = df_ml_udaq_initial_and_filled_ml[init_colname],
#             filled_col = df_ml_udaq_initial_and_filled_ml[filled_colname],
#             vocname = filled_vocname
#         )

#         plot_interpolation_2hr_time_series(
#             df_index = df_ml_udaq_initial_and_filled_ml.index,
#             init_col = df_ml_udaq_initial_and_filled_ml[init_colname],
#             filled_col = df_ml_udaq_initial_and_filled_ml[filled_colname],
#             vocname = filled_vocname,
#             interpolated_index = ds_filled_interp_final.index,
#             interpolated_col = ds_filled_interp_final
#         )
#         if 'ML_NoVar' in ml_name_init:
#             new_voc_interpolated_2hr_dict[filled_vocname] = ds_filled_interp_final
#         elif 'Formaldehyde' in ml_name_init:
#             new_voc_interpolated_2hr_dict['HCHO_CRDS'] = ds_filled_interp_final
#         else:
#             new_voc_interpolated_2hr_dict[ml_name_init] = ds_filled_interp_final
# df_new_voc_interpolated_2hr = pd.DataFrame(new_voc_interpolated_2hr_dict)
# print(df_new_voc_interpolated_2hr.index)
# print(df_new_voc_interpolated_2hr.columns.values)
# savepath = dirpath + '/Merge_scripts/calibration_adjustments/ml_udaq_initial_and_filled_ml_with_udaq_adjusted_then_interpolated.csv'
# df_new_voc_interpolated_2hr.to_csv(savepath)

#         # diurnal_fifteen_minute_avgs(
#         #     interpolated_index = ds_filled_interp_final.index, 
#         #     interpolated_col = ds_filled_interp_final, 
#         #     vocname = filled_voc_name
#         # )

#     else:
#         pass

# def subset_day(date_time_start, date_time_stop,file_subset_name, var_name):
#     """
#     This function is for if you want to interpolate a subset of the days for the campaign.

#     INPUTS:
#         date_time_start: A string that sets the beginning of your date and time subset, in format YYYY-MM-DD HH:MM:SS
#         date_time_stop: A string that sets the end of your date and time subset, in format YYYY-MM-DD HH:MM:SS. 
#                         NOTE: THIS IS AN INCLUSIVE VALUE so date_time_stop= "2024-08-08 23:30:00" includes the 23:30:00 value.
#         file_subset_name: A string with format separated by underscore:
#             date start
#             date end
#             time for averaging (such as parked data with 30 min averages)
#             CSL_mobile_lab
#             parked / driving
#             with_interp (to indicate that it includes the interpolation)
#         var_name: A string that sets the name of the MATLAB variable when you import the file into MATLAB.
#     """
#     df_subsetdays = df_ml_data.sort_index().loc[date_time_start:date_time_stop]
#     df_interp_subset=df_subsetdays.copy()

#     #Nell turned all benzaldehyde, styrene, HONO negative values into NaNs = 0
#     neg_species = ['Benzaldehyde_PTR','Styrene_PTR', 'HONO_CIMS']

#     for i,col in enumerate(vars2fill):
#         if col in neg_species:
#             df_interp_subset[col] = df_interp_subset[col].mask(df_interp_subset[col] < 0, 0)
#         elif col == 'Lon': #avoids making longitude points into NaNs
#             df_interp_subset['Lon'] = df_interp_subset['Lon'].interpolate(method='linear')
#         else:
#             # Set any negative values to NaN so we can interp them... 
#             df_interp_subset[col] = df_interp_subset[col].mask(df_interp_subset[col] < 0, np.nan)
#         #Benzaldehyde has all NaNs so substitute with zeros instead
#         if col == 'Benzaldehyde_PTR':
#             df_interp_subset[col] =  df_interp_subset[col].mask(np.isnan(df_interp_subset[col]), 0)
#         else:
#             pass
#         # Calc number of points that are negative or Nans: 
#         n_baddies= len([item for item in df_ml_data[col] if item <0 or np.isnan(item)]) 
        
#         if n_baddies > 0: 
#             #apply the linear interpolation
#             df_interp_subset[col] = df_interp_subset[col].interpolate(method='linear')

#     #get a ratio for jNO2 measured to TUV
#     df_interp_subset['jNO2_ratio'] = df_interp_subset['jNO2_meas']/df_interp_subset['jNO2']
#     #level out the inf values and values that are too high for the jNO2 ratio
#     msk = ((df_interp_subset['jNO2_ratio'] ==np.inf)  | (df_interp_subset['jNO2_ratio'] >10) )
#     df_interp_subset.loc[msk,'jNO2_ratio'] = 1.0

#     savepath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/' + 'F0AM_filled/' + file_subset_name + '.csv'
#     df_interp_subset.to_csv(savepath)
#     print('Saved CSV to:' + savepath)

#     # Convert the dataframe to a nested dictionary (so scipy can output to a matlab structure!) 
#     ddict=dataframe_to_nested_dict(df_interp_subset)

#     # Sort alphabetically so not annoying in MATLAB...  
#     ddict= OrderedDict(sorted(ddict.items())) 

#     # Save the USOS data in an output .mat file: 
#     outpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/F0AM-4.4.2/Campaign_Data/matlab_merge/parked/original/'
#     matfilename = file_subset_name + '.mat'
#     savemat(outpath+matfilename,{var_name: ddict})
#     print('Saved MATLAB file to:' + outpath + matfilename)



#     diurnal_time_series_with_gap_counts(noaa_species = voc_spec)
# #     #gap_measuring_total(noaa_species = voc_spec)
#     for date_i in campaign_dateslist:
# #         #gap_measuring_daily(noaa_voc = voc_spec, dates = date_i)
#         plot_both(noaa, udaq, 
#               noaa_species = voc_spec, 
#               udaq_species = mapping[voc_spec], 
#               dates = date_i)




                

# # print(mapping.keys())
# for spec in mapping.keys():
# 	print(mapping[spec])
# 	plot_both(noaa, udaq, spec, mapping[spec], date_i)

# 	preferred=input(f"Which is better for {mapping[spec]} on {date_i}? Press U for 'UDAQ', 'N' for NOAA, 'NA' for No Data")
# 	if preferred.lower()=='n':
# 		pref=f"NOAA.{spec}"
# 	elif preferred.lower() =='u':
# 		pref=f"UDAQ.{mapping[spec]}"
# 	elif preferred.lower()=='na':
# 		pref='No Data'
# 	else:
# 		print('try again')

# 	which_is_better(file, date_i, mapping[spec], pref)

def diurnal_time_series_with_gap_counts(noaa_species):
    #make average diurnal plot with overall counts of gaps
    if 'ML_NoVar' in noaa_species:
        nan_series = pd.Series(index=hour_range)
        noaa_var = nan_series

    else:
        noaa['hour']=noaa.index.hour
        # total_by_hour = noaa[noaa_species].groupby(noaa.index.hour).size()
        # print(total_by_hour)
        display(noaa[noaa_species])

        # Compute daily diurnal cycles
        diurnal_daily = noaa[noaa_species].groupby([noaa.index.date, noaa.index.hour]).mean().reset_index(level=0, drop=True)
        print(diurnal_daily)
        # Average across days
        diurnal_mean = diurnal_daily.groupby(level=0).mean()
        print(diurnal_mean)

        diurnal = noaa[noaa_species].groupby([noaa[noaa_species].index.floor("D"), noaa[noaa_species].index.hour]).mean().groupby(level=1).mean()
        print(diurnal)

        # hourly_diurnal = noaa[noaa_species].groupby(noaa.index.hour)
        # for name, group in hourly_diurnal:
        #     print(f"Group: {name}")
        #     print(group)
        #     print()
        mean_overall_diurnal = noaa[noaa_species].mean()

        fig_diurnal = plt.figure(figsize = (15,7), constrained_layout=True)
        widths = [15, 15]
        heights = [6, 2]
        spec_diurnal = fig_diurnal.add_gridspec(ncols=2, nrows=2, width_ratios=widths,
                                height_ratios=heights)
        
        ax1 = fig_diurnal.add_subplot(spec_diurnal[0,0])
        ax1.plot(hour_range, diurnal_mean, color='r', marker='.', label=f"Hrly Avg.")
        ax1.plot(hour_range, np.ones(len(hour_range))*mean_overall_diurnal, linestyle = 'dashed', color='r',label=f"Overall Avg.={mean_overall_diurnal:.2f} ppb")
        ax1.set_ylabel(f'{noaa_species}')
        ax1.set_xlabel('Hour (MDT)')
        ax1.set_xlim([0, 23])
        ax1.set_xticks(hour_range)
        # plt.title('')
        ax1.legend()
        
        # ax2 = fig_diurnal.add_subplot(spec_diurnal[1,0])
        # available_points = noaa[noaa_species].notna().groupby(noaa.index.hour).sum()
        # available_points.plot(ax=ax2, kind = 'bar', color = 'r')
        # ax2.plot(hour_range, np.ones(len(hour_range))*total_by_hour, linestyle = 'dashed', color='r',label=f"100%")
        # ax2.set_ylim([0, total_by_hour.max()+10])
        # ax2.set_yticks(np.arange(0,total_by_hour.max()+10,20))
        # ax2.set_xlabel('Hour (MDT)')
        # ax2.set_xlim([0, 23])
        # ax2.tick_params("x", rotation=0)
        # ax2.set_xticks(hour_range)

        # plt.show()

        # bihour_range = np.arange(0,24,2)

        # noaa['hour_2h'] = (noaa.index.hour // 2) * 2
        # bihourly_mean = noaa.groupby(['hour_2h'])[noaa_species].mean()
        # species_mean = noaa[noaa_species].mean()

        # plt.figure(figsize=(10,8))
        # plt.plot(bihour_range, bihourly_mean, marker='o', color = 'b')
        # plt.plot(bihour_range, np.ones(len(bihour_range))*species_mean, linestyle = 'dashed', color='b')
        # plt.xlabel('Hour (MDT)')
        # plt.ylabel(noaa_species)
        # plt.title('Bihourly Mean ' + noaa_species)
        # plt.grid(True)
        # plt.xticks(range(0,24))
        # plt.show()
def diurnal_fifteen_minute_avgs(interpolated_index, interpolated_col, vocname):
    if interpolated_col.isna().all():
        print(vocname, ' has no 15 min diurnal avg because it has no measurements.')
    
    else:
        #Extract time of day (HH:MM) as string for grouping
        time_of_day = interpolated_index.time

        #Show that the shape is roughly the same on all days
        expected_sample_count = 96          # 24h * 4 samples/hour
        minimum_coverage = 0.75       # 75% 
        #minimum counts needed to be considered a day with good coverage
        minimum_counts = expected_sample_count * minimum_coverage

        def good_coverage(day_series):
            # total non-NaN samples
            n_obs = day_series.notna().sum()

            if n_obs < minimum_counts:
                return False

            # check 6-hour blocks
            blocks = [(0, 6), (6, 12), (12, 18), (18, 24)]

            for start, end in blocks:
                block = day_series.between_time(f"{start:02d}:00", f"{end-1:02d}:59")
                if block.notna().sum() == 0:
                    return False

            return True
        good_coverage_day_flag = (interpolated_col.groupby(interpolated_index.normalize()).apply(good_coverage))
        good_coverage_days_indices = good_coverage_day_flag[good_coverage_day_flag].index
        series_good_coverage_days = interpolated_col[interpolated_index.normalize().isin(good_coverage_days_indices)]
        #print('Good coverage days: ', series_good_coverage_days)

        plt.figure(figsize=(12,8))
        series_good_coverage_days_july = series_good_coverage_days[series_good_coverage_days.index.month == 7]
        series_good_coverage_days_august = series_good_coverage_days[series_good_coverage_days.index.month == 8]

        for day, day_series in series_good_coverage_days_july.groupby(series_good_coverage_days_july.index.normalize()):
            print('day_series: ', day_series)
            each_day_overall_mean = day_series.mean()
            print('each_day_overall_mean: ', each_day_overall_mean)
            mean_15min_for_day_divided_by_overall_day_mean = day_series.divide(each_day_overall_mean)
            print('mean_15min_for_day_divided_by_overall_day_mean: ', mean_15min_for_day_divided_by_overall_day_mean)
            hour_and_minute_part_of_index = day_series.index.strftime("%H:%M")
            plt.plot(hour_and_minute_part_of_index, mean_15min_for_day_divided_by_overall_day_mean,  marker='o', label= str(day.strftime('%m-%d')))
            
        #Extract just the hour from index and use that for xticks
        hour_positions = [i for i, t in enumerate(hour_and_minute_part_of_index) if t.endswith(':00')]
        hour_labels = [str(int(hour_and_minute_part_of_index[i][:2])) for i in hour_positions]
        plt.ylim([0,9])
        plt.xticks(hour_positions, hour_labels)
        plt.legend(ncol = 3)

        plt.show()

        #Compute mean for 15-min time across all days
        diurnal_mean_15min = interpolated_col.groupby(interpolated_index.strftime("%H:%M")).mean()
        total_size = interpolated_col.groupby(interpolated_index.strftime("%H:%M")).size()

        index_oneday_15min = diurnal_mean_15min.index

        interpolated_mean_overall = interpolated_col.mean()

        fig_diurnal_15min = plt.figure(figsize = (15,7), constrained_layout=True)
        widths = [15]
        heights = [6, 2]
        spec_diurnal = fig_diurnal_15min.add_gridspec(ncols=1, nrows=2, width_ratios=widths,
                                height_ratios=heights)
        print('index_oneday_15min: ', index_oneday_15min)
        print('index_oneday_15min type: ', type(index_oneday_15min))
        print('index_oneday_15min values type: ', type(index_oneday_15min[0]))

        ax1 = fig_diurnal_15min.add_subplot(spec_diurnal[0,0])
        ax1.plot(index_oneday_15min, diurnal_mean_15min, color='c', marker='.', label=f"Hrly Avg.")
        ax1.plot(index_oneday_15min, np.ones(len(index_oneday_15min))*interpolated_mean_overall, linestyle = 'dashed', color='c',label=f"Overall Avg.={interpolated_mean_overall:.2f} ppb")
        ax1.set_ylabel(f'{vocname} (ppb)')
        ax1.set_xlabel('Hour (MDT)')

        #Extract just the hour from index and use that for xticks
        hour_positions = [i for i, t in enumerate(index_oneday_15min) if t.endswith(':00')]
        hour_labels = [str(int(index_oneday_15min[i][:2])) for i in hour_positions]
        ax1.set_xticks(hour_positions, hour_labels)
        ax1.legend()

        ax2 = fig_diurnal_15min.add_subplot(spec_diurnal[1,0])
        available_points = interpolated_col.notna().groupby(interpolated_index.strftime("%H:%M")).sum()
        available_points.plot(ax=ax2, kind = 'bar', color = 'r')
        ax2.tick_params("x", rotation=0)
        ax2.set_xticks(hour_positions, hour_labels)
        
        ax2.plot(index_oneday_15min, np.ones(len(index_oneday_15min))*total_size, linestyle = 'dashed', color='r',label='100%')
        ax2.annotate('100% (35/35)', xy = (1.01, 35), xycoords=('axes fraction', 'data'), va='center', ha='left')
        ax2.plot(index_oneday_15min, np.ones(len(index_oneday_15min))*total_size*0.75, linestyle = 'dashed', color='r',label='75%')
        ax2.annotate('75% (26.25/35)', xy = (1.01, 35*0.75), xycoords=('axes fraction', 'data'), va='center', ha='left')
        ax2.plot(index_oneday_15min, np.ones(len(index_oneday_15min))*total_size*0.50, linestyle = 'dashed', color='r',label='50%')
        ax2.annotate('50% (17.5/35)', xy = (1.01, 35*0.50), xycoords=('axes fraction', 'data'), va='center', ha='left')
        ax2.set_ylim([0, total_size.max()+10])
        ax2.set_yticks(np.arange(0,total_size.max()+10, 5))
        ax2.set_xlabel('Hour (MDT)')
        ax2.set_ylabel('Counts')
        # plt.show()

        plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/diurnal_15min/diurnal_15min_avg_after_interpolation_'+ str(vocname) + '.png', dpi =300)
        plt.show()

def diurnal_cycle_post_merge_check(df_index, init_col, filled_col, vocname):
    if 'ML_NoVar' in vocname:
        nan_series = pd.Series(index=hour_range)
        noaa_var = nan_series
        print('New diurnal cycle invalid for species ', vocname, 'due to not having ML data')
    
    else:
        # Compute daily diurnal cycles
        diurnal_daily = filled_col.groupby([df_index.date, df_index.hour]).mean().reset_index(level=0, drop=True)
        # print('diurnal_daily: ', diurnal_daily)

        # Average across days
        diurnal_mean = diurnal_daily.groupby(level=0).mean()
        diurnal_mean_size = diurnal_daily.groupby(level=0).size()
        # print(diurnal_mean)

        diurnal = filled_col.groupby([df_index.floor("D"), df_index.hour]).mean().groupby(level=1).mean()
        # print(diurnal)

        # hourly_diurnal = noaa[noaa_species].groupby(noaa.index.hour)
        # for name, group in hourly_diurnal:
        #     print(f"Group: {name}")
        #     print(group)
        #     print()
        mean_overall = filled_col.mean()

        fig_diurnal = plt.figure(figsize = (15,7), constrained_layout=True)
        widths = [15]
        heights = [6, 2]
        spec_diurnal = fig_diurnal.add_gridspec(ncols=1, nrows=2, width_ratios=widths,
                                height_ratios=heights)
        
        ax1 = fig_diurnal.add_subplot(spec_diurnal[0,0])
        ax1.plot(hour_range, diurnal_mean, color='c', marker='.', label=f"Hrly Avg.")
        ax1.plot(hour_range, np.ones(len(hour_range))*mean_overall, linestyle = 'dashed', color='c',label=f"Overall Avg.={mean_overall:.2f} ppb")
        ax1.set_ylabel(f'{vocname} (ppb)')
        ax1.set_xlabel('Time (MDT)')
        ax1.set_xlim([0, 23])
        ax1.set_xticks(hour_range)
        ax1.legend()

        ax2 = fig_diurnal.add_subplot(spec_diurnal[1,0])
        available_points = filled_col.notna().groupby(df_index.strftime("%H:%M")).sum()
        available_points.plot(ax=ax2, kind = 'bar', color = 'r')
        ax2.tick_params("x", rotation=0)
        ax2.set_xticks(hour_range)
        
        ax2.plot(hour_range, np.ones(len(hour_range))*diurnal_mean_size, linestyle = 'dashed', color='r',label='100%')
        ax2.annotate('100% (35/35)', xy = (1.01, 35), xycoords=('axes fraction', 'data'), va='center', ha='left')
        ax2.plot(hour_range, np.ones(len(hour_range))*diurnal_mean_size*0.75, linestyle = 'dashed', color='r',label='75%')
        ax2.annotate('75% (26.25/35)', xy = (1.01, 35*0.75), xycoords=('axes fraction', 'data'), va='center', ha='left')
        ax2.plot(hour_range, np.ones(len(hour_range))*diurnal_mean_size*0.50, linestyle = 'dashed', color='r',label='50%')
        ax2.annotate('50% (17.5/35)', xy = (1.01, 35*0.50), xycoords=('axes fraction', 'data'), va='center', ha='left')

        ax2.set_ylim([0, diurnal_mean_size.max()+10])
        ax2.set_yticks(np.arange(0,diurnal_mean_size.max()+10, 5))
        ax2.set_xlabel('Hour (MDT)')
        ax2.set_ylabel('Counts')

        plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/diurnal_vocs_after_filling_ml_with_udaq_plots/diurnal_avg_after_filling_ml_with_udaq_'+ str(vocname) + '.png', dpi =300)
        plt.show()