### Methods
* UDAQ: Hourly measurements at Hawthorne provided by Bart, removed values > 120 ppb, marked as unusable (NOTE: Need to check with Bart on color-coding for other quality controlled points, maybe say something about checking). Initially recorded in MST, shifted 1 hour to match MDT.
* ML: Mobile Lab data is 1 hour merged.
* Instrument is 2B Tech 106-L, utilized by both UDAQ and the Mobile Lab. Uncertainty is 1.5 ppb or 2%. https://2btech.io/items/ambient-ozone-monitors/model-106-l-ozone-monitor/

### Ozone comparison between UDAQ & ML with uncertainties
`hawthorne_udaq_mobilelab_o3_comparison_july_aug_uncertainties.png`

### Alternatives: 
* `hawthorne_udaq_mobilelab_o3_comparison_july_aug_with_inset_and_uncertainties.png`
* `hawthorne_udaq_mobilelab_o3_comparison_july_aug_no_uncertainties.png`
* `hawthorne_udaq_mobilelab_o3_comparison_first_week_with_grid_and_uncertainties.png`

**NOTES:** The uncertainty is so small that it doesn't show very clearly on the figure. The shaded uncertainty in the legend looks awkward (but can add it)


---
### MDA8 Ozone
`hawthorne_udaq_mobilelab_mda8_exceedances_o3_comparison.png`

UDAQ data shows 8 exceedance days during 34 days of the campaign, Mobile Lab data shows 9 exceedance days. The extra date that the Mobile Lab observes as an exceedance day is 07/28. However, due to holes in the Mobile Lab data (missing/invalid data), we should have more exceedance days than Mobile Lab shows.

We are trying to use a modified UDAQ O3 concentration to "fill" the holes in the Mobile Lab data during the F0AM runs. This should also create a new MDA8 ozone estimate, where we see more exceedances. For example, 07/30 has a hole of several hours for the Mobile Lab data. If we were to fill the Mobile Lab data on that date with modified UDAQ O3 concentrations, it should be an exceedance day that's both missed by the UDAQ and Mobile Lab instruments independently. 

---
**NOTES:** Negative indicates an underestimate.

### Mean Bias 
`hawthorne_udaq_mobilelab_o3_comparison_meanbias.png`

* also sometimes called Mean Bias Error to highlight that it is an error-based metric, think of as an absolute difference

* $\text{Mean bias} = \frac{1}{n} \sum_{i=1}^n (M_i - O_i) = \overline{M} - \overline{O}$

Overall average shows that UDAQ's instrument is underestimating ozone by 2.27 ppb compared to the Mobile Lab's instrument. When looking at the hourly mean bias, there are a few hours in which the UDAQ instrument overestimates ozone (9 AM to 12 PM). Morning hours of 8, 9, 11, and 12 have mean biases within the instrument's uncertainty range of 1.5 ppb, while all other hours show a more significant error. The largest mean bias occurs at Hour 20, where the UDAQ instrument underestimates by 6.82 ppb of ozone. 


### Mean Normalized Bias
`hawthorne_udaq_mobilelab_o3_comparison_mean_normalized_bias_percentage.png`

* think of as a relative difference
* $\text{MNB} = \frac{1}{n} \sum_{i=1}^n \left( \frac{M_i - O_i}{O_i} \right)$

Overall average shows that UDAQ's instrument is underestimating ozone by 4.45% compared to the Mobile Lab's instrument. Only the 12 PM hourly mean normalized bias falls within the instrument's uncertainty range of 2%, with all other hours showing a more significant error. In the morning hours of 9AM-12 PM, the UDAQ instrument overestimates ozone, with all other hours showing a consistent underestimation. The hours of 6 and 7 AM, as well as 8 PM, show the largest underestimates in MNB, where their error reaches over 10% and the most significant error underestimates by 12.25%.

---
### Average Diurnal Cycle of Ozone

* `hawthorne_udaq_mobilelab_o3_comparison_hrly_mean.png` - means for UDAQ and ML only
* `hawthorne_udaq_mobilelab_o3_comparison_hrly_mean_with_std.png` - means for UDAQ and ML with standard deviations shaded
* `hawthorne_udaq_mobilelab_o3_comparison_hrly_median_with_percentiles.png` - median with 25th and 75th percentile dotted (I didn't like the shading but not 100% happy with this either)

### Extra Statistics:
The UDAQ instrument overestimates ozone 22.3% of the time, with only 11.0% of the difference being larger than the instrument uncertainty. 

---
### VOC Speciation
Mobile Lab: `mobilelab_voc_speciation.png`
* 1 hr merges; with Isoprene, Benzene, Toluene, and Styrene from the PTR, Formaldehyde from Picarro, and the rest from iWAS
* Each slice represents the fractional contribution to total VOCs
* Formaldehyde generally has the largest contribution except isoprene starts to have a more significant contribution to total VOCs starting 9 PM
