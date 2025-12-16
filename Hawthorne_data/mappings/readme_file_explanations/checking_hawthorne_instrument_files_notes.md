Use CTRL+Shift+V to show this Markdown file with formatting.

---

This file describes the process for getting extra instrument data for the species UDAQ is measuring at Hawthorne.

---

1. Copied hawthorne-only instrument data from `given_data/AMP500_2321577-0.txt` to `manually_edited/hawthorne_all_instruments_udaq.txt`

2. Pulled Parameter Codes and Method Codes from the MM at Hawthorne from `hawthorne_all_instruments_udaq.txt`
Put into `hawthorne_instruments_with_parameter_code_and_method_code.csv`
The excel file is the same as the CSV but extra sheets as scrap for how I filtered all the Parameter and Method codes.

3. Was given a file by Bart called `methods_all_with_cnumber2 - Copy.csv` , renamed to `given_data/all_instrument_methods.csv`

4. Wrote code in `scripts/hawthorne_instruments_match_script.py` to match Parameter Codes and Method Codes from `hawthorne_instruments_with_parameter_code_and_method_code.csv` with `all_instrument_methods.csv`
and make new dataframe, copying over the Parameter Code, Species, Method Code, info about instrument used to take measurement, and units for all species measured at Hawthorne by UDAQ.

This is saved into `script_output/hawthorne_all_species_measured_by_udaq_instruments_info.csv`

From my Verbose file, only missing alpha and beta pinene from VOC Hawthorne UDAQ measurements. UDAQ also measured carbon monoxide, NOx (compare to Mobile Lab data?), sulfur dioxide (not measured by Mobile Lab) that we may able to use with F0AM.