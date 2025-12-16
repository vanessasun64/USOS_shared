Use CTRL+Shift+V to show this Markdown file with formatting.

---

This file describes Vanessa's process to obtain info on all VOCs measured at Hawthorne and begin mapping them to CRACMM, GEOS-Chem, and CB6r5h chemical mechanisms. It covers some of the files in USOS\_shared/Hawthorne\_data/mappings

---

**1\.** Got large EPA VOC file from UDAQ (Bart), called Verbose. Includes measurements for VOC species in many different locations in Utah.


**2\.** Downloaded CSV file containing EPA's Parameters Codes, species, units, CAS Numbers, etc. in `/given_data/epa_parameters_official.csv`


**3\.** Function `hawthorne_usos_extract_species` in `scripts/mapper_functions.py` matches species in the Verbose file to the EPA Parameter codes to give us an excel spreadsheet that specifies which species we have measurements from UDAQ. It is saved as `script_output/Hawthorne_EPA_parameters_match.xlsx` and .csv of the same name. We now have a CAS Number for each species and the Parameter code.


**4\.** Used the function `get_chem_info_epa` in this `script_output/Hawthorne_EPA_match.csv` file. This matches the CAS Number from PubChem to give us additional information, specifically we'll have columns for:

* EPA Species Name
* query\_results: confirms that our species was found by CID match
* iupac\_name
* synonyms
* cid
* molecular weight
* chemical formula
* inchi
* inchi\_key
* SMILES
* cas\_number

File is saved as `script_output/EPA_Hawthorne_pubchem_match.xlsx` and csv

Manually separated m/p xylene, resaved as `manually_edited/EPA_Hawthorne_pubchem_match_updated.xlsx` and csv



**5\.** Copied the InChI keys given from the excel file `manually_edited/EPA_Hawthorne_pubchem_match_updated.xlsx` and pasted into EPA Comptox Batch Search at: https://comptox.epa.gov/dashboard/batch-search

Checked off the Input Type as InChIKey

Click CHOOSE EXPORT OPTIONS

Click Excel under CHOOSE EXPORT FORMAT

Select the following:

* **Chemical Identifiers**
    * DTXSID
    * Chemical Name
    * InChIKey
    * IUPAC Name

* **Structures**
    * SMILES
    * InChI String

* **Intrinsic and Predicted Properties**
    * Molecular Formula
    * OPERA Model Predictions


Click DOWNLOAD EXPORT FILE


Output file is named `given_data/CCD-Batch-Search_Year-Month-Day_hour_minute_second.xlsx` where hour, min, and second are likely in UTC.  Opened the file in Excel. There should be one sheet called Main Data. This is the one we want. Deleted the unnecessary newly added species:

* ROW 5: Hydrocarbons, C7-C8, n-alkanes
* ROW 15: C10-12 Isoalkanes
* ROW 17: Alkanes, C16-20-iso-

Saved and renamed to `manually_edited/comptox_batch_search_all.xlsx` as CSV and Excel.

**6\.** Run function `comptox_extract` in `scripts/mapper_functions.py`. This will give us an excel and CSV file called `script_output/EPA_CRACMM_mapped` that gets the kOH and calculates a log10 C Star from the Comptox files we saved in the previous step. Then, it uses the CRACMM Mapper (needs SMILES, kOH, and log10C Star for each species) to tell us what species it thinks is the equivalent in CRACMM.

Manually, I added colors to the excel spreadsheet to indicate which ones didn't match my initial mapping by hand and may need to **ask Havala about.** Resaved as `manually_edited/EPA_CRACMM_mapped_updated.xlsx`


**7\.** GEOS-Chem species info is provided by Jessica, stored in `F0AM-4.3.0.1/Chem/GEOSChem/GC_database_final\_UPDATED.csv`
geoschem_mech_mapping function in `scripts/mapping_functions.py` matches the species measured by UDAQ by InChI to GEOS-Chem species. Saves output as `script_output/UDAQ_Hawthorne_CRACMM_GEOSCHEM_mapped.xlsx` and `script_output/UDAQ_Hawthorne_CRACMM_GEOSCHEM_mapped.csv`

**8\.** Made manual edits to `script_output/UDAQ_Hawthorne_CRACMM_GEOSCHEM_mapped.xlsx` by Vanessa:
* Jessica Haskins manually helped Vanessa fill the GEOS-Chem mappings
* Kelvin Bates helped to clarify some of the TMB mappings that Jessica didn't know
* m-xylene and p-xylene were recombined into one row (Row 45)
    * different values separated by semicolon delimiter (;) and no spaces
    * only one value used if characteristic is the same between the two species
* CB6r5h Mappings added, one column is labeled as INCORRECT and another labeled as maybe
    * INCORRECT represents Vanessa's initial attempt to map species but Unreactive Carbon (UNR) isn't actually a species in this version of CB
    * maybe represents Vanessa's subsequent attempt to map species with accurate carbon allocations. **Needs review by Greg Yarwood**
* USOS Mapping added, representing the corresponding Mobile Lab measurement for each species
* UDAQ_Variable added, representing the variable that will be used when exporting to F0AM.

Resaved to `UDAQ_Hawthorne_CRACMM_GEOSCHEM_CB6r5h_mapped_updated_date.xlsx` and csv, where date is the most recent date of update.