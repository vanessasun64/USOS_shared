Current Workflow



1\. Got large EPA VOC file from UDAQ (Bart), called Verbose. Includes measurements for VOC species in many different locations in Utah.

2\. Downloaded CSV file containing EPA's Parameters Codes, species, units, CAS Numbers, etc. in epa\_parameters\_official.csv



3\. Function hawthorne\_usos\_extract\_species in mapper\_functions matches species in the Verbose file to the EPA Parameter codes to give us an excel spreadsheet that specifies which species we have measurements from UDAQ. It is saved as Hawthorne\_EPA\_match.xlsx and .csv of the same name. We now have a CAS Number for each species and the Parameter code.



4\. Used the function get\_chem\_info\_epa on this Hawthorne\_EPA\_match.csv file. This matches the CAS Number from PubChem to give us additional information, specifically we'll have columns for:

 	- EPA Species Name

 	- query\_results: confirms that our species was found by CID match

 	- iupac\_name

 	- synonyms

 	- cid

 	- molecular weight

 	- chemical formula

 	- inchi

 	- inchi\_key

 	- SMILES

 	- cas\_number

File is saved as EPA\_Hawthorne\_pubchem\_match.xlsx and EPA\_Hawthorne\_pubchem\_match.csv



Manually separated m/p xylene.



Re-saved CSV file.



5\. Copied the InChI keys given from the excel file EPA\_Hawthorne\_pubchem\_match.xlsx and pasted into EPA Comptox Batch Search at: https://comptox.epa.gov/dashboard/batch-search

Checked off the Input Type as InChIKey

Click CHOOSE EXPORT OPTIONS

Click Excel under CHOOSE EXPORT FORMAT

Select the following:



**Chemical Identifiers**

DTXSID

Chemical Name

InChIKey

IUPAC Name



**Structures**

SMILES

InChI String



**Intrinsic and Predicted Properties**

Molecular Formula

OPERA Model Predictions



Click DOWNLOAD EXPORT FILE



Output file is named CCD-Batch-Search\_Year-Month-Day\_hour\_minute\_second.xlsx where hour, min, and second are likely in UTC. Renamed to comptox\_batch\_search\_all.xlsx



Opened the file in Excel. There should be one sheet called Main Data. This is the one we want.



Deleted the unnecessary newly added species:

ROW 5: Hydrocarbons, C7-C8, n-alkanes

ROW 15: C10-12 Isoalkanes

ROW 17: Alkanes, C16-20-iso-



Saved as CSV and Excel.



6\. Run function comptox\_extract. This will give us an excel and CSV file currently called EPA\_CRACMM\_mapped that gets the kOH and calculates a log10 C Star from the Comptox files we saved in the previous step. Then, it uses the CRACMM Mapper (needs SMILES, kOH, and log10C Star for each species) to tell us what species it thinks is the equivalent in CRACMM.



Manually, I added colors to the excel spreadsheet to indicate which ones didn't match my initial mapping by hand and may need to ask Havala about. Resaved as the same file name.





-------

GEOS-Chem mapped by database spreadsheet provided by Jessica, matching by inchi using new function in mapping\_functions.py, then manual inputs by Vanessa. Saved to EPA\_CRACMM\_GEOSCHEM\_mapped.csv and excel.

