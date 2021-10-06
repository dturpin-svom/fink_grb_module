# fink_grb_module
A project to develop a science module in Fink devoted to the direct search for LSST/ZTF transient candidates with known GRBs


1/ Two GRB catalogs can be loaded (currently manually updated) to perform the direct search : 

	- all_fermi-gbm.xml extracted from the Fermi GBM GRB catalog : https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3table.pl?tablehead=name%3Dfermigbrst&Action=More+Options.
  
	- all_swift_2019_2021.csv extracted from the Swift GRB catalog : https://swift.gsfc.nasa.gov/archive/grb_table/
	
2/ A grb_config.json file is provided to modify the different key parameters used to do the analysis
NOTE: The Fermi GRB210204270 is chosen as an example as it is associated with ZTF21aagwbjr discovered by the ZTFRest pipeline (https://github.com/growth-astro/ztfrest) and later robustly identified as the afterglow counterpart of GRB210204270. 

3/ A Jupyter notebook is provided as a template to run the different steps of the search for identifying a credible association between a chosen GRB and some ZTF candidates
