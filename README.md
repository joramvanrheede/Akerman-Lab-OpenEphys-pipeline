Akerman-Lab-OpenEphys-pipeline

Akerman lab in-house Matlab pipeline for OpenEphys data processing

The pipeline starts with the scripts in the 'Preprocessing' folder, which will save a series of .mat files with an ephys_data variable, in which extracted spike times are organised by condition / trial type, by channel and by trial, 


Data preprocessing...

...for multiunit / local field potential analysis:

1) Collect data using OpenEphys system and Akerman lab PulsePal stimulation protocols

2) Enter experiment metadata in a metadata spreadsheet (excel)

3) Run 'preprocess_multiunit' function (see 'help preprocess multiunit' for instructions on how to run). This will 



...for spike sorting using Kilosort:

1) Collect data using OpenEphys system and Akerman lab PulsePal stimulation protocols

2) Enter experiment metadata in a metadata spreadsheet (excel)

3) Use run_Kilosort_Akermanlab to prepare data for Kilosort and run Kilosort

4) curate data using Phy

5) Use sync_curated_data


Data visualisation

/Plotting has functions for visualising post-stimulus time histograms,


Data analysis





