# Akerman Lab Open Ephys pipeline

Akerman lab in-house Matlab pipeline for multi-unit and sorted spike data processing. Data is acquired using NeuroNexus 16-channel or 32-channel linear acute probes in S1, under different types of whisker and optogenetic stimulation.

The pipeline starts with the scripts in the 'Preprocessing' folder, which will save a series of .mat files with an ephys_data variable, in which extracted spike times are organised by condition, by channel and by trial. Conditions are automatically determined from trial synchronisation data contained in ADC channels 1-4 from the Open Ephys A2D board. 


# Data preprocessing...

## ...for multiunit / local field potential analysis:

1) Collect electrophysiology data using Open Ephys software and Akerman lab PulsePal stimulation protocols

2) Enter experiment and recording metadata in a metadata spreadsheet (excel); see /Metadata for a template of how to organise this

3) Run 'preprocess_multiunit' function (see 'help preprocess multiunit' for instructions on how to run). This will extract spike times and LFPs segments from each recording aligned with the trials, and organise the data by condition. The resulting preprocessed data files will be stored in a folder structure by experiment type and by date.

4) Visualise the data using individual plots in /Plotting

5) Visualise multiple plots and extract metrics from individual experiment data using functions in /TBC

6) Extract data from multiple experiments to look at entire population results or compare between populations using functions in /TBC



## ...for spike sorting using Kilosort:

1) Collect data using OpenEphys system and Akerman lab PulsePal stimulation protocols

2) Enter experiment and recording metadata in a metadata spreadsheet (excel); see /Metadata for a template of how to organise this

3) Use run_Kilosort_Akermanlab to prepare data for Kilosort and run Kilosort (this includes copying the relevant Kilosort files to a dedicated folder, )

4) Curate data using Phy

5) Use sync_curated_data function to distribute the sorted unit data into a folder structure organised by experiment type, by date

6) The data should now be compatible with all the same visualisation and analysis functions specified in steps 4), 5) and 6) of the multiunit analysis pipeline, except that data is organised by unit rather than by channel.


# Data visualisation

/Plotting has functions to easily make raster plots, post-stimulus time histograms, heat maps of spiking activity by channel over time, etc:

*raster_plot:* 		A raster plot that can be organised either by trial or by channel
*psth:* 		Post stimulus time histogram
*spike_density_plot:* 	Heat map of spiking activity by channel over time


# Data analysis

...

***...to be continued...***





Joram van Rheede