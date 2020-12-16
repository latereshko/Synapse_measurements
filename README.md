# Synapse_measurements with 'metamorph_synapse_G_R_B_colocalize.JNL'
Quantify intensity of fluorescence at colocalized puncta (synaptic markers) on GFP+ neurons with Metamorph journal.
'metamorph_synapse_G_R_B_colocalize.JNL' separates channels, uses granularity app to identify colocalization between channels, and logs intensity of all channels within ROI drawn. Final steps are optional to measure dendrite lengths within the ROI.
Clean data & plot synapse measurements made in Metamorph with R .

# Synapse_intensity_metamorph_to_R
Use this script to get plot intensity values from Metamorph summary data.
Fill in treatment/condition labels based on the order in which images were analysed. The script is designed to pull in image names from directory folder.
*Use this before the script to get calculate density per dendrite length

# Density_per_dendrite_lengths_metamorph_to_R
Use this script directly after synapse_intensity_metamorph_to_R to sum dendritic length measurements per cell and calculate "synapse density" by dividing "Count"/sum dendrite lengths
