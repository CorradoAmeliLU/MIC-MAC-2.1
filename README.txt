=============================================================================================================================================
Immunofluorescence staining protocol and large‐scale analysis to quantify microglial cell morphology in mouse brain at single-cell resolution
                   Frida Lind-Holm Mogensen, Corrado Ameli, Alexander Skupin, Alessandro Michelucci
=============================================================================================================================================

The files for the tutorial are available at: https://owncloud.lcsb.uni.lu/s/qoXE9SmgIZ34xvX/authenticate
The password is: micmac

Extract the Library.tar file before using the Matlab scripts.

The file tutorial_raw.ims is the raw acquisition. Open this file in Imaris and perform the segmentation as described. For this image we chose a threshold for the background subtraction of 0.8, and a filtering between 18.000 voxels and 200.000 voxels.

The file tutorial_segmented.ims is the output from Imaris Segmentation. This file should be placed into ./Tutorial/Code/segmented to be processed by the Matlab script FeatureExtract.m
This script will return two files:
- _ims_attributes.mat -> A table containing all the morphological features for each cell
- _ims_volumes.mat -> A structure containing the 3d volumes of each cell (in the same order as the _ims_attributes.mat file)

In order to unify the datasets coming from each sample, you can use the script Merge.m.
This script will return four files:
-data.txt -> contains the matrix with the numerical morphological values
-feat_names.txt -> contains the names of the morphological features (or name of the columns of data.txt)
-samples_attribute_name -> contains the names of the samples merged
-structsPerSample -> contains the number of cells present in each sample (sorted by samples_attribute_name)

Note that the ordering of the cells/samples is preserved throughout all the files.
The nth cell of the first sample will be in position n in the data.txt file.
The nth cell of the second sample will be in position (number of cells of sample 1) + n in the data.txt file.


