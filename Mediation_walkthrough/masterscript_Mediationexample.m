close all
clear all

%This script is explained in detail in the powerpoint titled
%Mediation_sample_data_walkthrough

%If first time running, go to folder of images and change
%con_00810028.img.txt to just con_00810028.img

%Step 1 of walkthrough

%cd('/Users/Erik/Desktop/Fall 2015/Principals of fMRI/fMRI_matlab/Mediation_example_data')
%Just drag the folder into the command window to get to the proper directory

%Step 2
imgs = filenames('Images_for_mediation/con_*img', 'char', 'absolute');

%Step 3
behav_dat=importdata('X_Y_data.txt')

%Step 4
mask=filenames('Images_for_mediation/gray_matter_mask.img','char','absolute')
mask(100:end)
canlab_results_fmridisplay(mask, 'compact2');

%Step 5
andir='Test_mediation_2017';
mkdir(andir)
cd(andir)

%Step 6
scn_map_image(mask,deblank(imgs(1,:)), 'write', 'resliced_mask.img');
ls

%Step 7 
x=behav_dat.data(:,1);
y=behav_dat.data(:,2);
names={'X:RVLPFC' 'Y:Reappraisal_Success' 'M:BrainMediator'};
results=mediation_brain(x,y,imgs,'names',names,'mask','resliced_mask.img','boot','pvals',5); 

%Step 8
    %Make yourself a cup of tea while the results are compiled as this is going to take a while

%Step 9

%Not sure what functions to run here? 

