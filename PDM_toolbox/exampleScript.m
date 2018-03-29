% This example scripts demonstrates the multivariate mediation analysis
% (PDM) for single-trial data as described in: 
% 
% Chén OY, Crainiceanu C, Ogburn EL, Caffo BS, Wager TD, Lindquist MA (2017) 
% High-dimensional multivariate mediation with application to neuroimaging 
% data. Biostatistics 
% Available at: https://academic.oup.com/biostatistics/article/doi/10.1093/biostatistics/kxx027/3868977/High-dimensional-multivariate-mediation-with [Accessed August 8, 2017].
% 
% Geuter S, Losin EAR, Roy M, Atlas LY, Schmidt L, Krishnan A, Koban L, 
% Wager TD, Lindquist MA (under review) Beyond the pain matrix - a 
% comprehensive analysis of brain mediators of pain.
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Because we can't share single-trial data on github you will need to
% create your own test data. You need an treatment variabele (X),
% single-trial brain images (M), and an outcome variable (Y). The example
% below assumes that each subject has their data in a separate exp*.mat
% file in the sub-directory 'data'. Each exp[subjectID].mat file has three
% variables named 'treatment', 'outcome', and 'mediator'. Alternatively,
% you can specify the data directly in 3 cell arrays (X,Y,M) and skip the
% first section below.
% 
%
%
% You will need SPM12 and the CanLabCore tools 
% (https://github.com/canlab/CanlabCore) installed
%

clear; close all;

% % load your data
% F = filenames('data/exp*mat');
% 
% for k=1:numel(F)
%    
%     load(F{k});
%     
%     xx{k} = treatment;
%     yy{k} = outcome;
%     mm{k} = mediator;
%         
% end
% 
% disp('min trials:');
% disp(min(cellfun(@numel,xx)));


%% Reduce the dimensionality of the brain-mediator data using PVD

% project onto lower dimensional space keeping th max number of components
pdm = multivariateMediation(xx,yy,mm,'noPDMestimation');

% same as above, but keep only 25 components
pdm = multivariateMediation(xx,yy,mm,'noPDMestimation','B',25);


%% Compute the multivariate brain mediators

% use previous PVD data and compute 5 PDMs [default]
pdm = multivariateMediation(pdm);


% compute 10 PDMs but not the jointPDM
pdm = multivariateMediation(pdm,'nPDM',10,'jpdm',false);


% redo the PVD with 20 components and compute 5 PDMs + jointPDM
pdm = multivariateMediation(xx,yy,mm,'B',20);



%% bootstrap voxel weights for significance

% bootstrap the first PDM with 100 samples
pdm = multivariateMediation(pdm,'bootPDM',1,'Bsamp',100);


% do PVD, compute 5 PDMs, and bootstrap weights for the first two PDMs 
% with 500 samples [use 5,000 or more samples for real studies]. 
% Save results to file 'PDMresults.mat'
pdm = multivariateMediation(xx,yy,mm,'B',20,'nPDM',5,'bootPDM',1:2,'Bsamp',500,'save2file','PDMresults.mat');


%% visualize the results
close all;

% visualize PDM1
mask = 'brainmask.nii'; % mask used for data extraction. needs to match the voxels in mm{i} and pdm.Wfull
dat = fmri_data(mask,mask,'noverbose');

% threshold PDM1 at p<0.01
dat.dat = pdm.Wfull{1}.*(pdm.boot.p{1}<0.01);

% display thresholded PDM
canlab_results_fmridisplay(dat);
orthviews(dat);


% and PDM2
dat2=dat;
dat2.dat = pdm.Wfull{2}.*(pdm.boot.p{2}<0.01);

% display thresholded PDM
figure(3);
canlab_results_fmridisplay(dat2);


