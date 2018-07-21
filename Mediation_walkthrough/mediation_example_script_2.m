
% mediation_example_script1 and 2 do the same analysis.
% ...script1 is very terse, and includes the essential commands only.
% ...script2 is longer and includes more checking that files are available, etc.

% Go to the directory this script is stored in and run the script from
% there.

disp('Mediation example script: Creating a new analysis subfolder and running')
disp('Sample mediation analysis here:');

basedir = pwd;
disp(basedir);

%% 1. Locate the sample data directory and files with sample behavioral data
% ---------------------------------------------------------------------
% Goal: Make sure we have folders for the data and mask on the Matlab path

% Locate pre-cooked results (and behavioral data in Mediation_SETUP.mat):
dir_to_find = 'mediation_Example_Data_Wager2008_Msearch_R_XisRIFGstim_norobust';

dinf = what(dir_to_find);
if isempty(dinf)
    disp('CANNOT FIND Pre-cooked mediation results')
    fprintf('Looking for: %s\n', dir_to_find);
    disp('Install Canlab_Core_Tools and Mediation Toolbox from canlab.github.com');
    disp('Run canlab_toolbox_setup or add toolboxes with subfolders to the Matlab path');
    return
end
beh_file_path = dinf.path;

behavior_file = fullfile(beh_file_path, 'mediation_SETUP.mat');

% Locate image data in Sample Datasets:
dir_to_find = 'CanlabCore/Sample_datasets/Wager_et_al_2008_Neuron_EmotionReg';

dinf = what(dir_to_find);
if isempty(dinf)
    disp('CANNOT FIND sample image data')
    fprintf('Looking for: %s\n', dir_to_find);
    disp('Install Canlab_Core_Tools and Mediation Toolbox from canlab.github.com');
    disp('Run canlab_toolbox_setup or add toolboxes with subfolders to the Matlab path');
    return
end
image_file_path = dinf.path;

% Optional: Specify a mask for the analysis
% You can run without this,but it's helpful.
% You can also use a mask to constrain multiple comparisons correction when
% you get results:

mask_name = which('gray_matter_mask.img');

if isempty(mask_name)
    disp('CANNOT FIND mask image')
    fprintf('Looking for: %s\n', 'gray_matter_mask.img');
    disp('Install Canlab_Core_Tools and Mediation Toolbox from canlab.github.com');
    disp('Run canlab_toolbox_setup or add toolboxes with subfolders to the Matlab path');
    return
end

% Alternative:  specify mask file from saved mediation output
% mask_name = fullfile(beh_file_path, 'mask.img');


%% 2. Locate the images and list file names
% ---------------------------------------------------------------------
% Goal: Generate a list of the image filenames you want to use in the mediation

% Load an object with the emotion regulation data:
% save the file names in img_names
data_obj = load_image_set('emotionreg');
img_names = data_obj.fullpath;  

% Alternative: List image names; may not work on all file systems
% img_names = filenames(fullfile(image_file_path, 'con*img'), 'char');


%% 3. Create a new test directory to run mediation
% ---------------------------------------------------------------------

% make a new analysis directory
andir = 'Test_mediation';
mkdir(andir)
cd(andir)

%% 4a. Run the mediation with OLS (quick):
% ---------------------------------------------------------------------
% Fast version: No bootstrapping. Good for a preliminary check:

results = mediation_brain(SETUP.X, SETUP.Y, img_names, 'names', {'IFG' 'ReappSuccess' 'BrainMediator'}, 'mask', mask_name);

%% 4b. Run the mediation with bootstrapping:
% ---------------------------------------------------------------------
% Slow version: Good for final results. For FINAL results, we recommend
% 10,000 bootstrap samples

results = mediation_brain(SETUP.X, SETUP.Y, img_names, 'names', {'IFG' 'ReappSuccess' 'BrainMediator'}, 'mask', mask_name, 'boot', 'pvals', 5, 'bootsamples', 1000);

%% 5. Get results and publish an HTML report:
% ---------------------------------------------------------------------

publish_mediation_report

%% 6. Get and save results (older format but nore complete):
% ---------------------------------------------------------------------

mediation_brain_results_all_script;

