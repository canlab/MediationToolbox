% Runs batch analyses and publishes HTML report with figures and stats to
% results/published_output in local study-specific analysis directory.

% Adapted from publish_mediation_report.m by Lukas Van Oudenhove
% May 2021
% NOTE: only tested on single-level multivariate mediation for now!
% pdm calculated using 
% https://github.com/labgas/proj-emosymp/blob/main/secondlevel/model_1_CANlab_classic_GLM/emosymp_m1_s6_mediation_NPS.m
% see also example script
% https://github.com/canlab/MediationToolbox/blob/master/PDM_toolbox/Multivariate_Mediation_ExampleScript.m

% Run this from the main mediation results directory (basedir)

close all
warning off, clear all, warning on 

resultsdir = pwd;
contrastname = split(resultsdir,"\");
contrastname = contrastname{end-1};
Yname = split(resultsdir,"\");
Yname = Yname{end};

fprintf('Creating HTML report for results in:\n%s\n', resultsdir);

% ------------------------------------------------------------------------
% Check that we have a valid multivariate mediation results dir
% ------------------------------------------------------------------------

is_mediation_dir = exist(fullfile(resultsdir, strcat('PDMresults_',contrastname,'_',Yname,'.mat')));

if ~is_mediation_dir
    fprintf('%s\nis not a mediation results directory because mediation_SETUP.mat is missing.\nSkipping report.\n', resultsdir);
    return
end

% ------------------------------------------------------------------------
% Set HTML report filename and options
% ------------------------------------------------------------------------

pubdir = fullfile(resultsdir, 'published_output');
if ~exist(pubdir, 'dir'), mkdir(pubdir), end

pubfilename = [strcat('multivariate_mediation_brain_report_',contrastname,'_',Yname,'_') scn_get_datetime];

p = struct('useNewFigure', false, 'maxHeight', 800, 'maxWidth', 1600, ...
    'format', 'html', 'outputDir', fullfile(pubdir, pubfilename), 'showCode', true);

% ------------------------------------------------------------------------
% Run and report status
% ------------------------------------------------------------------------

publish('multivariate_mediation_brain_results_report.m', p)

myhtmlfile = fullfile(pubdir, pubfilename, 'multivariate_mediation_brain_results_report.html');

if exist(myhtmlfile, 'file')
    
    fprintf('Saved HTML report:\n%s\n', myhtmlfile);
    
    web(myhtmlfile);
    
else
    
    fprintf('Failed to create HTML report:\n%s\n', myhtmlfile);
    
end

