% Runs batch analyses and publishes HTML report with figures and stats to
% results/published_output in local study-specific analysis directory.

% Run this from the main mediation results directory (basedir)

close all
warning off, clear all, warning on 

resultsdir = pwd;

fprintf('Creating HTML report for results in:\n%s\n', resultsdir);

% ------------------------------------------------------------------------
% Check that we have a valid mediation results dir
% ------------------------------------------------------------------------

is_mediation_dir = exist(fullfile(resultsdir, 'mediation_SETUP.mat'));
is_run = exist(fullfile(resultsdir, 'X-M-Y_effect.img'));

if ~is_mediation_dir
    fprintf('%s\nis not a mediation results directory because mediation_SETUP.mat is missing.\nSkipping report.\n', resultsdir);
    return
end

if ~is_run
    fprintf('Mediation does not appear to have run correctly because X-M-Y_effect.img is missing.\nSkipping report.\n');
    return
end

% ------------------------------------------------------------------------
% Set HTML report filename and options
% ------------------------------------------------------------------------

pubdir = fullfile(resultsdir, 'published_output');
if ~exist(pubdir, 'dir'), mkdir(pubdir), end

pubfilename = ['mediation_brain_report_' scn_get_datetime];

p = struct('useNewFigure', false, 'maxHeight', 800, 'maxWidth', 1600, ...
    'format', 'html', 'outputDir', fullfile(pubdir, pubfilename), 'showCode', false);

% ------------------------------------------------------------------------
% Run and report status
% ------------------------------------------------------------------------

publish('mediation_brain_results_report.m', p)

myhtmlfile = fullfile(pubdir, pubfilename, 'mediation_brain_results_report.html');

if exist(myhtmlfile, 'file')
    
    fprintf('Saved HTML report:\n%s\n', myhtmlfile);
    
    web(myhtmlfile);
    
else
    
    fprintf('Failed to create HTML report:\n%s\n', myhtmlfile);
    
end

