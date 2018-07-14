% Runs batch analyses and publishes HTML report with figures and stats to 
% results/published_output in local study-specific analysis directory.

% Run this from the main mediation results directory (basedir)

close all
clear all

resultsdir = pwd;

pubdir = fullfile(resultsdir, 'published_output');
if ~exist(pubdir, 'dir'), mkdir(pubdir), end

% ------------------------------------------------------------------------
pubfilename = ['mediation_brain_report_' scn_get_datetime];

p = struct('useNewFigure', false, 'maxHeight', 800, 'maxWidth', 1600, ...
    'format', 'html', 'outputDir', fullfile(pubdir, pubfilename), 'showCode', false);

publish('mediation_brain_results_report.m', p)