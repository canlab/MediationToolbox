% This is older and not complete
% see mediation_brain_results_all_script.m

error('Do not run this whole script; it''s just pieces of code.  See mediation_brain_results_all_script.m');

%% init results 

% [clpos, clneg, clpos_data, clneg_data] = mediation_brain_results('all','prune','overlay',EXPT.overlay, 'thresh', [.005 .01 .05], 'size', [3 10 10]);
ovl = '/Users/tor/Documents/Tor_Documents/CurrentExperiments/Lab_Emotion/Pre-appraisal/Spontaneous_vs_voluntary_reg/Analyses_and_images_June2007/meanwT1.img';

heightthresh = [.005 .01 .05];
sizethresh = [1 1 1];
whichtype = 'all';
prunestr = 'prune';
dosave = 1;

[clpos, clneg, clpos_data, clneg_data] = mediation_brain_results(whichtype,prunestr,'overlay',ovl, 'thresh', heightthresh, 'size', sizethresh);

load mediation_SETUP
% get data from peak clusters
% % clpos{1} = tor_extract_rois(SETUP.M, clpos{1});
% % clneg{1} = tor_extract_rois(SETUP.M, clneg{1});

% for X search
clpos{1} = tor_extract_rois(SETUP.X, clpos{1});
clneg{1} = tor_extract_rois(SETUP.X, clneg{1});


% overall size threshold
overallsize = 10;
s = cat(1, clpos_data.numVox); clpos_data(s < overallsize) = [];
s = cat(1, clneg_data.numVox); clneg_data(s < overallsize) = [];



doranks = strcmp(SETUP.rank_data, 'Yes');

%%
diary('mediation_network_level_output.txt');
fprintf('Height threshold(s) are: '); fprintf('%3.4f\t', heightthresh); fprintf('\n');
fprintf('Extent threshold(s) are: '); fprintf('%3.0f\t', sizethresh); fprintf('\n');
fprintf('Pruning string is: %s\n ', prunestr); 
fprintf('Blobs are based on effect: %s\n ', whichtype); 
fprintf('Rank data: %s\n ', SETUP.rank_data);

ynstr = {'No' 'Yes'};
fprintf('Saving results clusters in med_results.mat: %s\n ', ynstr{dosave + 1});

diary off
%% name contiguous clusters 
% use peak voxels at highest thresholds from regions!
% % clpos_peak_cl = cluster_names(clpos{1}, 1);
% % clneg_peak_cl = cluster_names(clneg{1}, 1);

clpos_data = cluster_names(clpos_data, 1);
clneg_data = cluster_names(clneg_data, 1);

%% make table


diary('mediation_network_level_output.txt');
[clpos_data, clneg_data] = mediation_brain_print_tables(clpos_data, clneg_data, doranks);
diary off

if dosave, save med_results -append cl*
end
%% figures of blobs

% show blobs for figures, no text
savename = sprintf('mediators_%3.3f_%3.3f_%3.3f_%02d_%02d_%02d_', heightthresh, sizethresh);
savename(savename == '.') = '-';

h = cluster_orthviews_showcenters([clpos_data clneg_data], 'saggital', ovl, 0, 0, 'k');
scn_export_papersetup(600);
saveas(gcf,[savename 'sagg'],'png')
h = cluster_orthviews_showcenters([clpos_data clneg_data], 'axial', ovl, 0, 0, 'k');
scn_export_papersetup(600);
saveas(gcf,[savename 'ax'],'png')
h = cluster_orthviews_showcenters([clpos_data clneg_data], 'coronal', ovl, 0, 0, 'k');
scn_export_papersetup(600);
saveas(gcf,[savename 'cor'],'png')

figh = findobj('Tag', 'Graphics'); figure(figh);
scn_export_papersetup(600);

for i = 1:length(clpos_data)
    spm_orthviews('Reposition', clpos_data(i).mm_center);
    saveas(gcf, [savename 'orth' num2str(i)], 'png');
end

for i = 1:length(clneg_data)
    spm_orthviews('Reposition', clneg_data(i).mm_center);
    saveas(gcf, [savename 'orth' num2str(length(clpos_data) + i)], 'png');
end

mkdir slice_png_images
!mv mediators_*png slice_png_images

%% other correlation plots, network analysis
% save clusters in med_results
load mediation_SETUP
mediation_nmds_results


%% surface figs

mediation_brain_surface_figs(clpos, clneg)
scn_export_papersetup(800);

%saveas(gcf,'mediation_surf_lat_all','png');
saveas(gcf,'mediation_surf_lat_all','fig');
