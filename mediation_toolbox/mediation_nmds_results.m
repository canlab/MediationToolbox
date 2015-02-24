% see mediation_brain_results_script_all for setup details

x = SETUP.X;
y = SETUP.Y;
yname = SETUP.names(2);
covs = SETUP.covariates;

clpos_peak_cl = clpos_data;
clneg_peak_cl = clneg_data;

%% name contiguous clusters
% use peak voxels at highest thresholds from regions!
% % clpos_peak_cl = cluster_names(clpos{1}, 1);
% % clneg_peak_cl = cluster_names(clneg{1}, 1);

if ~isempty(clneg_peak_cl)
    nnames = cat(2, {clneg_peak_cl.shorttitle});
    datn = cat(2, clneg_peak_cl(:).timeseries);
else
    nnames = [];
    datn = [];
end

if ~isempty(clpos_peak_cl)
    pnames = cat(2, {clpos_peak_cl.shorttitle});
    datp = cat(2, clpos_peak_cl(:).timeseries);
else
    pnames = [];
    datp = [];
end

if doranks
    y = rankdata(y);
    datp = rankdata(datp);
    datn = rankdata(datn);
    x = rankdata(x);
end

%% get average data



diary('mediation_network_level_output.txt');

% correlation matrix
disp('Partial correlations between regions and Y, controlling for X');
disp(['Ranks: ' SETUP.rank_data])
disp('_____________________________________________________________');

[partialrs, p] = partialcorr([y datp datn], x);
correlation_to_text(partialrs, p < .05, [yname pnames nnames]);

disp('_____________________________________________________________');
fprintf(1,'\n\n');

diary off
%% stepwise regressions


diary('mediation_network_level_output.txt');
disp('Stepwise regressions (Including X, but not controlling for it.');

STEP = stepwise_tor([x datp datn], y, [SETUP.names(1) pnames nnames]);
% % STEP = stepwise_tor([x mean(datp,2) mean(datn,2)], y, [SETUP.names(1) {'Pos' 'Neg'}]);

diary off

%% nmds

diary('mediation_network_level_output.txt');

disp('Multidimensional scaling output');
disp('_____________________________________________________________');

cl = [clpos_peak_cl clneg_peak_cl];
if ~isfield(SETUP, 'covariates'), SETUP.covariates = []; end
c = []; c.covs_nointerest = [x covs]; c.outcome = y;
for i = 1:length(cl), c.names{i} = cl(i).shorttitle; end
c = nmdsfig_tools('get_data_matrix',c,cl,'timeseries',1,[],doranks);
c = nmdsfig_tools('get_correlations',c);
[c.GroupSpace,c.obs,c.implied_dissim] = shepardplot(c.D,[]);

c = nmdsfig_tools('cluster_solution',c, c.GroupSpace, 2:5, 1000, []);
c = nmdsfig_tools('apply_clusters',c);

disp('_____________________________________________________________');
fprintf(1,'\n\n');

diary off

if dosave, save med_results -append c
end

%%
% ***************
% Change c.ClusterSolution.classes to plot in diff colors

%           1      2    3   4
%colors = {'yo'  'bo' 'go' 'ro' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};

% c_complete.ClusterSolution.classes = [3 3 1 1 1 1 1 1 1 2 2]';
% nmdsfig_tools('nmdsfig_plot', c, 0, 0, 'nofill');

% **********************************************************************
%%% MODIFY CLASSES AND COLORS HERE TO CHANGE colors of blobs and spheres
% :  classes(end-1:end) = 2
% **********************************************************************

classes = c.ClusterSolution.classes;
colors = {'yo' 'bo' 'go' 'ro' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};


%% NMDS plot
nmdsfig_tools('nmdsfig_plot',c, 0, 0, 0)
if dosave, saveas(gcf,'NMDS_flat_plot_orig_nooutcome','png'); end

[c.direct_mtx, c.mediated_mtx, c.mediators] = matrix_direct_effects(c.STATS.sigmat2, c.dat);
nmdsfig_tools('removelines');
nmdsfig_tools('drawlines',c.GroupSpace, c.direct_mtx);

if dosave, saveas(gcf,'NMDS_flat_plot_direct_nooutcome','png'); end

%% NMDS plot: With outcome
disp('Creating c_complete structure and plotting with outcomes and X (mediated)');
c_complete = nmdsfig_tools('nmdsfig_plot_withcovs',c, cl, doranks);

if dosave
    save med_results -append c_complete
    saveas(gcf,'NMDS_flat_plot_direct_withoutcomes','png');
    saveas(gcf,'NMDS_flat_plot_direct_withoutcomes','fig');
end


%% NMDS plot with OUTCOMES and with your chosen colors (re-assign classes above)

create_figure('nmdsfig');
nmdsfig(c_complete.GroupSpace,'classes', c_complete.ClusterSolution.classes, 'names', c_complete.names,'sig', ...
    c_complete.direct_mtx, ...
    'colors', {'yo' 'bo' 'go'}, 'sizes', 16);

h = findobj(gcf,'Type','Text')
set(h, 'FontSize', 24)
if dosave
    scn_export_papersetup(800);
    saveas(gcf,['NMDS_complete_flat_plot_colors'],'png')
end

%% glass brain plots
% SEPARATE CLASSES

classes = c.ClusterSolution.classes;
colors = {'yo' 'bo' 'go'};
type = 'blobs';     % blobs or spheres
sigmatrix = c.STATS.sigmat;  % or c.direct_mtx

for i = 1:max(classes)
    %subplot(2, 2, i);

    wh = classes == i;
    mycl = cl(wh);

    cluster_nmdsfig_glassbrain(mycl, ones(length(mycl), 1), colors(i), sigmatrix(wh, wh), [], type); %'samefig');
    material dull, lightRestoreSingle(gca), lighting gouraud, %camlight left

    if dosave
        scn_export_papersetup(600);
        saveas(gcf,['mediation_nmds_glass_' type '_class' num2str(i)],'png')
        saveas(gcf,['mediation_nmds_glass_' type '_class' num2str(i)],'fig')
    end
end

%% glass brain plots: ALL CLASSES TOGETHER, WITH SEED

if exist('seed_cl', 'var')

    % add cluster for X seed onto regions
    cl = [clpos_peak_cl clneg_peak_cl];
    cl2 = merge_clusters(cl, seed_cl);
    cl2 = cl2([end 1:end-1]);
    cl2(1).timeseries = SETUP.X;

    c2 = []; c2.covs_nointerest = [SETUP.covariates]; c.outcome = SETUP.Y;
    for i = 1:length(cl2), c2.names{i} = cl2(i).shorttitle; end

    c2 = nmdsfig_tools('get_data_matrix',c2,cl2,'timeseries',1,[],doranks);
    c2 = nmdsfig_tools('get_correlations',c2);


    %[c2.GroupSpace,c2.obs,c2.implied_dissim] = shepardplot(c2.D,[]);

    c2.ClusterSolution.classes = [3 ; ones(length(clpos_peak_cl), 1); 2 * ones(length(clneg_peak_cl), 1)];


    %       cluster_orthviews_classes(cl,c2.ClusterSolution.classes, EXPT.overlay, 'saggital', 0);
    %       nmdsfig_tools('nmdsfig_plot',c2, 0, 0, 'fill');
    %  c2_complete = nmdsfig_tools('nmdsfig_plot_withcovs',c2, cl2, doranks)

    [c2.direct_mtx, c2.mediated_mtx, c2.mediators] = matrix_direct_effects(c2.STATS.sigmat2, c2.dat);


    cluster_nmdsfig_glassbrain(cl2, c2.ClusterSolution.classes, colors, c2.direct_mtx, []);
    if dosave
        scn_export_papersetup(600);
        saveas(gcf,['mediation_nmds_glass_all_withseed' num2str(i)],'png')
    end


    % now with direct only
    %cluster_nmdsfig_glassbrain(cl, classes, colors, c.direct_mtx, []);

end

%% predictions of behavior
diary('mediation_network_level_output.txt');

% NOTE: THIS MAY BE DIFFERENT FROM ACTUALLY INCLUDING X IN THE REGRESSION!

if max(classes) > 1
    disp('Stepwise regressions for CLASSES, controlling for X:')
    disp('_____________________________________________________________');

    c = nmdsfig_tools('predict_behavior',c,'classes');
else
    disp('Stepwise regressions for REGIONS, controlling for X:')
    disp('_____________________________________________________________');

    c = nmdsfig_tools('predict_behavior',c,'regions');
end
fprintf(1,'\n\n');

f = findobj('Tag', 'class_scatterplots');
if ishandle(f) && dosave
    
    texth = findobj(f, 'Type', 'Text'); set(texth, 'FontSize', 24);
    axh = get(f,'Children');
    for hh = 1:length(axh), set(axh(hh), 'FontSize', 24); end

    disp('Saving class_scatterplots.fig');
    scn_export_papersetup; saveas(gcf,'class_scatterplots','fig');
    saveas(gcf,'class_scatterplots','png');
end

diary off

%% save
if dosave
    disp('Appending results and clusters to med_results.mat')
    save med_results -append clpos clneg clpos_peak_cl clneg_peak_cl c

end

fprintf(1,'\n\n');

%% slices of regions

cluster_orthviews_classes(cl,c.ClusterSolution.classes, EXPT.overlay, 'saggital', 0);
if dosave
    scn_export_papersetup(600)
    saveas(gcf, 'cluster_orthviews_classes_sagg', 'png');
end

%% mediation class results

diary('mediation_network_level_output.txt');

disp('Mediation results at class (network) level')
disp('Controlling for other class avgs')
disp('_____________________________________________________________');
fprintf(1,'\n\n');

% go back to original data (not adjusted data in c struct)
myclasses = c.ClusterSolution.classes;
mydat = [datp datn];

for i = 1:max(c.APPLY_CLUSTER.classes)

    wh = myclasses == i;
    c.raw_class_avg_dat(:, i) = nanmean(mydat(:, wh)')';

end

for i = 1:max(c.APPLY_CLUSTER.classes)
    wh = find(myclasses == i);

    fprintf('Class %3.0f, Regions =', i)
    fprintf(' %3.0f', wh);
    fprintf(1,'\n----------------------------------------\n');
    covdat = c.raw_class_avg_dat;
    covdat(:,i) = [];

    if doranks
        [paths, stats] = mediation(rankdata(SETUP.X), rankdata(SETUP.Y), rankdata(c.raw_class_avg_dat(:, i)), 'M', rankdata(covdat), 'plots', 'verbose', 'boot', 'bootsamples', 10000, 'save');
    else
        [paths, stats] = mediation(SETUP.X, SETUP.Y, c.raw_class_avg_dat(:, i), 'M', covdat, 'plots', 'verbose', 'boot', 'bootsamples', 10000, 'save');
    end

        
    fprintf(1,'\n----------------------------------------\n');
end

disp('_____________________________________________________________');
fprintf(1,'\n\n');

diary off
