function mediation_region_stepwise(clpos_data, clneg_data, SETUP, doranks)
%
% Single-level only.  Stepwise regression of clusters predicting outcome.

    diary mediation_region_stepwise.txt


    cl = [clpos_data clneg_data];
    if ~isfield(SETUP, 'covariates'), SETUP.covariates = []; end
    c = []; c.covs_nointerest = [SETUP.X SETUP.covariates]; c.outcome = SETUP.Y;
    for i = 1:length(cl), c.names{i} = cl(i).shorttitle; end
    c = nmdsfig_tools('get_data_matrix',c,cl,'timeseries',1,[],doranks);
    c = nmdsfig_tools('get_correlations',c);


    c = nmdsfig_tools('predict_behavior',c,'regions');
%%

    wh = find(c.indiv_STEPWISE.inmodel);
    nregions = length(wh);

    disp('Mediation results for regions')
    disp('Controlling for other class avgs')
    disp('_____________________________________________________________');
    fprintf(1,'\n\n');

    % go back to original data (not adjusted data in c struct)
    mydat = cat(2, cl.timeseries);
    if doranks
        for i = 1:size(mydat, 2), mydat(:,i) = rankdata(mydat(:,i)); end
    end

    mydat = mydat(:, wh);

    mynames = c.names(wh);

    for i = 1:nregions

        covdat = mydat;
        covdat(:,i) = [];
        covnames = mynames;
        covnames(i) = [];

        fprintf('Mediating Region: %s\n', c.names{wh(i)})
        fprintf('Controlling for: ');
        for j = 1:length(covnames), fprintf('%s ', covnames{j}); end

        fprintf(1,'\n----------------------------------------\n');


        [paths, stats] = mediation(SETUP.X, SETUP.Y, mydat(:, i), 'covs', covdat, 'plots', 'verbose', 'boot', 'bootsamples', 10000, 'names', [SETUP.names(1:2) c.names(wh(i))]);

        %[paths, stats] = mediation(SETUP.X, SETUP.Y, mydat(:, 1), 'covs', [], 'plots', 'verbose', 'boot', 'bootsamples', 10000, 'names', [SETUP.names(1:2) cl(1).shorttitle]);

        fprintf(1,'\n----------------------------------------\n');
    end

    disp('_____________________________________________________________');
    fprintf(1,'\n\n');

    diary off

end