function med = spm_config_mediation(varargin)
    % Configuration file for mediation analysis
    %_______________________________________________________________________

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    images.type = 'files';
    images.name = 'Images';
    images.tag  = 'images';
    images.num  = [1 Inf];
    images.filter = 'image';
    images.help   = {'Data for each '};

    clusters.type = 'files';
    clusters.name = 'ROIs/clusters';
    clusters.tag  = 'clusters';
    clusters.num  = [1 Inf];
    clusters.filter = '.*cl.*\.mat';
    clusters.help   = {'Data for each '};

    expression.type = 'entry';
    expression.name = 'Matlab expression';
    expression.tag  = 'expression';
    expression.strtype = 'e';
    expression.num  = [1 Inf];
    expression.help   = {'Data for each '};

    X.type = 'choice';
    X.name = 'X';
    X.tag  = 'X';
    X.values = {images, clusters, expression};
    X.val  = {expression};
    X.help = {'Data for the X component. Only one of X, Y, or M may be a set of images.'};

    Y.type = 'choice';
    Y.name = 'Y';
    Y.tag  = 'Y';
    Y.values = {images, clusters, expression};
    Y.val  = {expression};
    Y.help = {'Data for the Y component. Only one of X, Y, or M may be a set of images.'};

    M.type = 'choice';
    M.name = 'M';
    M.tag  = 'M';
    M.values = {images, clusters, expression};
    M.val  = {images};
    M.help = {'Data for the M component. Only one of X, Y, or M may be a set of images.'};


    data.type = 'branch';
    data.name = 'Data';
    data.tag  = 'data';
    data.val  = {X, Y, M};
    data.check = @check_data;
    data.help = {'Specify the data sources for X, Y, and M. The data can be either manually entered as a Matlab expression, selected from a cluster/ROI file, or entered as a set of images.' ...
        '' 'By default, M3 searches over the brain for mediators, but both X and Y can be searched for as well.'};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Options
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    noop.type = 'const';
    noop.name = 'None';
    noop.tag  = 'noop';
    noop.val  = {0};
    noop.help   = {'None'};

    mask.type = 'files';
    mask.name = 'Mask';
    mask.tag  = 'mask';
    mask.val = {''};
    mask.num  = [0 1];
    mask.filter = 'image';
    mask.help   = {'Mask image'};

    robust.type = 'menu';
    robust.name = 'Robust regression';
    robust.tag  = 'robust';
    robust.val  = {0};
    robust.values = {1, 0};
    robust.labels = {'Yes', 'No'};
    robust.help   = {'Uses robust regression to downweight outliers.'};

    weighting.type = 'menu';
    weighting.name = 'Weighting';
    weighting.tag  = 'weighting';
    weighting.val  = {0};
    weighting.values = {1, 0};
    weighting.labels = {'Yes', 'No'};
    weighting.help   = {'Weight second-level based on first-level variance.'};

    shift.type = 'entry';
    shift.name = 'Shifted correlations';
    shift.tag  = 'shift';
    shift.val  = {[-6 6]};
    shift.num = [1 2];
    shift.strtype = 'i';
    shift.help   = {'Shifts timeseries +/- the specified units, and then uses best fit.', '',...
        'Default is -6 to 6.'};

    latent.type = 'menu';
    latent.name = 'Latent activity';
    latent.tag  = 'latent';
    latent.val  = {0};
    latent.values = {1, 0};
    latent.labels = {'Yes', 'No'};
    latent.help   = {'Computes HRF and deconvolves it from timeseries in order to use latent, neural activity instead of BOLD activity.'};

    timeseries_options.type = 'choice';
    timeseries_options.name = 'Timeseries options';
    timeseries_options.tag  = 'timeseries_options';
    timeseries_options.values = {noop shift latent};
    timeseries_options.val = {noop};

    bootsamples.type = 'entry';
    bootsamples.name = 'Number of bootstrap iterations';
    bootsamples.tag  = 'bootsamples';
    bootsamples.val  = {1000};
    bootsamples.num = [1 1] ;
    bootsamples.strtype = 'n';
    bootsamples.help   = {'Number of bootstrap iterations to run.', '',...
        'Default is 1000.'};

    boot.type = 'branch';
    boot.name = 'Bootstrap top level';
    boot.tag  = 'boottop';
    boot.val  = {bootsamples};
    boot.help   = {'Uses bootstrapping on top level.'};

    bootfirst.type = 'branch';
    bootfirst.name = 'Bootstrap';
    bootfirst.tag  = 'bootfirst';
    bootfirst.val  = {bootsamples};
    bootfirst.help   = {'Uses bootstrapping on bottom level.'};

    signperm.type = 'menu';
    signperm.name = 'Sign permutation';
    signperm.tag  = 'signperm';
    signperm.val  = {0};
    signperm.values = {1, 0};
    signperm.labels = {'Yes', 'No'};
    signperm.help   = {'Uses sign permutation test on second level.'};

    multi_boot_perm_options.type = 'choice';
    multi_boot_perm_options.name = 'Bootstrapping/permutation options';
    multi_boot_perm_options.tag  = 'boot_perm_options';
    multi_boot_perm_options.values = {noop boot signperm};
    multi_boot_perm_options.val = {noop};
    multi_boot_perm_options.help = {'Choose bootstrapping on the bottom level, bootstrapping at the top-level, sign permutation testing at the top level, or ''None''. ', '',...
        'If ''None'' is chosen, ordinary least squares (OLS) will be used.'};

    single_boot_perm_options.type = 'choice';
    single_boot_perm_options.name = 'Bootstrapping options';
    single_boot_perm_options.tag  = 'boot_perm_options';
    single_boot_perm_options.values = {noop boot};
    single_boot_perm_options.val = {noop};
    single_boot_perm_options.help = {'Choose bootstrapping or ''None''.', '',...
        'If ''None'' is chosen, ordinary least squares (OLS) will be used.'};


    ARorder.type = 'entry';
    ARorder.name = 'Auto-regression order';
    ARorder.tag  = 'ARorder';
    ARorder.val  = {0};
    ARorder.num  = [1 1];
    ARorder.strtype = 'w';
    ARorder.help   = {'Auto-regressive noise model on first level data. An AR of 0 is equivalent to no AR noise model.'};

    covs.type = 'entry';
    covs.name = 'Covariates';
    covs.tag  = 'covs';
    covs.strtype = 'e';
    covs.num  = [1 Inf];
    covs.val  = {[]};
    covs.help   = {'Covariates to be controlled for in all regressions.'};

    singleoptions.type = 'branch';
    singleoptions.name = 'Options';
    singleoptions.tag  = 'options';
    singleoptions.val  = {mask, robust, single_boot_perm_options, covs};

    multioptions.type = 'branch';
    multioptions.name = 'Options';
    multioptions.tag  = 'options';
    multioptions.val  = {mask, robust, weighting, ARorder, timeseries_options, multi_boot_perm_options, covs};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Overall
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    singlelevel.type = 'branch';
    singlelevel.name = 'Single-level mediation';
    singlelevel.tag  = 'singlelevel';
    singlelevel.val = {data, singleoptions};
    singlelevel.prog = @run_mediation;
    singlelevel.help = {'Single-level mediation'};

    multilevel.type = 'branch';
    multilevel.name = 'Multi-level mediation';
    multilevel.tag  = 'multilevel';
    multilevel.val = {data, multioptions};
    multilevel.prog = @run_mediation;
    multilevel.help = {'Multi-level mediation'};

    med.type = 'repeat';
    med.name = 'Mediation';
    med.tag  = 'mediation';
    med.values = {singlelevel, multilevel};
    med.num  = [1 Inf];
    med.help = {'Add single-level and multi-level mediation jobs here.'};
    med.modality = {'PET','FMRI'};
end

function result = check_data(data)
    result = [];
    [X, Y, M] = get_XYM_data(data);

    if(isnumeric(X) && ~isvector(X))
        result = 'X is not a vector';
    elseif(isnumeric(Y) && ~isvector(Y))
        result = 'Y is not a vector';
    elseif(isnumeric(M) && ~isvector(M))
        result = 'M is not a vector';
    elseif(input_length(X) ~= input_length(Y) || input_length(Y) ~= input_length(M))
        result = 'Mismatch in length between X, Y and M.';
    end
end

function l = input_length(in)
    if(ischar(in))
        l = size(in, 1);
    else
        l = length(in);
    end
end


function run_mediation(job)
    [X, Y, M] = get_XYM_data(job.data);
    
    args = {};

    if(~isempty(job.options.mask{1}))
        args{end+1} = 'mask';
        args{end+1} = job.options.mask{1};
    end

    if(~isempty(job.options.covs))
        args{end+1} = 'covs';
        args{end+1} = job.options.covs;
    end

    if(job.options.robust)
        args{end+1} = 'robust';
    else
        args{end+1} = 'norobust';
    end

    if(job.options.ARorder > 0)
        args{end+1} = 'arorder';
        args{end+1} = job.options.ARorder;
    end

    if(isfield(job.options.timeseries_options, 'shift'))
        args{end+1} = 'shiftrange';
        args{end+1} = job.options.timeseries_options.shift;
    end

    if(isfield(job.options.timeseries_options, 'latent'))
        args{end+1} = 'latent';
    end

    if(isfield(job.options.boot_perm_options, 'boottop'))
        args{end+1} = 'boottop';
    end

    if(isfield(job.options.boot_perm_options, 'bootfirst'))
        args{end+1} = 'bootstrapfirst';
    end

    if(isfield(job.options.boot_perm_options, 'signperm'))
        args{end+1} = 'signperm';
    end

    mediation_brain(X, Y, M, args{:});
end

function [X, Y, M] = get_XYM_data(data)
    X = data.X.(char(fieldnames(data.X)));
    Y = data.Y.(char(fieldnames(data.Y)));
    M = data.M.(char(fieldnames(data.M)));

    if iscellstr(M), M = char(M); end
    if iscellstr(X), X = char(X); end
    if iscellstr(Y), Y = char(Y); end
end


