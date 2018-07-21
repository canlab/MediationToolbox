% med_results = mediation_brain(X, Y, M, ['thresholds', thresholds], ['mask', maskname], ['names', names], ...
%   ['robust'|'norobust'], ['signperm'], ['arorder', arorder], ['boot'|'noboot'], ['multilevel'|'summarystats'], ['covs', covariates], ...
%   ['covnames', covnames])
%
% This function does a robust/nonrobust search for mediators between two vars X and Y
%
% For each contrast, it saves con, t, p, filtered t imgs, and clusters.mat
% NB: p-value images are 2-tailed
%
% Argument list:
% ==============================================================
% data options
% -----------------------------------
% 'mask', maskname = varargin{i+1};
% 'names', names = varargin{i+1};
% 'thresholds', thresholds = varargin{i+1};
% {'covs', 'covariates', 'mediation_covariates'}, mediation_covariates = varargin{i+1};
% 
% analysis choice options
% -----------------------------------
% 'robust', robust_opt = 'robust';
% 'norobust', robust_opt = 'norobust';
% 'rank', dorank = 1;
% {'boot', 'boot1'}, boot_opt = 'boot1';
% 'noboot', boot_opt = 'noboot';
% 
% multilevel options (if data is multilevel)
% -----------------------------------
% 'signperm', boot_opt = 'signperm';
% {'multilevel', 'hierarchical'}, multilev_opt = 'hierarchical';
% 'summarystats', multilev_opt = 'summarystats';
%
% process control options
% -----------------------------------
% 'start' or 'startslice' : followed by slice number to start with (for
% restarts)
%
% E.g.:
% Do robust regression search over matrix M for mediators
% ---------------------------------------------------------------------
% names = {'Lneg-Lneu_corr_pca1' 'Neg. emotion report' 'Neg - Neu Stimulation'}
% X = EXPT.cov(:,2); Y = EXPT.cov(:,1); M = EXPT.SNPM.P{6};
% mask = 'rob_tmap_0002_filt_t_2-76_k10_both.img';
% med_results = mediation_brain(X, Y, M, 'names', names, 'mask', mask, 'robust');
%
% Do search within anticipation-responsive regions for effects linked to
% stim period mediator (from pca)
% ---------------------------------------------------------------------
% names = {'Neg - Neu Anticip' 'Neg. emotion report' 'Lneg-Lneu_corr_pca1' }
% M = EXPT.cov(:,2); Y = EXPT.cov(:,1); X = EXPT.SNPM.P{8};
% mask = '../robust0099/rob_tmap_0001_filt_t_2-76_k10_both.img';
% med_results = mediation_brain(X, Y, M, 'names', names, 'mask', mask);
%
% SEE CODE FOR MORE INPUT OPTIONS: ** THESE SHOULD BE DOCUMENTED
% results = mediation_brain(X, Y, M, 'names', names, 'mask', mask, 'rank', 'startslice', 2); save mediation_results_struct results
%
% Load an existing structure and re-run, searching for mediators
% ---------------------------------------------------------------------
% load mediation_SETUP
% med_results = mediation_brain(SETUP.X, SETUP.Y, SETUP.M, 'names', ...
% SETUP.names, 'covs', SETUP.covariates, 'mask', SETUP.mask, 'boot', 'pvals', 5);

%Programmers' notes
% ----------------------------
% Created by Tor Wager, Feb 2007
% Modified Dec 2009 : Add functionality to save paths for additional
% mediators

function med_results = mediation_brain(X, Y, M, varargin)
    % --------------------------------------
    % * Defaults
    % --------------------------------------

    spm_get_defaults('modality','fmri');
    
    robust_opt = 'norobust';
    boot_opt = 'boot1';
    multilev_opt = 'hierarchical';
    arorder = 0;
    names = {'X' 'Y' 'M'};
    thresholds = .05;
    dorank = 0;
    startslice = 1;
    mediation_covariates = [];
    covnames = {};
    additionalM = [];
    
    % --------------------------------------
    % * Setup
    % --------------------------------------

    % Get which type of search to do based on which input is an image set
    % E.g., 'Search for mediators'; 'Search for indirect influences';  'Search for mediated outcomes';
    [cmdstring, maskname] = get_cmdstring(X, Y, M);

    
    % parse inputs
    for i = 1:length(varargin)
        arg = varargin{i};
        if ischar(arg)
            switch lower(arg)
                case 'mask', maskname = varargin{i+1};
                case 'names', names = varargin{i+1};
                case {'thr', 'thresh', 'thresholds'}, thresholds = varargin{i+1};
                case {'covs', 'covariates', 'mediation_covariates'}, mediation_covariates = varargin{i+1};
    
                case {'covnames', 'covariate names'}, covnames = varargin{i+1};
                    
                % Analysis choice options
                case 'robust', robust_opt = 'robust';
                case 'norobust', robust_opt = 'norobust';
                case 'rank', dorank = 1;
                case {'boot', 'boot1'}, boot_opt = 'boot1';
                case 'noboot', boot_opt = 'noboot';
                    
                % Multilevel options (if data is multilevel)    
                case 'signperm', boot_opt = 'signperm';
                case {'multilevel', 'hierarchical'}, multilev_opt = 'hierarchical';
                case 'summarystats', multilev_opt = 'summarystats';
                    
                % Timeseries options
                case 'arorder', arorder = varargin{i+1};
                case {'start', 'startslice'}, startslice = varargin{i+1}; load med_results
                
                 % Additional mediators, other optional inputs
                case {'m', 'M'}, additionalM = varargin{i+1}; % just need number, to get output size; that's it
            end
        end
    end

    % Number of additional mediators
    if iscell(additionalM)
        num_additionalM = size(additionalM{1}, 2); 
    else
        num_additionalM = size(additionalM, 2);
    end
    
    ynstr = {'No' 'Yes'};
    fprintf('Mask: %s\n', maskname);
    fprintf('Rank data: %s\n', ynstr{dorank + 1});
    fprintf('\n')
    
    
    % --------------------------------------
    % * Get flags for which imgs to save
    %   Initialize SETUP structure for saving options
    % --------------------------------------
    med_flags = get_med_flags(cmdstring);
    
    SETUP = struct('cmdstring', cmdstring, 'names', {names}, 'mask', maskname, 'X', X, 'Y', Y, 'M', M, 'robust_opt', robust_opt, ...
        'boot_opt', boot_opt, 'multilev_opt', multilev_opt, 'rank_data', ynstr{dorank + 1}, 'covariates', mediation_covariates, ...
        'covnames', {covnames}, 'arorder', arorder);

    % --------------------------------------
    % * Load data
    % --------------------------------------
    str = display_string('Reading data.');

    dolegacy = false;
    if ~exist('fmri_data.m', 'file') || ~isa(fmri_data(), 'fmri_data')
        dolegacy = true;
    end
    
    if ~dolegacy
        % We have object-oriented tools installed.
        % Object-oriented version.  Now reslices mask automatically.
        
        mask_obj = fmri_data(maskname, 'noverbose');

        switch cmdstring
            case 'Search for mediators'
                data_obj = fmri_data(M, maskname, 'noverbose');
                M = double(data_obj.dat');
                
            case 'Search for indirect influences'
                data_obj = fmri_data(X, maskname, 'noverbose');
                X = double(data_obj.dat');
                
            case 'Search for mediated outcomes'
                data_obj = fmri_data(Y, maskname, 'noverbose');
                Y = double(data_obj.dat');
            otherwise
                error('Unknown cmd string "%s".', cmdstring);
        end
        
        % resample mask to make sure in same space
        % 
                
        mask_obj = resample_space(mask_obj, data_obj, 'noverbose');
        maskInfo = mask_obj.volInfo;  % uses n_inmask, xyzlist
        
    else
        % Legacy version: Does not use object-oriented tools.
        % No object-oriented tools/CANlab core installed. May cause other
        % problems
        warning('CANlab Core Tools not installed. We recommend you install from https://github.com/canlab');
        
        switch cmdstring
            case 'Search for mediators'
                [M, maskInfo] = iimg_get_data(maskname, M);
            case 'Search for indirect influences'
                [X, maskInfo] = iimg_get_data(maskname, X);
            case 'Search for mediated outcomes'
                [Y, maskInfo] = iimg_get_data(maskname, Y);
            otherwise
                error('Unknown cmd string "%s".', cmdstring);
        end

    end % Object-oriented or legacy
    
    SETUP.maskInfo = maskInfo;
    save mediation_SETUP SETUP

    erase_string(str);

    % --------------------------------------
    % * Initialize Output
    % --------------------------------------
    npaths = 5 + 3 * num_additionalM;
    
    med_results.paths = NaN .* zeros(maskInfo.n_inmask, npaths);
    med_results.ste = NaN .* zeros(maskInfo.n_inmask, npaths);
    med_results.pvals = NaN .* zeros(maskInfo.n_inmask, npaths);
 
   % maskInfo.dim(4) = 16;   % float

    if startslice > 1
        str = display_string('Using existing output images.');
    else
        str = display_string('Initializing output images.');
        write_output(med_results, names, maskInfo, med_flags);
    end

    erase_string(str);

    % --------------------------------------
    % * Rank all data, if requested
    % --------------------------------------
    if dorank
        str = display_string('Ranking data for nonparametric tests.');
        for i = 1:size(X, 2), X(:,i) = rankdata(X(:,i)); end
        for i = 1:size(Y, 2), Y(:,i) = rankdata(Y(:,i)); end
        for i = 1:size(M, 2), M(:,i) = rankdata(M(:,i)); end
        erase_string(str);
    end
    
    % --------------------------------------
    % * Computation
    % --------------------------------------
    %   1 a   X -> M relationship
    %   2 b   M -> Y relationship
    %   3 cp  unmediated X -> Y relationship (residual)
    %   4 c   X -> Y relationship
    %   5 ab  mediated X -> Y by M (a * b)

    if dolegacy

        z = maskInfo.xyzlist(:,3);
        
    else
        
        z = data_obj.volInfo.xyzlist(:, 3); % for in-mask voxels only
        
%         z = double(maskInfo.image_indx);
%         z(maskInfo.wh_inmask) = maskInfo.xyzlist(:, 3);
        
    end

    nslices = max(z);
     
    str = display_string('Statistics.');
    spm_progress_bar('init');

    for i = startslice:nslices
        % calculates and writes images for this slice
        process_slice();

        med_results;

        save med_results med_results
        spm_progress_bar('set', i/nslices);
    end

    spm_progress_bar('clear');
    erase_string(str);



    % END MAIN FUNCTION CODE


    % -------------------------------------------------------------------
    %
    %
    % INLINE FUNCTIONS
    %
    %
    %
    % -------------------------------------------------------------------


    % --------------------------------------
    % * Do one slice
    % --------------------------------------
    function process_slice()
        % i is slice
        whvox = (z == i);

        fprintf('\n\n===============================================\n\nSlice %3.0f : %3.0f voxels in-mask', i, sum(whvox));
        t1 = clock;

        switch cmdstring
            case 'Search for mediators'
                slice_results = mediation_search('M', X, Y, M(:,whvox), varargin{:});

                %slice_results = mediation_search('M', X, Y, M(:,whvox), 'names', names, 'thresholds', thresholds, robust_opt, boot_opt, multilev_opt, 'covs', mediation_covariates, varargin{:});
                %med_results = mediation_M_search(X, Y, M, 'names', names, 'thresholds', thresholds, robust_opt, boot_opt, multilev_opt);

            case 'Search for indirect influences'
                slice_results = mediation_search('X', X(:,whvox), Y, M, varargin{:});
                
                %slice_results = mediation_search('X', X(:,whvox), Y, M, 'names', names, 'thresholds', thresholds, robust_opt, boot_opt, multilev_opt, 'covs', mediation_covariates, varargin{:});

                %med_results = mediation_X_search(X, Y, M, 'names', names, 'thresholds', thresholds, robust_opt, boot_opt, multilev_opt);

            case 'Search for mediated outcomes'
                slice_results = mediation_search('Y', X, Y(:,whvox), M, varargin{:}); % Wani added this line 01/24/15
                %slice_results = mediation_search('Y', X, Y(:,whvox), M, 'names', names, 'thresholds', thresholds, robust_opt, boot_opt, multilev_opt, 'covs', mediation_covariates, varargin{:});
                %med_results = mediation_Y_search(X, Y, M, 'names', names, 'thresholds', thresholds, robust_opt, boot_opt, multilev_opt);

            otherwise
                error('Unknown cmdstring: "%s".', cmdstring);
        end

        % write slices of output images
        write_output(slice_results, names, maskInfo, med_flags, i);

        % add results onto total results
        % (first get indices in masked output list for voxels of this slice)
        vox_this_slice = sum(z == i);  % in-mask voxels in this slice
        st = sum(z < i) + 1;      % first index in list for this slice
        en = st + vox_this_slice - 1;  % last index in list for this slice

        med_results.paths(st:en,:) = slice_results.paths;
        med_results.ste(st:en,:) = slice_results.ste;
        med_results.pvals(st:en,:) = slice_results.pvals;

        fprintf('Done in %3.0f s', etime(clock, t1));
    end
end     % End main function





% --------------------------------------
%
%
% * Sub-functions
%
%
% --------------------------------------

function [cmdstring, maskname] = get_cmdstring(X, Y, M)
    if sum([ischar(X) ischar(Y) ischar(M)]) > 1
        error('Only one of X, Y, and M can be an image set.');
    end

    if ischar(M)
        cmdstring = 'Search for mediators';
        maskname = deblank(M(1,:));
    elseif ischar(X)
        cmdstring = 'Search for indirect influences';
        maskname = deblank(X(1,:));
    elseif ischar(Y)
        cmdstring = 'Search for mediated outcomes';
        maskname = deblank(Y(1,:));
    else
        error('None of X, Y, or M is a list of filenames');
    end
end

% --------------------------------------
% * Get flags for which images to save
% --------------------------------------
function med_flags = get_med_flags(cmdstring)
    switch cmdstring
        case 'Search for mediators'
            med_flags = struct('XtoM', 1, 'MtoY', 1, 'unmediatedXtoY', 0, 'XtoY', 0, 'mediatedXtoY', 1);

        case 'Search for indirect influences'
            med_flags = struct('XtoM', 1, 'MtoY', 0, 'unmediatedXtoY', 1, 'XtoY', 1, 'mediatedXtoY', 1);

        case 'Search for mediated outcomes'
            med_flags = struct('XtoM', 0, 'MtoY', 1, 'unmediatedXtoY', 1, 'XtoY', 1, 'mediatedXtoY', 1);

        otherwise
            error('Unknown cmdstring: "%s".', cmdstring);
    end
end




function str = display_string(str)
    fprintf(str);
end



% ------------------------------------------------------------------
% Write output images/slices for whole brain or slice
% ------------------------------------------------------------------
function write_output(med_results, names, maskInfo, med_flags, varargin)
    X_to_M = 1;
    M_to_Y = 2;
    UNMED_X_to_Y = 3;
    X_to_Y = 4;
    MED_X_to_Y = 5;

    % varargin: slice number, for slice mode
    % only does anything if you enter slice number; otherwise does whole
    % brain
    slicestr = 'fullbrain'; % default
    slice_number = 0;
    
    if ~isempty(varargin)
        slicestr = 'slice';
        slice_number = varargin{1};
    end

    if med_flags.XtoM
        descrip = sprintf('X->M, %s to %s', names{1}, names{3});
        iimg_reconstruct_vols(med_results.paths(:,X_to_M), maskInfo, 'outname', 'X-M_effect.img', 'descrip', descrip, slicestr, slice_number);
        iimg_reconstruct_vols(med_results.pvals(:,X_to_M), maskInfo, 'outname', 'X-M_pvals.img', 'descrip', descrip, slicestr, slice_number);
    end

    if med_flags.MtoY
        descrip = sprintf('M->Y, %s to %s', names{3}, names{2});
        iimg_reconstruct_vols(med_results.paths(:,M_to_Y), maskInfo, 'outname', 'M-Y_effect.img', 'descrip', descrip, slicestr, slice_number);
        iimg_reconstruct_vols(med_results.pvals(:,M_to_Y), maskInfo, 'outname', 'M-Y_pvals.img', 'descrip', descrip, slicestr, slice_number);
    end

    if med_flags.unmediatedXtoY
        descrip = sprintf('X->Y direct, %s to %s', names{1}, names{2});
        iimg_reconstruct_vols(med_results.paths(:,UNMED_X_to_Y), maskInfo, 'outname', 'X-Y_direct_effect.img', 'descrip', descrip, slicestr, slice_number);
        iimg_reconstruct_vols(med_results.pvals(:,UNMED_X_to_Y), maskInfo, 'outname', 'X-Y_direct_pvals.img', 'descrip', descrip, slicestr, slice_number);
    end

    if med_flags.XtoY
        descrip = sprintf('X->Y total, %s to %s', names{1}, names{2});
        iimg_reconstruct_vols(med_results.paths(:,X_to_Y), maskInfo, 'outname', 'X-Y_total_effect.img', 'descrip', descrip, slicestr, slice_number);
        iimg_reconstruct_vols(med_results.pvals(:,X_to_Y), maskInfo, 'outname', 'X-Y_total_pvals.img', 'descrip', descrip, slicestr, slice_number);
    end

    if med_flags.mediatedXtoY
        descrip = sprintf('X->M->Y, %s thru %s to %s', names{1}, names{3}, names{2});
        iimg_reconstruct_vols(med_results.paths(:,MED_X_to_Y), maskInfo, 'outname', 'X-M-Y_effect.img', 'descrip', descrip, slicestr, slice_number);
        iimg_reconstruct_vols(med_results.pvals(:,MED_X_to_Y), maskInfo, 'outname', 'X-M-Y_pvals.img', 'descrip', descrip, slicestr, slice_number);
    end
end
