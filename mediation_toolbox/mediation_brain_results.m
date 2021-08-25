function [clpos, clneg, clpos_data, clneg_data, clpos_data2, clneg_data2] = mediation_brain_results(meth, varargin)
%
% [clpos, clneg, clpos_data, clneg_data, clpos_data2, clneg_data2] = mediation_brain_results(meth, varargin)
%
% This is a results-printing utility that will get thresholded results
% from three types of specialized analysis directories:
% M3 (mediation), robust regression, and igls_brain analysis directories.
%
% Features are:
% - Flexibly specify height and extent thresholds
% - Optional FDR correction
% - Multiple thresholds: Default is 3 thresholds, which can be 'pruned'
%   to show only voxels contiguous with those that pass the most
%   stringent threshold
% - Can input mask at results stage. This mask will be automatically resliced to the image space.
% - Optional printing of tables
% - Optional saving of figures to disk
% - Orthviews display
% - Optional slice montages with 'slices'
% - Extracts data from clusters, and returns clusters for positive and
% negative effects in clusters structure
% - Clusters can be easily passed to other functions to render surface
% figures and other results:
%    a. Pass to mediation_brain_surface_figs or cluster_surf for surface rendering
%    b. Pass to cluster_orthviews for orthviews display
%    c. Pass data in cl(x).timeseries to mediation.m for detailed
%    follow-up mediation analyses, mediation plots, and
%    multiple-mediator analyses
% - Called by mediation_brain_results_all_script, which is a batch mode
% for this function.
%
% by Tor Wager, Feb 2007
%
% Optional arguments and methods
% --------------------------------------------
% - Optional arguments are entered as keyword/value pairs, e.g., the
% string specifying a keyword (e.g., 'mask') followed by the name of
% the mask as the next argument.
% - Many of the options below are single-keyword options; just enter
% the string as an optional argument.  These control which specific
% images/effects in the analysis directories mediation_brain_results
% uses.
%
% Mediation toolbox methods
% ----------------------------
% case 'all'
%   Do each eligible mediation results image in current directory
%   on separate orthviews
%   Clusters are returned for the LAST image
%
% % case {'xy', 'c', 'total'}
%         Total effect
%         pname = 'X-Y_total_pvals.img'; effname = 'X-Y_total_effect.img';
%
%     case {'xydirect', 'c''', 'direct'}
%         Direct effect
%         pname = 'X-Y_direct_pvals.img'; effname = 'X-Y_direct_effect.img';
%
%     case {'xm', 'a', 'x2m', 'xtom'}
%         A path
%         pname = 'X-M_pvals.img'; effname = 'X-M_effect.img';
%
%     case {'xmy', 'ab', 'mediation'}
%         Indirect effect (mediation)
%         pname = 'X-M-Y_pvals.img'; 	effname = 'X-M-Y_effect.img';
%
%     case {'my', 'b'} {'xmy', 'ab', 'mediation'}
%         B path
%         pname = 'M-Y_pvals.img'; 	effname = 'M-Y_effect.img';
%
%     case 'l2mod'
%         Second-level moderators
%         Show all paths for first L2 moderator.
%
%     case 'al2mod'
%         pname = 'ap_L2mod.img'; 	effname = 'a_L2mod.img';
%         Z_descrip2 = 'Level 2 moderator of a path';
%
%     case 'bl2mod'
%         pname = 'bp_L2mod.img'; 	effname = 'b_L2mod.img';
%         Z_descrip2 = 'Level 2 moderator of b path';
%
%     case 'abl2mod'
%         pname = 'abp_L2mod.img'; 	effname = 'ab_L2mod.img';
%         Z_descrip2 = 'Level 2 moderator of ab effect';
%
%
% Robust regression toolbox methods
% ----------------------------
%     case 'rob0'
%       robust results intercept
%
%     case 'rob1'
%       robust results covariate 1
%
%     case {'rob' 'robfit'}
%       robust results, intercept + all valid covariates in directory
%
% IGLS toolbox methods
% ----------------------------
%     case 'igls'
%     [pname, effname, nimgs] = draw_underlay_images_igls(overlay, doconj);
%
%     case 'igls slope'
%         pname = 'slope_p.img';
%         effname = 'slope_b.img';
%         Z_descrip2 = 'IGLS fixed effect of slope (group).';
%
%     case 'igls rfx'
%         pname = 'slope_LRTrfxvar_p.img';
%         effname = 'slope_rfxvar_b.img';
%         Z_descrip2 = 'IGLS random effect of slope (indiv diffs).';
%
%     case 'igls cov slope'
%         pname = 'cov_slope_p.img';
%         effname = 'cov_slope_b.img';
%         Z_descrip2 = 'IGLS fixed effect of covariate (2nd level, indiv diffs).';
%
%     case 'igls cov intercept'
%        pname = 'cov_int_p.img';
%        effname = 'cov_int_b.img';
%        Z_descrip2 = 'IGLS fixed effect of covariate (2nd level, indiv diffs).';
%
% OTHER OPTIONAL INPUTS
%
% case 'subclusters', dosubclusters = 1;
% case 'dosubclusters' dosubclusters = varargin{i+1};
%
% additional options
% --------------------------------------------
% case 'prune', pruneopt = 'prune';
% case 'add', addopt = 'add';
%
% case 'thresh', thresh = varargin{i+1};
% case 'size', szthresh = varargin{i+1};
%
% case 'overlay', overlay = varargin{i+1}; varargin{i+1} = [];
% case 'mergeclusters', domergeclusters = varargin{i+1};
% case 'poscolors', poscolors = varargin{i+1};
% case 'negcolors', negcolors = varargin{i+1};
%
% case {'conj', 'conjunction'}, doconj = 1; Conjunction across images
% of interest
%
% case 'tables'
% prints tables of clpos_data and clneg_data if requested
% with partial correlations between brain and y (***M searchonly!***)
%
% case 'slices'
% Show montages of axial, coronal, and sagittal slices of results
%
% case 'show centers only'
% with slices, show montage for slices with cls only (using
% cluster_orthviews_showcenters) - this is good if there are only a few
% clusters. otherwise, use slices, which shows a fixed number of slices
% (good for when there are many many clusters)
%
% case 'whitebackground'
% Orthviews will have white background
%
% case 'names'
% name clusters one-by-one before returning cluster/table output
% *This option is deprecated now, because this function attempts to autolabel
% using autolabel_regions_by_atlas region object method.*
%
% case 'save'
% save output in a text file in the current directory
%
%  case 'subclusters', dosubclusters = 1;
%  case 'dosubclusters', dosubclusters = varargin{i+1};
%  Separate contiguous clusters into sub-clusters around local maxima
%  (useful for tables)
%
%  case 'handle_offset', handle_offset = varargin{i+1};
%  Show results on different orthviews axis than main one
%  Useful for displaying multiple sets of results side by side
%
% output variables
% --------------------------------------------
% Most output is saved to disk when you run mediation_brain or mediation_brain_multilevel.
% The "clusters" structure, and it's newer object-oriented version--the region 
% class object--is a data format that is designed to store information about 
% contiguous suprathreshold (significant) regions ('blobs'). The cluster or 
% -class variable is a vector, with one element per brain region.  Within 
% each element, a struct-format data structure describes the region, including 
% its voxel locations, mapping from voxel to "world" (mm brain) space, statistics 
% associated with each voxel and the region average, and potentially data for 
% voxels or the region average that has been extracted from images and attached.
% By default (if output variables are requested and it can find valid data 
% images), mediation_brain_results returns clusters for both positive and 
% negative effects, with image data extracted for each region and averaged 
% over in-region voxels. This allows you to examine and plot the data, and 
% -run a descriptive mediation analysis and mediation plots for any given region.
%
% Clusters are defined relative to a threshold, and are produced by mediation_brain_results. When we ran the results above, it generated these cluster structure variables:
% clpos             Cell array with one cell per threshold entered. Clusters showing positive effects only
% clneg             A cell array with one cell per threshold entered in the analysis. Negative effects only
% clpos_data     Cell array with one cell per participant. Each cell has clusters structure with extracted mediation effects. Positive effects only. 
% clneg_data    Same as above, but negative effects only.
% clpos_data2   Same as clpos_data, but extracted data in the .timeseries field have been reorganized into a cell array with one cell per subject, suitable for direct input into mediation.m to re-run region-level analyses. 
% clneg_data2   Same as clpos_data2,  but negative effects only.
% 
% Within each cluster, the field .shorttitle contains autolabeled names and 
% .mm_center contains the coordinates of the region center. For clpos_data 
% and clneg_data, clpos_data2, clneg_data2 only, .timeseries contains 
% extracted data averaged across the region, and .all_data contains extracted 
% data from each voxel for a given participant.  The fields .a_effect_Z, 
% .a_effect_p, and so on give Z-values and P-values for each effect for each 
% voxel in the region (for the group, saved in the first subject's cell only). 
% The clusters with extracted data are also saved to the hard drive, since 
% the 'save' option was specified. For example, because we entered the 'save' 
% option above, mediation_brain_results has saved files like this one: 				
% load('cl_M-Y_pvals_003_k5_noprune.mat')
% whos cl*
% M-Y indicates the Path b effect,  003 indicates a P < 0.003 threshold, and 
% k5 indicates a 5-contiguous-voxel minimum extent threshold. Noprune indicates 
% that we did not ask for multi-threshold pruning but instead entered a single 
% threshold (see above). This yields the clpos, clneg, etc. variables 
% described above.   
%
% data extraction
% --------------------------------------------
% ...is done automatically IF the following conditions are met:
% 1) in valid mediation directory containing mediation_SETUP.mat
% 2) which contains a SETUP variable
% 3) with valid image file names in SETUP.X / Y / or M
%
% Examples:
% ---------------------------------------------
% [clpos, clneg, clpos_data, clneg_data] = ...
% mediation_brain_results('all', 'thresh', [.005 .01 .05], 'size', [5 10 10]);
%
% For robfit directory, covariate 1:
% [clpos, clneg, clpos_data, clneg_data] = mediation_brain_results('rob1', 'thresh', [.005 .01 .05], 'size', [3 1 1], 'prune', 'overlay', EXPT.overlay);
%
% To make tables after already having clusters:
% cluster_table(clpos_data, 1, 0, 'num_sig_voxels')
%
% Thresholding with FDR
% [clpos, clneg, clpos_data, clneg_data] = mediation_brain_results('all', 'thresh', [Inf .005 .01], 'size', [1 3 10], 'fdrthresh', .05, 'overlay', overlay, 'prune', 'conj');
%
% Display results and set up interactive viewing plots:
% mediation_results_interactive_view_init
%
% [clpos, clneg, clpos_data, clneg_data] = mediation_brain_results('all', 'thresh', [.005 .05 .05], 'size', [1 5 10], 'fdrthresh', .05, 'overlay', overlay, 'prune');
% iv = InteractiveViewer(args{:}, 'UseExistingGraphicsWindow', 1, 'LoadDataOnDemand', 1, 'IVObserver', mediationIVObserver('X', X, 'Y', Y));
%
% Then get all data and re-plot manually for selected voxel:
% vox_data = cell(1, length(X));
% for i = 1:length(X), vox_data{i} = getDataVolProp(iv, i, 'Ts'); end
% [paths, stats2, wistats] = mediation(X, Y, vox_data, 'boot',
% 'bootsamples', 10000, 'plots', 'verbose');
%
% OR re-run, changing scaling:
% load mediation_SETUP
% for i = 1:N, Xcentered{i} = scale(SETUP.data.X{i}, 1); Yzscore{i} = scale(SETUP.data.Y{i}); end
% [paths, stats, wistats] = mediation(Xcentered, Yzscore, stats.inputOptions.M, 'boot', 'plots', 'verbose', 'bootsamples', 10000);
%
% This example uses FDR correction across a set of images
% It also makes tables and saves clusters and figures of slices
% automatically.
% SETUP = mediation_brain_corrected_threshold('fdr');
% [clp, cln, clpp, clnn] = mediation_brain_results('b', 'thresh', ...
% SETUP.fdr_p_thresh, 'size', 1, 'overlay', overlay, 'mask', mask, ...
% 'slices', 'tables', 'names', 'save');
%
% Enter a mask at the results stage:
% This mask will be automatically resliced to the image space.
% [clpos, clneg, clpos_data, clneg_data, clpos_data2, clneg_data2] = mediation_brain_results('xmy', 'tables', 'mask', which('brainmask.nii'));

% Programmers' notes:
% Update: July 07, works with robfit directories too!
% 7/2014: Tor:  % Enter a mask at the results stage:  This mask will be automatically resliced to the image space.
%               % bug fix in iimg_multi_threshold

% -----------------------------------
% defaults
% -----------------------------------

switch spm('Ver')
    case 'SPM2'
        % spm_defaults is a script
        disp('WARNING: spm defaults not set for spm2. Make sure your defaults are set correctly');
        
    case {'SPM5' 'SPM8' 'SPM12'}
        % spm_defaults is a function
        spm_defaults()
        
    otherwise
        % unknown SPM
        disp('Unknown version of SPM!');
        spm_defaults()
end


% Initialize all variables in nested functions intended to be available in main function
clpos = [];
clpos_extent = [];
clneg = [];
clneg_extent = [];
%     cl = [];
%     text_tag = [];

[pruneopt, addopt, overlay, meth, thresh, szthresh, mask, ...
    doconj, domergeclusters, dosubclusters, handle_offset, poscolors, negcolors, ...
    dotables, doslices, doshowcenters, whitebackground, donames, dosave, nimgs, ...
    colors, pname, effname, Z_descrip2, analysistype] ...
    = defaults_and_arguments(meth);

% The body of the main function is below
% -------------------------------------------
if dosave
    text_tag = get_text_tag(pname, thresh, szthresh, pruneopt);
    logfilename = [text_tag '_log.txt'];
    outfilename = [text_tag '_results.txt'];
    
    fprintf('Saving log file: %s\nOther key results in: %s\n', ...
        logfilename, outfilename);
    
    fprintf('Will also save figures of slice montages if ''slices'' requested.\n');
    diary(logfilename)
end

% plot the blobs or blobs from all the images
% clpos and clneg are returned for the LAST image in the set
% -------------------------------------------

for wh_handle = 1:nimgs
    doorthviews();
end

% ****Need to check if no results, and put empty orthviews up if so

% put names (back) on
add_axis_names(pname, nimgs);

% legend
makelegend([thresh thresh], colors);


% check for clpos
if ~exist('clpos', 'var')
    error('No valid images found.  Are you in the wrong directory?');
end

% Extract data
% -------------------------------------------

if nargout > 2 || dotables || dosave
    clpos_data = extract_data(clpos_extent, meth, domergeclusters, dosubclusters);
    clneg_data = extract_data(clneg_extent, meth, domergeclusters, dosubclusters);
    
    [clpos_data, clneg_data] = get_table_relevant_info(clpos_data, clneg_data, overlay, Z_descrip2, donames, thresh);
end

if nargout > 4 || dotables || dosave
    % Save multi-level cluster results in easy-to-use format
    % Do nothing for single-level clusters
    % Only works/only necessary if you have data extracted
    clpos_data2 = [];
    clneg_data2 = [];
    
    if iscell(clpos_data) % aug 2010: do even if no data, because needed for tables
        clpos_data2 = mediation_multilev_reformat_cl(clpos_data);
    elseif iscell(clpos_data)
        clpos_data2 = clpos_data{1};  % no data, but just copy clusters for convenience (table printing,etc.)
    end
    
    if iscell(clneg_data)
        clneg_data2 = mediation_multilev_reformat_cl(clneg_data);
    elseif iscell(clneg_data)
        clneg_data2 = clneg_data{1};  % no data, but just copy clusters for convenience (table printing,etc.)
    end
    
end

% Output control
% -------------------------------------------
if dosave
    diary off
    diary(outfilename)
end

disp('Summary of output images:');
disp(pname)
disp('Results clusters clpos and clneg are returned for the LAST image in this set.');

if dotables
    try
        print_tables;
    catch
        disp('Error printing tables.');
    end
else
    disp('Table printing was not requested.');
end

if dosave
    diary off
    diary(logfilename)
    
    clfilename = ['cl_' text_tag];
    
    if exist([clfilename '.mat'], 'file')
        disp(['Appending clusters in ' clfilename])
        save(clfilename, '-append', 'cl*')
    else
        disp(['Saving clusters in ' clfilename])
        save(clfilename, 'cl*')
    end
    
end

% Montages of slices
if doslices
    do_slice_view;
end

if dosave
    diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Nested functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------------------
%
% Inline:
% Set up default options and process input arguments
%
% ---------------------------------------------------

    function     [pruneopt, addopt, overlay, meth, thresh, szthresh, mask, ...
            doconj, domergeclusters, dosubclusters, handle_offset, poscolors, negcolors, ...
            dotables, doslices, doshowcenters,whitebackground, donames, dosave, nimgs, ...
            colors, pname, effname, Z_descrip2, analysistype] ...
            = defaults_and_arguments(meth)
        
        pruneopt = 'noprune';
        addopt = 'noadd';
        
        if isempty(meth), meth = 'all'; end
        
        defaultoverlay = 'SPM8_colin27T1_seg.img'; %'spm2_single_subj_T1_scalped.img';
        overlay = which(defaultoverlay);
        
        thresh = [.005 .01 .05];
        szthresh = [1 1 1];
        mask = [];
        doconj = 0;
        domergeclusters = 1;
        dosubclusters = 0;
        
        handle_offset = 0;  % shift plots by this many orthviews
        
        poscolors = { [1 1 0] [1 .5 0] [1 .3 .3] };
        negcolors = { [0 0 1] [0 .5 1] [.3 .3 1] };
        
        dotables = 0;
        doslices = 0;
        doshowcenters = 0;
        donames = 0;
        dosave = 0;
        whitebackground = 0;
        
        nimgs = 1;
        
        % -------------------------------------------------------------
        % Define necessary variables based on input arguments
        % -------------------------------------------------------------
        
        for i = 1:length(varargin)
            if ischar(varargin{i}) && ~isempty(varargin{i})
                switch varargin{i}
                    % reserved keywords
                    case 'prune', pruneopt = 'prune';
                    case 'add', addopt = 'add';
                        
                    case 'thresh', thresh = varargin{i+1};
                    case 'size', szthresh = varargin{i+1};
                        
                    case 'mask', mask = varargin{i+1}; varargin{i+1} = [];
                    case 'overlay', overlay = varargin{i+1}; varargin{i+1} = [];
                        %case 'mask', mask = varargin{i+1}; varargin{i+1} = [];
                    case 'poscolors', poscolors = varargin{i+1};
                    case 'negcolors', negcolors = varargin{i+1};
                        
                    case 'mergeclusters', domergeclusters = varargin{i+1};
                        
                    case 'subclusters', dosubclusters = 1;
                        
                    case 'dosubclusters', dosubclusters = varargin{i+1};
                        
                    case {'conj', 'conjunction'}, doconj = 1;
                        
                    case 'handle_offset', handle_offset = varargin{i+1};
                        
                    case {'table', 'tables'}, dotables = 1; disp('Will print tables.')
                        
                    case 'slices', doslices = 1; disp('Will show slices.')
                        % clear existing windows
                        figh = findobj('Tag', 'montage_axial');
                        if ishandle(figh), clf(figh); end
                        
                        figh = findobj('Tag', 'montage_coronal');
                        if ishandle(figh), clf(figh); end
                        
                        figh = findobj('Tag', 'montage_sagittal');
                        if ishandle(figh), clf(figh); end
                        
                    case 'show centers only', doshowcenters = 1; disp('Will only show slices with cls.')
                        
                    case 'whitebackground', whitebackground = 1; disp('Will print with white background')
                        
                    case 'names', donames = 1;
                        
                    case 'save', dosave = 1;
                        
                    otherwise, warning('mediation:badInput', ['Unknown input string option:' varargin{i}]);
                end
            end
        end
        
        colors = {poscolors{:} negcolors{:}};
        
        % Check overlay
        if isempty(overlay)  % Maybe the user entered an empty overlay image
            overlay = which(defaultoverlay);
        end
        if isempty(overlay)
            fprintf('The default underlay image is %s\n', defaultoverlay)
            error('You have not specified an anatomical underlay image, and I cannot find the default one.');
        end
        
        %is this a robfit directory?
        [isrob,b,c]=regexp(meth, '^rob\d+$');
        if isrob
            robnum = str2num(meth(b:end)) + 1; %Yoni: not sure why +1, following previous code
            meth = 'robx';
        end
        
        switch lower(meth)
            case {'xy', 'c', 'total'}
                pname = 'X-Y_total_pvals.img'; effname = 'X-Y_total_effect.img';
                Z_descrip2 = 'Mediation c (total) effect';
                analysistype = 'mediation';
                
            case {'xydirect', 'c''', 'direct'}
                pname = 'X-Y_direct_pvals.img'; effname = 'X-Y_direct_effect.img';
                Z_descrip2 = 'Mediation direct (c'') effect';
                analysistype = 'mediation';
                
            case {'xm', 'a', 'x2m', 'xtom'}
                pname = 'X-M_pvals.img'; effname = 'X-M_effect.img';
                Z_descrip2 = 'Mediation a effect';
                analysistype = 'mediation';
                
            case {'xmy', 'ab', 'mediation'}
                pname = 'X-M-Y_pvals.img'; 	effname = 'X-M-Y_effect.img';
                Z_descrip2 = 'Mediation ab effect';
                analysistype = 'mediation';
                
            case {'my', 'b'}
                pname = 'M-Y_pvals.img'; 	effname = 'M-Y_effect.img';
                Z_descrip2 = 'Mediation b effect';
                analysistype = 'mediation';
                
            case 'all'
                analysistype = 'mediation';
                [pname, effname, nimgs] = draw_underlay_images(overlay, doconj);
                addopt = 'add';
                if doconj
                    Z_descrip2 = 'Mediation conjunction (max p-value for a, b, ab)';
                else
                    Z_descrip2 = 'Mediation p-values (ab)';
                end
                
            case 'l2mod'
                analysistype = 'mediation';
                [pname, effname, nimgs] = draw_underlay_images(overlay, doconj, 'l2mod');
                addopt = 'add';
                if doconj
                    Z_descrip2 = 'Mediation conjunction (max p-value for a, b, ab)';
                else
                    Z_descrip2 = 'Mediation p-values (ab)';
                end
                
            case 'al2mod'
                analysistype = 'mediation';
                pname = 'ap_L2mod.img'; 	effname = 'a_L2mod.img';
                Z_descrip2 = 'Level 2 moderator of a path';
                
            case 'bl2mod'
                analysistype = 'mediation';
                pname = 'bp_L2mod.img'; 	effname = 'b_L2mod.img';
                Z_descrip2 = 'Level 2 moderator of b path';
                
            case 'abl2mod'
                analysistype = 'mediation';
                pname = 'abp_L2mod.img'; 	effname = 'ab_L2mod.img';
                Z_descrip2 = 'Level 2 moderator of ab effect';
                
            case 'robx'
                analysistype = 'robfit';
                pname = sprintf('rob_p_00%02i.img', robnum);
                effname = sprintf('rob_tmap_00%02i.img', robnum);
                Z_descrip2 = 'Robust regression predictor.';
                
            case {'rob' 'robfit'}
                analysistype = 'robfit';
                [pname, effname, nimgs] = draw_underlay_images_robfit(overlay, doconj);
                addopt = 'add';
                if doconj
                    Z_descrip2 = '***CONJ not implemented yet. Robust log(1/p)';
                else
                    Z_descrip2 = 'Robust log(1/p)';
                end
                
            case 'igls'
                analysistype = 'igls';
                [pname, effname, nimgs] = draw_underlay_images_igls(overlay, doconj);
                
                addopt = 'add';
                if doconj
                    Z_descrip2 = '***CONJ not implemented yet. Returning Log(1/p)';
                else
                    Z_descrip2 = 'Log(1/p)';
                end
                
                
            case 'igls slope'
                analysistype = 'igls';
                pname = 'slope_p.img';
                effname = 'slope_b.img';
                Z_descrip2 = 'IGLS fixed effect of slope (group).';
                
            case 'igls rfx'
                analysistype = 'igls';
                pname = 'slope_LRTrfxvar_p.img';
                effname = 'slope_rfxvar_b.img';
                Z_descrip2 = 'IGLS random effect of slope (indiv diffs).';
                
            case 'igls cov slope'
                analysistype = 'igls';
                pname = 'cov_slope_p.img';
                effname = 'cov_slope_b.img';
                Z_descrip2 = 'IGLS fixed effect of covariate (2nd level, indiv diffs).';
                
            case 'igls cov intercept'
                analysistype = 'igls';
                pname = 'cov_int_p.img';
                effname = 'cov_int_b.img';
                Z_descrip2 = 'IGLS fixed effect of covariate (2nd level, indiv diffs).';
                
            otherwise
                error('Unknown command option.')
        end
        
        
        % Mask: resample to image space if necessary
        if ~isempty(mask)
            mask = scn_map_image(mask, effname);
            mask = mask(:);
        end
        
    end


% ---------------------------------------------------
%
% Nested:
% Show orthviews
%
% ---------------------------------------------------

    function doorthviews
        
        % Note: if spm graphics window was closed, will not display correctly
        % on the next loading.  (this bug could be fixed someday...)
        % need to initialize graphics window if it doesn't exist
        
        % positive
        
        this_p_img = deblank(pname(wh_handle, :));
        this_eff_img = deblank(effname(wh_handle, :));
        
        [clpos,tmp, clpos_extent] = iimg_multi_threshold(this_p_img, ...
            'thresh', thresh, 'size', szthresh, ...
            'p', 'overlay', overlay, 'colors', poscolors, ...
            'mask', mask, ...
            'pos', this_eff_img, addopt, pruneopt, 'handle', wh_handle + handle_offset);
        
        
        if ~isempty(clpos_extent)
            for ii = 1:length(clpos_extent)
                clpos_extent(ii).Z = log(1./clpos_extent(ii).Z);
                clpos_extent(ii).Z_descrip = 'Log(1/p)';
            end
            
        end
        
        if ~isempty(cat(2, clpos{:})), addopt = 'add'; end
        
        % negative
        
        [clneg,tmp, clneg_extent] = iimg_multi_threshold(this_p_img, ...
            'thresh', thresh, 'size', szthresh, ...
            'p', 'overlay', overlay, 'colors', negcolors, ...
            'mask', mask, ...
            'neg', this_eff_img, addopt, pruneopt, 'handle', wh_handle + handle_offset);
        
        if ~isempty(clneg_extent)
            for ii = 1:length(clneg_extent)
                clneg_extent(ii).Z = -log(clneg_extent(ii).Z);
                clneg_extent(ii).Z_descrip = 'Log(1/p)';
            end
        end
        
        
        if isempty(cat(2, clpos{:})) && isempty(cat(2, clneg{:}))
            fprintf('\n----------------------\nNo significant results\n----------------------\n');
            if ~strcmp(addopt, 'add'), spm_check_registration(overlay); end
            
        end
        
    end



% ---------------------------------------------------
%
% Nested:
% Show and save montages of slices
%
% ---------------------------------------------------

    function do_slice_view
        
        if ~exist('clpos_data', 'var') || ~exist('clneg_data', 'var')
            disp('You must request cluster output (clpos_data and clneg_data) to show slices.');
            return
        end
        
        text_tag = get_text_tag(pname, thresh, szthresh, pruneopt);
        
        if iscell(clpos_data) || iscell(clneg_data)
            if iscell(clpos_data), cl = clpos_data{1}; end
            
            if iscell(clpos_data) && iscell(clneg_data)
                cl = [cl clneg_data{1}];
                
            elseif iscell(clneg_data)
                cl = clneg_data{1};
                
            else
                % no clusters
            end
            
        else
            cl = [clpos_data clneg_data];
        end
        
        if isempty(cl)
            disp('No clusters to show slices for.')
            return
        end
        
        % activate the correct axis
        whichorth = size(pname, 1);
        if whitebackground
            spm_orthviews_white_background
        end
        if doshowcenters
            cluster_orthviews_showcenters(cl, 'coronal', overlay, 0, 1, [1 1 1] , 'whichorth', whichorth);
        else
            cluster_orthviews_montage(10, 'coronal', overlay, 'onerow');
        end
        %
        h = findobj(gcf,'Type', 'text','Color', 'g');
        delete(h)
        
        if dosave
            disp(['Saving image of slices: orth_' text_tag '_coronal.png'])
            %saveas(gcf,['orth_' text_tag '_coronal'], 'png')
            %print(gcf, '-dpng', '-r400', '-zbuffer', ['orth_' text_tag '_coronal'])
            print(gcf, '-dpng', '-r400', ['orth_' text_tag '_coronal'])
        end
        
        if doshowcenters
            cluster_orthviews_showcenters(cl, 'sagittal', overlay, 0, 1,'whichorth', whichorth);
        else
            cluster_orthviews_montage(10, 'sagittal', overlay, 'onerow');
        end
        
        h = findobj(gcf,'Type', 'text','Color', 'g');
        delete(h)
        
        if dosave
            disp(['Saving image of slices: orth_' text_tag '_sagittal.png'])
            %saveas(gcf,['orth_' text_tag '_sagittal'], 'png')
            %print(gcf, '-dpng', '-r400', '-zbuffer', ['orth_' text_tag '_sagittal'])
            print(gcf, '-dpng', '-r400', ['orth_' text_tag '_sagittal'])
        end
        
        if doshowcenters
            cluster_orthviews_showcenters(cl, 'axial', overlay, 0,1, 'whichorth', whichorth);
        else
            
            cluster_orthviews_montage(6, 'axial', overlay, 'onerow');
        end
        
        h = findobj(gcf,'Type', 'text','Color', 'g');
        delete(h)
        
        if dosave
            disp(['Saving image of slices: orth_' text_tag '_axial.png'])
            %saveas(gcf,['orth_' text_tag '_axial'], 'png')
            %print(gcf, '-dpng', '-r400', '-zbuffer', ['orth_' text_tag '_axial'])
            print(gcf, '-dpng', '-r400', ['orth_' text_tag '_axial'])
            
        end
        
    end

% ---------------------------------------------------
%
% Inline:
% Print tables
%
% ---------------------------------------------------
    function print_tables
        disp(' ')
        disp('Printing Tables.')
        
        if iscell(clpos_data) || iscell(clneg_data)
            
            if isfield(clpos_data2, 'timeseries') && isfield(clneg_data2, 'timeseries')
                switch analysistype
                    case 'mediation'
                        
                        wh_effect = 'ab';
                        switch meth
                            case {'xm', 'a', 'x2m', 'xtom', 'al2mod'}, wh_effect = 'a';
                            case {'my', 'b', 'bl2mod'}, wh_effect = 'b';
                            case {'xmy', 'ab', 'mediation', 'abl2mod'}, wh_effect = 'ab';
                            otherwise
                        end
                        
                        [clpos_data2, clneg_data2] = mediation_brain_print_tables(clpos_data2, clneg_data2, 'nosubpeaks', wh_effect);
                    case {'robfit', 'igls'}
                        print_simple_tables(clpos_data2, clneg_data2);
                    otherwise
                        error('unknown analysis type input.');
                end
            else
                disp('Full tables can only be printed if valid data are extracted.')
                disp('I did not find valid images to extract from, so I''m printing')
                disp('abbreviated tables.');
                disp(' ')
                
                print_simple_tables(clpos_data2, clneg_data2);
            end
            
        elseif ~isempty(clpos_data) || ~isempty(clneg_data)
            disp('NOTE: FULL MEDIATION TABLES DO NOT WORK YET FOR SINGLE-LEVEL RESULTS');
            disp('Printing abbreviated tables.');
            
            print_simple_tables(clpos_data, clneg_data);
            
            %[clpos_data, clneg_data] = mediation_brain_print_tables(clpos_data, clneg_data, 'nosubpeaks');
        end
        
        % update names if we didn't already have them in clpos_data{x}
        % we always have clpos_data2 if we're doing tables.
        %if exist('clpos_data2', 'var')
        for i = 1:length(clpos_data2)
            clpos_data{1}(i).shorttitle = clpos_data2(i).shorttitle;
        end
        
        for i = 1:length(clneg_data2)
            clneg_data{1}(i).shorttitle = clneg_data2(i).shorttitle;
        end
        
        %end
        
    end

end % main function



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------------------
%
% Subfunction:
% Extract data from contiguous regions
%
% ---------------------------------------------------
function clall = extract_data(clin, meth, domergeclusters, dosubclusters)
fprintf('Extracting image data for significant clusters. ');

% note: as of 8/1/09, pass in clpos_extent and clneg_extent
clall = [];

if isempty(clin) || ( iscell(clin) && isempty(cat(2, clin{:})) )
    fprintf('No clusters to extract.\n');
    return
end

if iscell(clin)
    % only for OLD way; pass in cell clpos{iii}
    % note: as of 8/1/09, pass in clpos_extent and clneg_extent, not cells
    if(domergeclusters)
        clall = clin{1};
        for i = 2:length(clin)
            clall = merge_clusters(clall, clin{i});
        end
        
        clall = clusters2CLU(clall);
    else
        clall = clin;
    end
    
    if isempty(clall)
        disp('Clusters is empty. Doing nothing.');
        return
    end
    
    % get clusters
    clall = tor_extract_rois([], clall, clall);
else
    % new way: pass in extent region
    clall = clin;
end

if(dosubclusters)
    clall = subclusters_from_local_max(clall, 10);
end


% load SETUP
% ----------------------------------------
switch meth
    case {'rob' 'rob0', 'rob1', 'rob2' 'rob3' 'rob4'}
        % it's a robfit directory
        [SETUP, imgs] = load_robfit_setup();
        
    case {'igls' 'igls slope' 'igls rfx' 'igls cov' 'igls cov slope' 'igls cov intercept'}
        [SETUP, imgs] = load_igls_setup();
        
    otherwise
        % it's a mediation directory
        [SETUP, imgs] = load_mediation_setup();
end
if isempty(SETUP) || isempty(imgs), return, end

% extract data
% ----------------------------------------
% multilevel: imgs is a cell array of names
if iscell(imgs)
    fprintf('Extracting data for multilevel mediation\n');
    clsubjects = cell(1, length(imgs));
    
    for i=1:length(imgs)
        fprintf(' %03d', i);
        
        filename = remove_commas_and_spaces(imgs{i}(1, :));
        
        if ~exist(filename, 'file')
            disp('Image files are not valid files! Cannot extract data');
            clsubjects{i} = tor_extract_rois([], clall); %, %clall);
        else
            clsubjects{i} = tor_extract_rois(imgs{i}, clall); % , clall);
        end
    end
    
    fprintf('\n');
    clall = clsubjects;
    
elseif exist(remove_commas_and_spaces(imgs(1, :)), 'file')
    % single level
    clall = tor_extract_rois(imgs, clall);
else
    % no files; just return clusters
    clall = tor_extract_rois([], clall);
    fprintf('Images are not valid files.\n');
    return
end

fprintf('\n');
end


% ---------------------------------------------------
%
% Subfunction:
% Support for naming
%
% ---------------------------------------------------
function filename = remove_commas_and_spaces(filename)
% remove , if exists
filename = deblank(filename);
wh = find(filename == ',');
if ~isempty(wh), filename(wh:end) = []; end
end



% ---------------------------------------------------
%
% Subfunction:
% Load SETUP file with info (specific to mediation)
%
% ---------------------------------------------------
function [SETUP, imgs, wh_is_image, name] = load_mediation_setup()
SETUP = [];
imgs = [];

fname = [pwd filesep 'mediation_SETUP.mat'];
if exist(fname, 'file')
    load(fname);
    
    % try to find names (single level)
    if exist('SETUP', 'var') && isfield(SETUP, 'M') && ischar(SETUP.M)
        imgs = SETUP.M;
        name = 'From mediation_SETUP SETUP.M';
        wh_is_image = 'M';
        
    elseif exist('SETUP', 'var') && isfield(SETUP, 'X') && ischar(SETUP.X)
        imgs = SETUP.X;
        name = 'From mediation_SETUP SETUP.X';
        wh_is_image = 'X';
        
    elseif exist('SETUP', 'var') && isfield(SETUP, 'Y') && ischar(SETUP.Y)
        imgs = SETUP.Y;
        name = 'From mediation_SETUP SETUP.Y';
        wh_is_image = 'Y';
        
    end
    
    % try to find names: multi-level
    if isfield(SETUP, 'data')
        switch SETUP.cmdstring
            case 'Search for mediators'
                imgs = SETUP.data.M;
                
            case 'Search for indirect influences'
                imgs = SETUP.data.X;
                
            case 'Search for mediated outcomes'
                imgs = SETUP.data.Y;
                
            otherwise
                error('Unknown cmdstring: "%s".', cmdstring);
        end
        
        if ~iscell(imgs) || ~ischar(imgs{1})
            imgs = []; % invalid data here
        end
        
        name = 'Multilevel, from SETUP.data.';
        
    end
    
    if isempty(imgs)
        fprintf('Could not find image list.\n');
    end
    
else
    fprintf('Go to valid mediation directory to extract image data.\n');
end
end

% ---------------------------------------------------
%
% Subfunction:
% Load SETUP file with info (specific to robust reg)
%
% ---------------------------------------------------
function [SETUP, imgs, name] = load_robfit_setup()
SETUP = [];
imgs = [];

fname = [pwd filesep 'SETUP.mat'];
if exist(fname, 'file')
    load(fname);
    
    % try to find names
    if exist('SETUP', 'var') && ischar(SETUP.files)
        imgs = SETUP.files;
        name = 'From SETUP.mat SETUP.files';
    end
    
    if isempty(imgs)
        fprintf('Could not find image list.\n');
    end
    
else
    fprintf('Go to valid robfit directory with SETUP.mat to extract image data.\n');
end
end

% ---------------------------------------------------
%
% Subfunction:
% Load SETUP file with info (specific to IGLS)
%
% ---------------------------------------------------
function [SETUP, imgs, name] = load_igls_setup()
SETUP = [];
imgs = [];

fname = [pwd filesep 'igls_SETUP.mat'];
if exist(fname, 'file')
    load(fname);
    
    % try to find names
    if exist('SETUP', 'var') && isfield(SETUP, 'data')
        imgs = SETUP.data.Y;
        name = 'From igls_SETUP.mat SETUP.data.Y';
    end
    
    if isempty(imgs)
        fprintf('Could not find image list.\n');
    end
    
    if ~exist(deblank(imgs{1}(1, :)), 'file')
        disp('Images do not appear to exist or cannot be found!')
        disp(['Sample image name: ' deblank(imgs{1}(1, :))])
        
        imgs = [];
    end
    
else
    fprintf('Go to valid IGLS results directory with igls_SETUP.mat to extract image data.\n');
end
end


% ---------------------------------------------------
%
% Subfunction:
% Draw set of underlay anatomical images in SPM orth
%
% ---------------------------------------------------

function [pname, effname, nimgs] = draw_underlay_images(overlay, doconj, varargin)
pname = []; effname = [];

dol2mod = 0;
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            case 'l2mod', dol2mod = 1;
                
            otherwise, warning('mediation:badInput', ['Unknown input string option:' varargin{i}]);
        end
    end
end

if dol2mod
    img = fullfile(pwd, 'ap_L2mod.img');
else
    img = fullfile(pwd, 'X-M_pvals.img');
end
if exist(img, 'file'), pname = char(pname, img); end

if dol2mod
    img = fullfile(pwd, 'a_L2mod.img');
else
    img = fullfile(pwd, 'X-M_effect.img');
end
if exist(img, 'file'), effname = char(effname, img); end

if dol2mod
    img = fullfile(pwd, 'bp_L2mod.img');
else
    img = fullfile(pwd, 'M-Y_pvals.img');
end
if exist(img, 'file'), pname = char(pname, img); end

if dol2mod
    img = fullfile(pwd, 'b_L2mod.img');
else
    img = fullfile(pwd, 'M-Y_effect.img');
end
if exist(img, 'file'), effname = char(effname, img); end

if dol2mod
    img = fullfile(pwd, 'c1p_L2mod.img');
else
    img = fullfile(pwd, 'X-Y_direct_pvals.img');
end
if exist(img, 'file'), pname = char(pname, img); end

if dol2mod
    img = fullfile(pwd, 'c1_L2mod.img');
else
    img = fullfile(pwd, 'X-Y_direct_effect.img');
end
if exist(img, 'file'), effname = char(effname, img); end

if dol2mod
    img = fullfile(pwd, 'abp_L2mod.img');
else
    img = fullfile(pwd, 'X-M-Y_pvals.img');
end
if exist(img, 'file'), pname = char(pname, img); end

if dol2mod
    img = fullfile(pwd, 'ab_L2mod.img');
else
    img = fullfile(pwd, 'X-M-Y_effect.img');
end
if exist(img, 'file'), effname = char(effname, img); end

% Additional mediators!

if ~dol2mod
    img1 = fullfile(pwd, 'X-M_pvals.img');
    img1 = expand_4d_filenames(img1);
    
    img2 = fullfile(pwd, 'M-Y_pvals.img');
    img2 = expand_4d_filenames(img2);
    
    img3 = fullfile(pwd, 'X-M-Y_pvals.img');
    img3 = expand_4d_filenames(img3);
    
    eimg1 = fullfile(pwd, 'X-M_effect.img');
    eimg1 = expand_4d_filenames(eimg1);
    
    eimg2 = fullfile(pwd, 'M-Y_effect.img');
    eimg2 = expand_4d_filenames(eimg2);
    
    eimg3 = fullfile(pwd, 'X-M-Y_effect.img');
    eimg3 = expand_4d_filenames(eimg3);
    
    if size(img1, 1) > 1 && size(img2, 1) > 1
        for i = 2:size(img1, 1)
            pname = char(pname, img1(i, :));  % a path
            pname = char(pname, img2(i, :)); % b path
            pname = char(pname, img3(i, :)); % ab path
            
            effname = char(effname, eimg1(i, :));  % a path
            effname = char(effname, eimg2(i, :)); % b path
            effname = char(effname, eimg3(i, :)); % ab path
        end
    end
    
end

% Add conjunction as last image
if doconj
    
    [SETUP] = load_mediation_setup();
    
    switch SETUP.cmdstring
        case 'Search for mediators' % Search for M
            p1 = 'X-M_pvals.img';
            p2 = 'M-Y_pvals.img';
            p3 = 'X-M-Y_pvals.img';
            
        case 'Search for indirect influences' % Search for X
            p1 = 'X-M_pvals.img';
            p2 = 'X-Y_total_pvals.img';
            p3 = 'X-M-Y_pvals.img';
            
        otherwise
            p1 = []; p2 = []; p3 = [];
            
    end
    
    if exist(p1, 'file') && exist(p2, 'file') && exist(p3, 'file')
        
        cnames = char(p1, p2, p3);
        img = fullfile(pwd, 'XMY_conjunction.img');
        if exist(img, 'file'), disp('Warning! Overwriting XMY_conjunction.img'); end
        
        createXMY(cnames);
        
        % add conjunction names to list
        pname = char(pname, img);
        effname = char(effname, 'X-M-Y_effect.img');
        
    else
        disp('Warning! Cannot find all required images, and thus cannot make conjunction image.')
    end
    
end

if size(pname, 1) > 1 && size(effname, 1) > 1
    pname = pname(2:end, :);
    effname = effname(2:end, :);
end

nimgs = size(pname, 1);
if size(effname, 1) ~= nimgs, error('p-value images and effect images in directory do not match!'), end

spm_check_registration(repmat(overlay, nimgs, 1));

% put names on
add_axis_names(pname, nimgs)
end


% ---------------------------------------------------
%
% Subfunction:
% Draw underlay images for robust reg
%
% ---------------------------------------------------

function [pname, effname, nimgs] = draw_underlay_images_robfit(overlay, doconj)
pname = []; effname = [];

dd = dir('rob_p_0*img');
if isempty(dd), error('No rob_p*img images.  You must be in robfit directory.'); end
pname = strvcat(dd(:).name);

dd = dir('rob_beta_0*img');
effname = strvcat(dd(:).name);

disp('Images:');
disp(pname);
disp(effname);


% Add conjunction as last image
if doconj
    % Not done yet, but could easily add
    % See mediation
end

nimgs = size(pname, 1);
if size(effname, 1) ~= nimgs, error('p-value images and effect images in directory do not match!'), end

spm_check_registration(repmat(overlay, nimgs, 1));

% put names on
add_axis_names(pname, nimgs)
end

% ---------------------------------------------------
%
% Subfunction:
% Draw underlay anatomical images for IGLS
%
% ---------------------------------------------------
function [pname, effname, nimgs] = draw_underlay_images_igls(overlay, doconj)

pname = []; effname = [];

img = fullfile(pwd, 'slope_p.img');
if exist(img, 'file'), pname = char(pname, img); end

img = fullfile(pwd, 'slope_LRTrfxvar_p.img');
if exist(img, 'file'), pname = char(pname, img); end

img = fullfile(pwd, 'cov_p.img');
if exist(img, 'file'), pname = char(pname, img); end

img = fullfile(pwd, 'slope_b.img');
if exist(img, 'file'), effname = char(effname, img); end

img = fullfile(pwd, 'slope_LRTrfxvar_p.img'); % no effect image now
if exist(img, 'file'), effname = char(effname, img); end

img = fullfile(pwd, 'cov_b.img');
if exist(img, 'file'), effname = char(effname, img); end

if size(pname, 1) > 1 && size(effname, 1) > 1
    pname = pname(2:end, :);
    effname = effname(2:end, :);
end

nimgs = size(pname, 1);
if size(effname, 1) ~= nimgs, error('p-value images and effect images in directory do not match!'), end

spm_check_registration(repmat(overlay, nimgs, 1));

% put names on
add_axis_names(pname, nimgs)

end


% ---------------------------------------------------
%
% Subfunction:
% Support for creating conjunction image
%
% ---------------------------------------------------
function createXMY(cnames)
disp('Creating image XMY_conjunction.img with max p-values for relevant paths');

[volInfo, dat] = iimg_read_img(cnames);
idat = max(dat, [], 2);
zeros = any(dat == 0, 2);
idat(zeros) = 0;

iimg_reconstruct_vols(idat, volInfo, 'outname', 'XMY_conjunction.img');
end



% ---------------------------------------------------
%
% Subfunction:
% Support for counting numbers of significant voxels for tables
%
% ---------------------------------------------------
function cl = count_sig_voxels(cl, thresh)
if ~isempty(cl) && isfield(cl, 'Z_descrip') && strcmp(cl(1).Z_descrip, 'Log(1/p)')
    for i = 1:length(cl)
        pvals = 1 ./ exp(cl(i).Z);
        
        cl(i).thresholds = thresh;
        for j = 1:length(thresh)
            cl(i).num_sig_voxels(j) = sum(pvals < thresh(j));
        end
    end
else
    % cannot count voxels; information in unknown format or empty clusters; skip.
end
end

% ---------------------------------------------------
%
% Subfunction:
% Get info relevant for tables
%
% ---------------------------------------------------
function [clpos_data, clneg_data] = get_table_relevant_info(clpos_data, clneg_data, overlay, Z_descrip2, donames, thresh)
if ~isempty(clpos_data)
    % count significant voxels at each threshold
    % add names
    % Tor: Edit 8/1/09; save Z field descrip as Zdescrip instead of Zdescrip2
    if(iscell(clpos_data))
        
        % Tor - Aug 2021 - make region object and autolabel
        % Some fields in clusters struct are  illegal in region
        % object, so preserve older-format cluster structure
        r = cluster2region(clpos_data{1});
        try
            warning off % we will get them otherwise...
            r = autolabel_regions_using_atlas(r);
            warning on
            for i = 1:length(clpos_data{1}), clpos_data{1}(i).shorttitle = r(i).shorttitle; end
        catch
            disp('Error autolabeling regions. Check that Neuroimaging_Pattern_Masks repository is installed on Matlab path.')
        end
        
        %             if donames
        %                 clpos_data{1} = cluster_names(clpos_data{1}, 1, overlay);
        %             else
        %                 for i = 1:length(clpos_data{1}), clpos_data{1}(i).shorttitle = ['R' num2str(i)]; end
        %             end
        
        clpos_data{1}(1).Z_descrip = Z_descrip2;
        clpos_data{1} = count_sig_voxels(clpos_data{1}, thresh);
        
        
        clpos_data{1} = get_effect_zvals(clpos_data{1});
    else
        %Single-level
        
        % Tor - Aug 2021 - make region object and autolabel
        % Some fields in clusters struct are  illegal in region
        % object, so preserve older-format cluster structure
        r = cluster2region(clpos_data);
        try
            warning off % we will get them otherwise...
            r = autolabel_regions_using_atlas(r);
            warning on
            for i = 1:length(clpos_data), clpos_data(i).shorttitle = r(i).shorttitle; end
        catch
            disp('Error autolabeling regions. Check that Neuroimaging_Pattern_Masks repository is installed on Matlab path.')
        end
        
        %             if donames
        %                 clpos_data = cluster_names(clpos_data, 1, overlay);
        %             else
        %                 for i = 1:length(clpos_data), clpos_data(i).shorttitle = ['R' num2str(i)]; end
        %             end
        
        clpos_data(1).Z_descrip = Z_descrip2;
        clpos_data = count_sig_voxels(clpos_data, thresh);
        
        
        clpos_data = get_effect_zvals(clpos_data);
    end
end

if ~isempty(clneg_data)
    % count significant voxels at each threshold
    % add cluster names
    % Tor: Edit 8/1/09; save Z field descrip as Zdescrip instead of Zdescrip2
    
    if(iscell(clneg_data))
        
        % Tor - Aug 2021 - make region object and autolabel
        % Some fields in clusters struct are  illegal in region
        % object, so preserve older-format cluster structure
        r = cluster2region(clneg_data{1});
        try
            warning off % we will get them otherwise...
            r = autolabel_regions_using_atlas(r);
            warning on
            for i = 1:length(clneg_data{1}), clneg_data{1}(i).shorttitle = r(i).shorttitle; end
        catch
            disp('Error autolabeling regions. Check that Neuroimaging_Pattern_Masks repository is installed on Matlab path.')
        end
        %
        %             if donames
        %                 clneg_data{1} = cluster_names(clneg_data{1}, 1, overlay);
        %             else
        %                 for i = 1:length(clneg_data{1}), clneg_data{1}(i).shorttitle = ['R' num2str(i)]; end
        %             end
        
        clneg_data{1}(1).Z_descrip = Z_descrip2;
        clneg_data{1} = count_sig_voxels(clneg_data{1}, thresh);
        
        
        clneg_data{1} = get_effect_zvals(clneg_data{1});
    else
        % Single-level
        
        r = cluster2region(clneg_data);
        try
            warning off % we will get them otherwise...
            r = autolabel_regions_using_atlas(r);
            warning on
            for i = 1:length(clneg_data), clneg_data(i).shorttitle = r(i).shorttitle; end
        catch
            disp('Error autolabeling regions. Check that Neuroimaging_Pattern_Masks repository is installed on Matlab path.')
        end
        
        %             if donames
        %                 clneg_data = cluster_names(clneg_data, 1, overlay);
        %             else
        %                 for i = 1:length(clneg_data), clneg_data(i).shorttitle = ['R' num2str(i)]; end
        %             end
        
        clneg_data(1).Z_descrip = Z_descrip2;
        clneg_data = count_sig_voxels(clneg_data, thresh);
        
        
        clneg_data = get_effect_zvals(clneg_data);
    end
end
end


% ---------------------------------------------------
%
% Subfunction:
% Support for getting z-scores for tables
%
% ---------------------------------------------------
function cl = get_effect_zvals(cl)
img = fullfile(pwd, 'X-M_pvals.img');
eimg = fullfile(pwd, 'X-M_effect.img');

if exist(img, 'file') && exist(eimg, 'file')
    V = spm_vol(img);
    eV = spm_vol(eimg);
    
    for i = 1:length(cl)
        XYZ = cl(i).XYZ;
        [cl(i).a_effect_Z, cl(i).a_effect_p] = get_z_p_data(XYZ, V, eV);
    end
    
end

img = fullfile(pwd, 'M-Y_pvals.img');
eimg = fullfile(pwd, 'M-Y_effect.img');

if exist(img, 'file') && exist(eimg, 'file')
    V = spm_vol(img);
    eV = spm_vol(eimg);
    
    for i = 1:length(cl)
        XYZ = cl(i).XYZ;
        [cl(i).b_effect_Z, cl(i).b_effect_p] = get_z_p_data(XYZ, V, eV);
    end
    
end

img = fullfile(pwd, 'X-M-Y_pvals.img');
eimg = fullfile(pwd, 'X-M-Y_effect.img');

if exist(img, 'file') && exist(eimg, 'file')
    V = spm_vol(img);
    eV = spm_vol(eimg);
    
    for i = 1:length(cl)
        XYZ = cl(i).XYZ;
        [cl(i).ab_effect_Z, cl(i).ab_effect_p] = get_z_p_data(XYZ, V, eV);
    end
end
end

function [Z, pvals] = get_z_p_data(XYZ, V, eV)
XYZ(4, :) = 1;

pvals = spm_get_data(V, XYZ);
Z = norminv(1 - pvals ./ 2);  % because p-values are two-tailed
dirxn = sign( spm_get_data(eV, XYZ) );
Z = Z .* dirxn;
end

% ---------------------------------------------------
%
% Subfunction:
% Support for naming orthviews axes
%
% ---------------------------------------------------
function add_axis_names(pname, nimgs)
global st
for i = 1:nimgs
    axes(st.vols{i}.ax{3}.ax)
    [tmp, myname] = fileparts(pname(i, :));
    myname(myname == '_') = ' ';
    wh = findstr(myname, 'pvals');
    myname(wh:wh+4) = [];
    
    title(myname, 'FontSize', 24, 'Color', 'k')
end
end

% ---------------------------------------------------
%
% Subfunction:
% Support for getting descriptive save file names
%
% ---------------------------------------------------
function text_tag = get_text_tag(pname, thresh, szthresh, pruneopt)
% get descriptive name for save files

[dd, ff] = fileparts(pname(end, :));
if ~isempty(ff), text_tag = ff; else text_tag = dd; end

thrstr = '_';

for i = 1:length(thresh)
    if isinf(thresh), thrstr = [thrstr 'FDR_'];
    else
        fmts = ['%0.' num2str(ceil(log10(1/thresh(i)))) 'f'];
        mystr = sprintf(fmts, thresh(i));
        wh = find(mystr == '.');
        if ~isempty(wh) && wh(1) ~= 1, mystr(wh(1)-1:wh(1)) = []; end % remove zero_dot
        thrstr = [thrstr mystr '_'];
    end
end

thrstr = [thrstr 'k'];
for i = 1:length(szthresh)
    thrstr = [thrstr num2str(szthresh(i)) '_'];
end

text_tag = [text_tag thrstr pruneopt];

text_tag(text_tag == '.') = [];

end


function print_simple_tables(clp, cln)
disp('Positive effects')
cluster_table(clp, 0, 0);

disp('Negative effects');
cluster_table(cln, 0, 0);
end

