function mediation_brain_multilevel(X, Y, M, SETUP, varargin)
    % mediation_brain_multilevel(X, Y, M, SETUP, [mediation optional inputs])
    %
    % Multilevel mediation on a set of brain images
    %
    % Inputs
    % ------------------------------------------------
    % X                     data matrix of t timepoints x N subjects
    % Y                     outcome variable; data matrix of t timepoints x N subjects
    % M                     mediating variable; data matrix of t timepoints x N subjects
    %
    % SETUP.(fields)
    % .mask                 name of mask image
    % .preprocX             flag for whether to HP filter X data
    % .preprocY             flag for whether to HP filter Y data
    % .preprocM             flag for whether to HP filter M data
    %
    % .TR                   repetition time of volume (image) acquisition
    % .HPlength             high-pass filter length, in s
    % .scans_per_session    vector of # volumes in each run, e.g., [128 128 128 128 128]
    % .dummyscans           indices of images in each run that will be modeled
    %                       with separate dummy variables
    %
    % mediation OPTIONAL INPUTS:
    % Any of the options in mediation.m
    % i.e., 'L2M' followed by data for a level-2 moderator variable, to
    % include the 2nd-level moderation maps, or 'covs' followed by
    % covariates
    % Type "help mediation" at the matlab prompt for more details.
    %
    % To add covariates:
    % mediation_brain_multilevel( ... , 'covs', {subj1covmatrix subj2covmatrix ... subjncovmatrix})
    %
    % Other OPTIONAL INPUTS: 
    % 'nopreproc' to skip preprocessing (i.e., for trial-level inputs)
    % 'custompreproc', followed by custom preprocessing handle
    %
    % Tor Wager, Oct 2007
    % Updated Jan 2010 to add search for X and Y
    %
    % Examples:
    %
    % mediation_brain_multilevel(tempvec, rating_events, trial_mags, struct('mask', '/Volumes/SCNBeta/IMAGING_DATA/NSF/brainmask.img'), 'boot', 'nopreproc')
    % mediation_brain_multilevel(tempvec, rating_events, trial_mags, struct('mask', '/Volumes/SCNBeta/IMAGING_DATA/NSF/brainmask.img', 'startslice', 3), 'boot', 'nopreproc')
    % mediation_brain_multilevel(tempvec, rating_events, trial_mags, struct('mask', 'mask.img', 'startslice', 10), 'boot', 'nopreproc')
    %
    % Temporary example:
    % CHeck input images for NSF1 study
    % imgs = trial_mags{1}(2,:); for i = 2:length(trial_mags), imgs = char(imgs, trial_mags{i}(2,:)); end
    % spm_check_registration(imgs);
    % spm_orthviews('Window', 1:length(trial_mags), [-20 20])
%
% Another example: Take existing mediation directory and use it to run a
% new one.  
% load mediation_SETUP
% cd ..
% mkdir Multilev_mediation-try3_10k
% cd Multilev_mediation-try3_10k/
% mediation_brain_multilevel(SETUP.data.X, SETUP.data.Y, SETUP.data.M, struct('mask', spm_get(1), 'startslice', 7), 'boot', 'nopreproc','bootsamples', 10000);
%
% Analysis with custom "preprocessing" : conversion to z-scores within each
% subject:
% 
% mediation_brain_multilevel(SETUP.data.X, SETUP.data.Y, SETUP.data.M,
% struct('mask', spm_get(1), 'startslice', 7, 'preprocX', 1, 'preprocY', 1,
% 'preprocM', 1), 'boot', 'custompreproc',
% @(data)scale(data),'bootsamples', 10000);
    

    % Code below is specific to searching for mediators (not X or Y search)

    
    % Get which type of search to do based on which input is an image set
    % E.g., 'Search for mediators'; 'Search for indirect influences';  'Search for mediated outcomes';
    [SETUP.cmdstring, SETUP.wh_is_mediator] = get_cmdstring(X, Y, M);
    
    
    % ---------------------------------------------------------------------
    % Set up preprocessing
    % To skip, enter 'nopreproc' as var. arg.
    % ---------------------------------------------------------------------

    [preprochandle, SETUP] = filter_setup(SETUP, X, varargin{:});

    % SETUPional: preproc Y?  Should do mostly only if not brain
    % Should not do if using trial-level estimates

    N = length(X);
    if SETUP.preprocX, for i = 1:N, X{i} = preprochandle(X{i}); end, end
    if SETUP.preprocY,  for i = 1:length(Y), Y{i} = preprochandle(Y{i}); end, end
    
    if ~SETUP.preprocM 
        preprochandle = []; 
    
    else
        tmp = cell(1, N);
        for i = 1:N, tmp{i} = preprochandle; end
        preprochandle = tmp;
    end

    % ---------------------------------------------------------------------
    % Set up analysis
    % ---------------------------------------------------------------------
    switch SETUP.wh_is_mediator
        case 'X'
            fhandle = @(X) mediation_brain_multilev_wrapper(X, Y, M, varargin{:});
        case 'Y'
            fhandle = @(Y) mediation_brain_multilev_wrapper(X, Y, M, varargin{:});
        case 'M'
            fhandle = @(M) mediation_brain_multilev_wrapper(X, Y, M, varargin{:});
        otherwise
            error('This should never happen.  Contact tor wager and ask him what''s up.');
    end

    SETUP.names = {'X' 'Y' 'M'};
    SETUP.outputnames = {'X-M_effect.img' 'M-Y_effect.img' 'X-Y_direct_effect.img' 'X-Y_total_effect.img' 'X-M-Y_effect.img' ...
        'X-M_ste.img' 'M-Y_ste.img' 'X-Y_direct_ste.img' 'X-Y_total_ste.img' 'X-M-Y_ste.img' ...
        'X-M_pvals.img' 'M-Y_pvals.img' 'X-Y_direct_pvals.img' 'X-Y_total_pvals.img' 'X-M-Y_pvals.img' ...
        'X-M_indiv_effect.img' 'M-Y_indiv_effect.img' 'X-Y_direct_indiv_effect.img' 'X-Y_total_indiv_effect.img' 'X-M-Y_indiv_effect.img' ...
        'X-M_indiv_ste.img' 'M-Y_indiv_ste.img' 'X-Y_direct_indiv_ste.img' 'X-Y_total_indiv_ste.img' 'X-M-Y_indiv_ste.img'};

    SETUP.preprochandle = preprochandle;
    SETUP.fhandle = fhandle;

    SETUP.data.descrip = 'Data after any preprocessing specified (for X and Y)';
    SETUP.data.X = X;
    SETUP.data.Y = Y;
    SETUP.data.M = M;
    
    % Check for some key inputs
    % -------------------------------------------------
    fprintf('First-level covariates: ')
    
    wh = find(strcmp('covs', varargin));
    if ~isempty(wh) && length(varargin) >= wh(end) + 1
        fprintf(' %3.0f Found\n', size(varargin{wh(end) + 1}, 2));
        
        SETUP.data.covs = varargin{wh(end) + 1};
    else
        fprintf(' None\n')
        SETUP.data.covs = [];
    end
    
    fprintf('Second-level moderators: ')
    
    wh = find(strcmp('L2M', varargin));
    if ~isempty(wh) && length(varargin) >= wh(end) + 1
        fprintf(' %3.0f Found: Creating additional L2mod output images.\n', size(varargin{wh(end) + 1}, 2));
        
        SETUP.data.L2M = varargin{wh(end) + 1};
        
        % add names to output names, because we're writing additional
        % images
        SETUP.outputnames = [SETUP.outputnames ...
            {'a_L2mod.img' 'b_L2mod.img' 'c1_L2mod.img' 'c_L2mod.img' 'ab_L2mod.img' ...
            'ap_L2mod.img' 'bp_L2mod.img' 'c1p_L2mod.img' 'cp_L2mod.img' 'abp_L2mod.img'}];
        
        
    else
        fprintf('None\n')
        
        SETUP.data.L2M = [];
    end
    
    
    SETUP.inputOptions = varargin;
    
    save mediation_SETUP SETUP

    % ---------------------------------------------------------------------
    % Run preprocessing and analysis
    % ---------------------------------------------------------------------
    if ~isfield(SETUP, 'startslice') || isempty(SETUP.startslice), SETUP.startslice = 1; end

    switch SETUP.wh_is_mediator
        case 'X'
            image_eval_function_multisubj(X, fhandle, 'mask', SETUP.mask, 'preprochandle', preprochandle, 'outnames', SETUP.outputnames, 'start', SETUP.startslice);
        case 'Y'
            image_eval_function_multisubj(Y, fhandle, 'mask', SETUP.mask, 'preprochandle', preprochandle, 'outnames', SETUP.outputnames, 'start', SETUP.startslice);
        case 'M'
            image_eval_function_multisubj(M, fhandle, 'mask', SETUP.mask, 'preprochandle', preprochandle, 'outnames', SETUP.outputnames, 'start', SETUP.startslice);
        otherwise
            error('This should never happen.  Contact tor wager and ask him what''s up.');
    end
    
    




end    % End Main Function




    % Setup command string
    % ---------------------------------------------------------------------
    function [cmdstring, wh_is_mediator] = get_cmdstring(X, Y, M)

        if ~iscell(X) || ~iscell(Y) || ~iscell(M)
            error('X, Y, and M must be cell arrays (one cell per subject) for multilevel analysis.');
        end

        if sum([ischar(X{1}) ischar(Y{1}) ischar(M{1})]) > 1
            error('Only one of X, Y, and M can be an image set.');
        end

        if ischar(M{1})
            cmdstring = 'Search for mediators';
            wh_is_mediator = 'M';
        elseif ischar(X{1})
            cmdstring = 'Search for indirect influences';
            wh_is_mediator = 'X';
        elseif ischar(Y{1})
            cmdstring = 'Search for mediated outcomes';
            wh_is_mediator = 'Y';
        else
            error('None of X, Y, or M appears to be a list of filenames');
        end
        
        disp(cmdstring); pause(1)
    end




% Set up preprocessing
function [preprochandle, SETUP] = filter_setup(SETUP, X, varargin)

    preprochandle = [];
    %wh_elim = [];
    hpflag = 1;  % only does it if requested, though

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % reserved keywords
                case 'custompreproc'
                    preprochandle = varargin{i + 1};  % e.g., 'custompreproc', @(data) scale(data) for z=scores;
                
                    hpflag = 0;
                    SETUP.TR = NaN;
                    SETUP.HPlength = [];
                    SETUP.dummyscans = [];
                    %wh_elim = i;
                    
                case {'nopreproc'}
                    hpflag = 0;
                    SETUP.preprocX = 0; SETUP.preprocY = 0; SETUP.preprocM = 0;
                    SETUP.TR = NaN;
                    SETUP.HPlength = [];
                    SETUP.dummyscans = [];
                    %wh_elim = i;
                        
                    % We need to allow mediation SETUPions here, so eliminate this from list and do not error check here.
                    %otherwise, warning(['Unknown input string SETUPion:' varargin{i}]);
            end
        end
    end

    %varargin(wh_elim) = [];

    N = fieldnames(SETUP);
    for i = 1:length(N)
        if ~isfield(SETUP, N{i}) || isempty(SETUP.(N{i}))
            switch N{i}
                case {'TR', 'mask', 'scans_per_session', 'preprocY'}
                    error(['Enter SETUP.' N{i}]);

                case 'HPlength'
                    SETUP.(N{i}) = [];

                case 'dummyscans'
                    SETUP.(N{i}) = 1:2;

                otherwise
                    disp('Warning! Unrecognized field in SETUPions structure SETUP.');
            end
        end
    end

    SETUP.preproc_any = SETUP.preprocX || SETUP.preprocY || SETUP.preprocM;

    if SETUP.preproc_any && hpflag

        [tmp, I, S] = hpfilter(X(:,1), SETUP.TR, SETUP.HPlength, SETUP.scans_per_session, SETUP.dummyscans); % creates intercept and smoothing matrices

        preprochandle = @(Y) hpfilter(Y,  [], S, SETUP.scans_per_session, I);  % function handle with embedded fixed inputs

    end

end

