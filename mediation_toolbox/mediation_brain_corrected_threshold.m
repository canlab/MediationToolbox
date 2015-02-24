function SETUP = mediation_brain_corrected_threshold(corr_type, varargin)
% SETUP = mediation_brain_corrected_threshold(corr_type, ['mask', mask image defining regions], ['images', p-value images])
%
% Multiple comparisons correction tool for mediation_brain
% Calculates corrected threshold and saves in mediation_SETUP.mat, SETUP variable
% Returns SETUP with corrected threshold.
% You can then use this threshold in mediation_brain_results to get results
% at corrected p-value thresholds.
%
% input: corr_type
% Must be one of the following:
% -------------------------------
% 'fdr'
%  False Discovery Rate correction across multiple effects
% (multiple images).  The idea is that you can find the threshold that
% controls the overall FDR in a mediation analysis, including first-level
% and second-level images in a multi-level mediation.  This threshold is a
% single threshold that provides control across the component tests (e.g.,
% a, b, and ab for a search for mediators).
%
% P-VALUE IMAGES TO USE
% ------------------------
% This function automatically uses certain p-value images from a mediation
% directory, depending on whether X, M, or Y is the "search" variable.
% You can enter your own image or set of images as well by using the 
% optional inputs:
% 'images', followed by a string matrix of p-value images
%
% SEARCH AREA
% ------------------------
% The search area is defined by default as the area in mask.img in the
% mediation results directory.  This image is written automatically when
% mediation_brain analyses are run.
% You can enter your own mask by using the optional inputs:
% 'mask', followed by the mask image defining the regions
% This provides facility for ROI-based correction.
%
% Examples:
% ------------------------
% In a multilevel mediation directory, type:
% SETUP = mediation_brain_corrected_threshold('fdr');
%
% To calculate FDR for only a single effect, e.g., the Path b effect in the
% example below, enter the p-value image name you'd like to use.
% You can actually do this for ANY p-value image, not just mediation analysis images:
% SETUP = mediation_brain_corrected_threshold('fdr', 'mask', mask, 'imgs', 'M-Y_pvals.img');
%
% Tor Wager, Dec 3, 2008

% Set up optional inputs and defaults
mask = 'mask.img';
imgs = [];

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            % functional commands
            case 'mask', mask = varargin{i+1}; varargin{i + 1} = [];
            case {'images', 'imgs'}, imgs = varargin{i+1}; varargin{i + 1} = [];

            otherwise, warning('mediation:badInput', ['Unknown input string option:' varargin{i}]);
        end
    end
end

switch lower(corr_type)
    case 'fdr'
        SETUP = run_fdr();

    otherwise
        error('Unknown correction type. See help.');
end

% do not save, because we may want to run results with different masks,
% etc.
%save mediation_SETUP -append SETUP


% --------------------------
% *
% INLINE FUNCTIONS
% *
% --------------------------

    function SETUP = run_fdr

        % Load SETUP
        % Choose sensible images to get FDR correction over
        % (Ignore images that do not vary across mediation tests.)

        load mediation_SETUP SETUP

        if ~isfield(SETUP, 'cmdstring'), error('SETUP : Not a valid mediation SETUP object.'); end

        if isempty(imgs)

            switch SETUP.cmdstring
                case 'Search for mediators'
                    imgs = char({'X-M_pvals.img', 'M-Y_pvals.img', 'X-M-Y_pvals.img'});

                    if isfield(SETUP, 'data') && isfield(SETUP.data, 'L2M') && ~isempty(SETUP.data.L2M)
                        imgs = char(imgs, 'abp_L2mod.img', 'ap_L2mod.img', 'bp_L2mod.img');
                    end

                case 'Search for indirect influences'
                    imgs = char({'X-M_pvals.img', 'X-Y_direct_pvals.img', 'X-M-Y_pvals.img'});

                    if isfield(SETUP, 'data') && isfield(SETUP.data, 'L2M') && ~isempty(SETUP.data.L2M)
                        imgs = char(imgs, 'ap_L2mod.img', 'c1p_L2mod.img', 'abp_L2mod.img');
                    end

                case 'Search for mediated outcomes'
                    imgs = char({'X-Y_direct_pvals.img', 'M-Y_pvals.img', 'X-M-Y_pvals.img'});

                    if isfield(SETUP, 'data') && isfield(SETUP.data, 'L2M') && ~isempty(SETUP.data.L2M)
                        imgs = char(imgs, 'c1p_L2mod.img', 'bp_L2mod.img', 'abp_L2mod.img');
                    end

                otherwise
                    error('Unknown mediation type in cmdstring field.');
            end

        end

        disp('Calculating FDR threshold across family of tests in these images:')
        disp(imgs)

        if ~exist(mask, 'file'), error('Cannot find mask image file.'); end
        maskInfo = iimg_read_img(mask, 2);

        pvals = iimg_get_data(maskInfo, imgs);
        fdr_p_thresh = FDR(pvals(:), .05);
        
        if isempty(fdr_p_thresh), fdr_p_thresh = -Inf; end
        
        fprintf('Total p-values: %7d\n', length(pvals(:)));

        fprintf('FDR threshold is %7f\n', fdr_p_thresh);

        disp('Saving in SETUP.fdr_p_thresh');

        SETUP.fdr_p_thresh = fdr_p_thresh;
        
        
%         create_figure('FDR'); plot(sort(pvals))
%         plot_horizontal_line(.001)
%         plot( .05 * (1:length(pvals(:))) ./ length(pvals(:)), 'r')

    end


end
