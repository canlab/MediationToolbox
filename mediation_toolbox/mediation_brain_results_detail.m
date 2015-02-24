function [paths, stats] = mediation_brain_results_detail(cl, wh_cluster, varargin)
    % [paths, stats] = mediation_brain_results_detail(cl, wh_cluster, [optional keywords])
    %
    % Tor Wager, Feb 2007
    %
    % First run mediation_brain_results, then try:
    % [paths, stats] = mediation_brain_results_detail(clpos_data, 1)
    %
    % Works on one region of interest
    % standardizes all variables so path coeffs are correlation coeffs
    %
    % Keywords (optional)
    % 'center' : center all variables
    % 'zscore' : z-score all variables
    % 'centerx' : center X.  This can be useful for plots if x has a scale
    %               not around zero
    % {'save' 'dosave'}
    %
    % Example:
    % -----------------------------------------------------------
    % %% get ab effect map
    % [clpos, clneg, clpos_data, clneg_data] = mediation_brain_results('ab',[.005 .01 .05],'overlay',overlay, 'prune');
    % clpos_data = cluster_names(clpos_data, 1);
    % save mediation_clusters2 clpos clneg clpos_data clneg_data
    %
    % %% run detail plots for each of 3 regions
    % cl_of_interest = clpos_data(1);
    % spm_orthviews('Reposition', cl_of_interest(1).mm_center);
    % [paths, stats] = mediation_brain_results_detail(cl_of_interest, 1);
    %
    %
    % cl_of_interest = clpos_data(2);
    % spm_orthviews('Reposition', cl_of_interest(1).mm_center);
    % [paths, stats] = mediation_brain_results_detail(cl_of_interest, 1);
    %
    %
    % cl_of_interest = clneg_data(1);
    % spm_orthviews('Reposition', cl_of_interest(1).mm_center);
    % [paths, stats] = mediation_brain_results_detail(cl_of_interest, 1);
    %
    % %% print raw correlations
    % names = {'DMPFCseed(X)' 'STS' 'rDMPFC' 'VMPFC' 'Neg-NeuRep'};
    % data = [SETUP.data.X clpos_data(1).timeseries clpos_data(2).timeseries clneg_data(1).timeseries];
    % y = SETUP.data.Y;
    % print_matrix(corrcoef([data y]), names, names)
    %
    % %% print partial correlations, controlling for x and y
    % x = data(:,1);
    % data = data(:,2:end);
    % [partialrs, p] = partialcorr(data,[x y])
    % print_correlation(partialrs, .4, names(2:end-1));

    docenter = 0;
    docenteronly = 0;
    centerx = 0;
    centery = 0;
    centerm = 0;
    savestr = [];
    l2mstr = [];  % level 2 moderator
    l2mdat = [];
    
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}

                % functional commands
                case 'center', docenter = 1; docenteronly = 1; centerx = 1; centery = 1; centerm = 1;
                case 'zscore', docenter = 1; docenteronly = 0; centerx = 1; centery = 1; centerm = 1;
                    
                case 'centerx', docenter = 1; docenteronly = 1; centerx = 1; 
                case 'centerm', docenter = 1; docenteronly = 1; centerm = 1; 
                case 'centery', docenter = 1; docenteronly = 1; centery = 1; 
            
                case {'save' 'dosave'}, savestr = 'dosave';
                    
                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end

    % Move orthviews
    ff = findobj(gcf, 'Tag', 'Graphics');
    if ishandle(ff)
        try
            if iscell(cl)
                spm_orthviews('Reposition', cl{wh_cluster}.mm_center);
            else
                spm_orthviews('Reposition', cl(wh_cluster).mm_center);
            end
        catch
            disp('Problem moving orthviews to cluster center.')
        end
    end
        
    [SETUP, imgs, wh_is_image, names] = load_mediation_setup;

    if iscell(cl)
        % multi-level data clusters, clpos_data format
        data = cell(1, length(cl));

        for i = 1:length(cl), data{i} = cl{i}(wh_cluster).timeseries; end
        cl_of_interest = cl{1}(wh_cluster);
        
    else
        % single-level data clusters, or clpos_data2 format with
        % timeseries values for multi-level data

        cl_of_interest = cl(wh_cluster);
        data = cat(2,cl_of_interest.timeseries);
    end

    % Level 2 moderators, for multi-level data only
    if iscell(data) && isfield(SETUP.data, 'L2M') && ~isempty(SETUP.data.L2M)
        disp('Adding 2nd-level moderators');
        l2mstr = 'L2M';
        l2mdat = SETUP.data.L2M;
    end
        
    getnames;
    names;

    if docenter
        do_centering();
    end
    

    switch wh_is_image
        case 'M'
            [paths, stats] = mediation(SETUP.data.X, SETUP.data.Y, data, 'plots', 'verbose', 'names', names, savestr, ...
                SETUP.inputOptions{:});

        case 'X'
            [paths, stats] = mediation(data, SETUP.data.Y, SETUP.data.M, 'boot', 'plots', 'verbose', 'names', names, savestr);

        case 'Y'
            [paths, stats] = mediation(SETUP.data.X, data, SETUP.data.M, 'boot', 'plots', 'verbose', 'names', names, savestr);

        otherwise
            error('Can''t determine which of X/Y/M is in cl.');

    end

    
    % ---------------------------------------
    % ---------------------------------------
    
    % Inline functions
    
    % ---------------------------------------
    % ---------------------------------------

    function getnames
        names = SETUP.names;

        if ~isfield(cl_of_interest, 'shorttitle') || isempty(cl_of_interest.shorttitle)
            cl_of_interest = cluster_names(cl_of_interest);
        end

        switch wh_is_image
            case 'M', wh_modify = 3;
            case 'X', wh_modify = 1;
            case 'Y', wh_modify = 2;

            otherwise
                error('Can''t determine which of X/Y/M is in cl.');
        end
        names{wh_modify} = sprintf('%s\n%s',names{wh_modify}, cl_of_interest.shorttitle);

    end

    function do_centering
        if iscell(SETUP.data.X)
            % Multi-level data
            N = length(SETUP.data.X);
            for i = 1:N
                if ~ischar(SETUP.data.X{1}) && centerx, SETUP.data.X{i} = scale(SETUP.data.X{i}, docenteronly); end
                if ~ischar(SETUP.data.Y{1}) && centery, SETUP.data.Y{i} = scale(SETUP.data.Y{i}, docenteronly); end
                if ~ischar(SETUP.data.M{1}) && centerm, SETUP.data.M{i} = scale(SETUP.data.M{i}, docenteronly); end
                switch wh_is_image
                    case 'M'
                        if centerm, data{i} = scale(data{i}, docenteronly); end

                    case 'X'
                        if centerx, data{i} = scale(data{i}, docenteronly); end

                    case 'Y'
                        if centery, data{i} = scale(data{i}, docenteronly); end

                    otherwise
                        error('Can''t determine which of X/Y/M is in cl.');
                end
            end
        else
            % single-level
            if ~ischar(SETUP.data.X) && centerx, SETUP.data.X = scale(SETUP.data.X, docenteronly); end
            if ~ischar(SETUP.data.Y) && centery, SETUP.data.Y = scale(SETUP.data.Y, docenteronly); end
            if ~ischar(SETUP.data.M) && centerm, SETUP.data.M = scale(SETUP.data.M, docenteronly); end
            switch wh_is_image
                    case 'M'
                        if centerm, data = scale(data, docenteronly); end

                    case 'X'
                        if centerx, data = scale(data, docenteronly); end

                    case 'Y'
                        if centery, data = scale(data, docenteronly); end

                    otherwise
                        error('Can''t determine which of X/Y/M is in cl.');
                end
            
        end
    end

    
end  % end main function


% ---------------------------------------
% ---------------------------------------

% Sub-functions

% ---------------------------------------
% ---------------------------------------

function [SETUP, imgs, wh_is_image, name] = load_mediation_setup

    SETUP = [];
    imgs = [];

    fname = [pwd filesep 'mediation_SETUP.mat'];
    if exist(fname,'file')
        load(fname);

        if ~exist('SETUP','var'), disp('No SETUP.'); return, end
        
        if isfield(SETUP, 'M') && ~isfield(SETUP, 'data')
            % single level.  rename things so it works.
            SETUP.data.M = SETUP.M;
            SETUP.data.X = SETUP.X;
            SETUP.data.Y = SETUP.Y;
        end
        
        if ~isfield(SETUP, 'inputOptions')
            disp('Warning: Missing SETUP.inputOptions.  Using default mediation options.');
            SETUP.inputOptions = {};
        end
        
        % try to find names
        if exist('SETUP','var') && (ischar(SETUP.data.M) || iscell(SETUP.data.M))
            imgs = SETUP.data.M;
            name = 'From mediation_SETUP SETUP.data.M';
            wh_is_image = 'M';

        elseif exist('SETUP','var') && (ischar(SETUP.data.X) || iscell(SETUP.data.X))
            imgs = SETUP.data.X;
            name = 'From mediation_SETUP SETUP.data.X';
            wh_is_image = 'X';

        elseif exist('SETUP','var') && (ischar(SETUP.data.Y) || iscell(SETUP.data.Y))
            imgs = SETUP.data.Y;
            name = 'From mediation_SETUP SETUP.data.Y';
            wh_is_image = 'Y';

        end

        if isempty(imgs)
            fprintf(1,'Could not find image list.\n');
        end

    else
        fprintf(1,'Go to valid mediation directory to extract image data.\n');
    end

end