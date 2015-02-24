function cl = mediation_extract_data(cl, meth, domergeclusters)
    %cl = mediation_extract_data(cl)
    %
    % Works with robust dir too, if you enter:
    % cl = mediation_extract_data(cl, 'rob');
    %
    %help not done -- see medation_+brain_results.m
    % 
    % Start in mediation or robfit directory with a valid SETUP file.
    % This function will extract data from images in SETUP.data.M (for
    % mediation dir) or SETUP.files (for robfit directory).
    %
    % You can specify any clusters cl variable you want, e.g., one from a
    % different analysis.
    %

    if nargin < 3, domergeclusters = 1; end
    if nargin < 2, meth = 'mediation'; end
    
    if ~iscell(cl), cl = {cl}; end
    
    cl = extract_data(cl, meth, domergeclusters);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function clall = extract_data(clin, meth, domergeclusters)
    fprintf('Extracting image data for significant clusters. ');

    clall = []; 
    
    if isempty(clin) || ( iscell(clin) && isempty(cat(2, clin{:})) )
        fprintf('No clusters to extract.\n');
        return
    end

    if(domergeclusters)
        clall = clin{1};
        for i = 2:length(clin)
            clall = merge_clusters(clall, clin{i});
        end

        clall = clusters2CLU(clall);
    else
        clall = clusters2CLU(clin{1});
    end
    
    if isempty(clall)
        disp('Clusters is empty. Doing nothing.');
        return
    end

    % load SETUP
    % ----------------------------------------
    switch meth
        case {'rob' 'rob0', 'rob1', 'rob2' 'rob3' 'rob4'}
            % it's a robfit directory
            [SETUP, imgs] = load_robfit_setup();
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
            if ~exist(deblank(imgs{i}(1, :)), 'file'), error('Not valid files!'); end
            clsubjects{i} = tor_extract_rois(imgs{i}, clall, clall);
        end
        fprintf('\n');
        clall = clsubjects;
        
    elseif exist(deblank(imgs(1, :)), 'file')
        % single level
        clall = tor_extract_rois(imgs, clall, clall);
    else
        % no files; just return clusters
        clall = tor_extract_rois([], clall, clall);
        fprintf('Images are not valid files.\n');
        return
    end

    fprintf('\n');
end



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
