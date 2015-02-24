function [clpos_data, clneg_data] = mediation_brain_print_tables(clpos_data, clneg_data, varargin)
%
% [clpos_data, clneg_data] = mediation_brain_print_tables(clpos_data, clneg_data, [do ranks flag],[other optional inputs])
%
% optional: doranks, ranks, rank : rank data before calculating partial
%                                   correlations
%           nosubpeaks :            suppress sub-peaks in output
%           a          :            report stats for peak voxel in a effect
%           b          :            report stats for peak voxel in b effect
%           ab         :            report stats for peak voxel in ab effect
% tor wager, july 2007
%
% NOTES: 
% This function returns statistics for a, b, and a*b effects for the peak
% voxel in each region.
% You must define which effect determines the peak voxel: a, b, or ab
% Either enter as input argument, or do not enter and you will prompted.
%
% Examples:
% mediation_brain_print_tables(clpos_data2, clneg_data2, 'nosubpeaks');
%

    if nargin < 3, doranks = 0; else doranks = varargin{1}; end
    
    dosubpeaks = 1;
    
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % reserved keywords
                case 'nosubpeaks', dosubpeaks = 0;
                    
                case {'ranks', 'rank', 'dorank'}, doranks = 1;

                case 'a'
                    wh_effect = 'a';
                    
                case 'b'
                    wh_effect = 'b';
                    
                case 'ab'
                    wh_effect = 'ab';
                    
                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end
    
    if ~exist('wh_effect', 'var')
        disp('This function returns statistics for a, b, and a*b effects for the peak voxel in each region.')
        disp('You must define which effect determines the peak voxel: a, b, or ab')
        wh_effect = input('Enter your choice: ', 's');
    end
    
    [SETUP, imgs, wh_is_image, name] = load_mediation_setup;

    if iscell(clpos_data) || iscell(clneg_data)
        disp('You seem to have a cell array of clusters.')
        disp('This is the format for a multi-level analysis.')
        disp('You can use this table program, but you must pass in')
        disp('clpos_data2 and clneg_data2 from mediation_brain_results.')
        disp('or use mediation_multilev_reformat_cl.m to reformat your clpos_data/neg_data.')
        
    end

    if isfield(clpos_data, 'thresholds')
        fprintf(1,'Thresholds: ');
        fprintf(1,'%3.4f ', clpos_data(1).thresholds);
        fprintf(1,'\n');

    end

    [clpos_data, clneg_data] = print_tables(clpos_data, clneg_data, SETUP, doranks, dosubpeaks, wh_effect);


end



function [clpos_data, clneg_data] = print_tables(clpos_data, clneg_data, SETUP, doranks, dosubpeaks, wh_effect)
    %% name contiguous clusters
    clpos_data = cluster_names(clpos_data, 1);
    clneg_data = cluster_names(clneg_data, 1);



    SETUP = load_mediation_setup;
    
    % Partial correlations with outcome
    % --------------------------------------------------------------
    if ~isempty(clpos_data)
        datp = cat(2, clpos_data(:).timeseries);

        if ~isfield(SETUP, 'Y')
            % we have multi-level analysis.  partial corr would be within
            % subjects. skip it.

            partialrs = NaN * zeros(1, length(clpos_data) + 1);
            p = NaN * zeros(1, length(clpos_data) + 1);

        else

            if doranks
                [partialrs, p] = partialcorr(rankdata([SETUP.Y datp]), rankdata(SETUP.X));
            else
                [partialrs, p] = partialcorr([SETUP.Y datp], SETUP.X);
            end
        end
        
        for i = 1:length(clpos_data)
            clpos_data(i).partialr_with_y = partialrs(1, i+1);
            clpos_data(i).partialr_p = p(1, i+1);
        end



        % Get max statistics for each effect
        % --------------------------------------------------------------
        myfields = {};
        
                        
        switch wh_effect
            case 'a', my_field = 'a_effect_p';
            case 'b', my_field = 'b_effect_p';
            case 'ab', my_field = 'ab_effect_p';
        end

        if ~isfield(clpos_data, my_field)
            disp(['I''m looking for the field ' my_field ' to use to define the peak voxel.'])
            error(['clusters structure does not have required field: ' my_field]);
        end
                        
        if isfield(clpos_data, 'a_effect_p') && isfield(clpos_data, 'a_effect_Z')
            for i = 1:length(clpos_data)

                [minp, wh] = min(clpos_data(i).(my_field));

                clpos_data(i).Za_max = clpos_data(i).a_effect_Z(wh(1));
                clpos_data(i).Pa_max = clpos_data(i).a_effect_p(wh(1));
            end
            
            myfields = [myfields {'Za_max' 'Pa_max'}];
        end

        if isfield(clpos_data, 'b_effect_p') && isfield(clpos_data, 'b_effect_Z')
            for i = 1:length(clpos_data)
                [minp, wh] = min(clpos_data(i).(my_field));

                clpos_data(i).Zb_max = clpos_data(i).b_effect_Z(wh(1));
                clpos_data(i).Pb_max = clpos_data(i).b_effect_p(wh(1));
            end
            
            myfields = [myfields {'Zb_max' 'Pb_max'}];
        end

        if isfield(clpos_data, 'ab_effect_p') && isfield(clpos_data, 'ab_effect_Z')
            for i = 1:length(clpos_data)
                [minp, wh] = min(clpos_data(i).(my_field));

                clpos_data(i).Zab_max = clpos_data(i).ab_effect_Z(wh(1));
                clpos_data(i).Pab_max = clpos_data(i).ab_effect_p(wh(1));
            end
            
            myfields = [myfields {'Zab_max' 'Pab_max'}];
        end

        % Print table
        % --------------------------------------------------------------
        cluster_table(clpos_data, dosubpeaks, 0, 'partialr_with_y', 'partialr_p', ...
            myfields{:}, 'num_sig_voxels');

    end


    if ~isempty(clneg_data)
        datn = cat(2, clneg_data(:).timeseries);

        if ~isfield(SETUP, 'Y')
            % we have multi-level analysis.  partial corr would be within
            % subjects. skip it.

            partialrs = NaN * zeros(1, length(clneg_data) + 1);
            p = NaN * zeros(1, length(clneg_data) + 1);

        else

            if doranks
                [partialrs, p] = partialcorr(rankdata([SETUP.Y datn]), rankdata(SETUP.X));
            else
                [partialrs, p] = partialcorr([SETUP.Y datn], SETUP.X);
            end
        end
        
        for i = 1:length(clneg_data)
            clneg_data(i).partialr_with_y = partialrs(1, i+1);
            clneg_data(i).partialr_p = p(1, i+1);
        end
        

        % Get max statistics for each effect
        % --------------------------------------------------------------
        myfields = {};
        
        switch wh_effect
            case 'a', my_field = 'a_effect_p';
            case 'b', my_field = 'b_effect_p';
            case 'ab', my_field = 'ab_effect_p';
        end

        if ~isfield(clneg_data, my_field)
            disp(['I''m looking for the field ' my_field ' to use to define the peak voxel.'])
            error(['clusters structure does not have required field: ' my_field]);
        end
        
        if isfield(clneg_data, 'a_effect_p') && isfield(clneg_data, 'a_effect_Z')
            for i = 1:length(clneg_data)

                [minp, wh] = min(clneg_data(i).(my_field));

                clneg_data(i).Za_max = clneg_data(i).a_effect_Z(wh(1));
                clneg_data(i).Pa_max = clneg_data(i).a_effect_p(wh(1));
            end
            
            myfields = [myfields {'Za_max' 'Pa_max'}];
        end

        if isfield(clneg_data, 'b_effect_p') && isfield(clneg_data, 'b_effect_Z')
            for i = 1:length(clneg_data)
                [minp, wh] = min(clneg_data(i).(my_field));

                clneg_data(i).Zb_max = clneg_data(i).b_effect_Z(wh(1));
                clneg_data(i).Pb_max = clneg_data(i).b_effect_p(wh(1));
            end
            
            myfields = [myfields {'Zb_max' 'Pb_max'}];
        end

        if isfield(clneg_data, 'ab_effect_p') && isfield(clneg_data, 'ab_effect_Z')
            for i = 1:length(clneg_data)
                [minp, wh] = min(clneg_data(i).(my_field));

                clneg_data(i).Zab_max = clneg_data(i).ab_effect_Z(wh(1));
                clneg_data(i).Pab_max = clneg_data(i).ab_effect_p(wh(1));
            end
            
            myfields = [myfields {'Zab_max' 'Pab_max'}];
        end

        % Print table
        % --------------------------------------------------------------
        cluster_table(clneg_data, dosubpeaks, 0, 'partialr_with_y', 'partialr_p', ...
            myfields{:}, 'num_sig_voxels');
    end
end


function [SETUP, imgs, wh_is_image, name] = load_mediation_setup()
    SETUP = [];
    imgs = [];

    wh_is_image = [];
    
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
                    wh_is_image = 'M';

                case 'Search for indirect influences'
                    imgs = SETUP.data.X;
                    wh_is_image = 'X';

                case 'Search for mediated outcomes'
                    imgs = SETUP.data.Y;
                    wh_is_image = 'Y';

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

% 
% function [SETUP, imgs, wh_is_image, name] = load_mediation_setup
% 
%     SETUP = [];
%     imgs = [];
% 
%     fname = [pwd filesep 'mediation_SETUP.mat'];
%     if exist(fname,'file')
%         load(fname);
% 
%         % try to find names
%         if exist('SETUP','var') && ischar(SETUP.M)
%             imgs = SETUP.M;
%             name = 'From mediation_SETUP SETUP.M';
%             wh_is_image = 'M';
% 
%         elseif exist('SETUP','var') && ischar(SETUP.X)
%             imgs = SETUP.X;
%             name = 'From mediation_SETUP SETUP.X';
%             wh_is_image = 'X';
% 
%         elseif exist('SETUP','var') && ischar(SETUP.Y)
%             imgs = SETUP.Y;
%             name = 'From mediation_SETUP SETUP.Y';
%             wh_is_image = 'Y';
% 
%         end
% 
%         if isempty(imgs)
%             fprintf(1,'Could not find image list.\n');
%         end
% 
%     else
%         fprintf(1,'Go to valid mediation directory to extract image data.\n');
%     end
% 
% end