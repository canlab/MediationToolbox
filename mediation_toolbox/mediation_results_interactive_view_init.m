function mediation_results_interactive_view_init
    %     mediation_results_interactive_view_init
    %
    % Works only for search for mediators now!  Needs updating.

    % need to do only once

    disp('Initializing interactive viewer')
    disp('------------------------------------------')
    
    disp('Loading mediation_SETUP.mat')
    [SETUP, imgs, wh_is_image, name] = load_mediation_setup;

    disp('Mapping mask volume')
    %V = spm_vol(SETUP.mask);  % do not use this because it may have wrong
    %dims
    V = spm_vol('mask.img');

    %imgs = SETUP.data.M;  % trial magnitues, one imgs per trial, one cell per subject
    Xdat = SETUP.data.X;
    Ydat = SETUP.data.Y;

    disp('Mapping image volumes')
    N = length(imgs);
    VV = cell(1, N);
    for i = 1:N, VV{i} = spm_vol(imgs{i}); end

    disp('Setting graphics callback')

        callback_handle = @(str, pos, reg, hReg) mediation_interactive_callback_wrapper(str, pos, reg, hReg);

     hSpmFig = spm_figure('GetWin', 'Graphics');
     
     hReg = uicontrol(hSpmFig, 'Style', 'Text', 'String', 'InteractiveViewer hReg', ...
            'Position', [100 200 100 025], 'Visible', 'Off', ...
            'FontName', 'Times', 'FontSize', 14, 'FontWeight', 'Bold', ...
            'HorizontalAlignment', 'Center');
        hReg = spm_XYZreg('InitReg', hReg, V.mat, V.dim(1:3)');
        spm_XYZreg('Add2Reg', hReg, 0, callback_handle);
        spm_orthviews('Register', hReg);
    

    %fh = findobj('Tag', 'Graphics');
    %set(fh, 'WindowButtonUpFcn', callback_handle);

    disp('Ready!')

    
    % inline

    function mediation_interactive_callback_wrapper(str, pos, reg, hReg)

        switch str
            case 'SetCoords'
                mediation_interactive_callback(pos, Xdat, Ydat, VV, V);
                
            otherwise
                disp('Unknown callback command from spm_XYZreg');
        end

    end

end





function mediation_interactive_callback(pos, Xdat, Ydat, VV, V)
    %% Done each time you click

    %pos = spm_orthviews('Pos');
    vox = (mm2voxel(pos', V))';
    vox(4) = 1;

    % Load data for this voxel
    N = length(VV);
    for i = 1:N, braindata{i} = spm_get_data(VV{i}, vox); end

    % Run
    [paths, stats, wistats] = mediation(Xdat, Ydat, braindata, 'boot', 'plots', 'verbose', 'bootsamples', 10000);

    % %
    % % clear Xcentered Ycentered Mcentered
    % % for i = 1:N, Xcentered{i} = Xdat{i} - nanmean(Xdat{i}); Ycentered{i} = Ydat{i} - nanmean(Ydat{i}); Mcentered{i} = braindata{i} - nanmean(braindata{i}); end
    % % [paths, stats, wistats] = mediation(Xcentered, Ycentered, Mcentered, 'boot', 'plots', 'verbose', 'bootsamples', 10000);
    % %
    % % for i = 1:N, brainzscore{i} = braindata{i} - nanmean(braindata{i}); end
    % %
    % % for i = 1:N, Yzscore{i} = scale(Ydat{i}); end
    % % [paths, stats, wistats] = mediation(Xdat, Yzscore, braindata, 'boot', 'plots', 'verbose', 'bootsamples', 10000);

    disp('-------------------------------------------------------');
    disp('Assigning paths, stats, and wistats in base workspace.');
    disp('Data from mediation is in stats.inputOptions.X, Y, M');
    disp('-------------------------------------------------------');
    disp(' ')

    assignin('base', 'paths', paths);
    assignin('base', 'stats', stats);
    assignin('base', 'wistats', wistats);

end



function [SETUP, imgs, wh_is_image, name] = load_mediation_setup

    SETUP = [];
    imgs = [];

    fname = [pwd filesep 'mediation_SETUP.mat'];
    if exist(fname,'file')
        load(fname);

        % try to find names (single level)
        if exist('SETUP','var') && isfield(SETUP, 'M') && ischar(SETUP.M)
            imgs = SETUP.M;
            name = 'From mediation_SETUP SETUP.M';
            wh_is_image = 'M';

        elseif exist('SETUP','var') && isfield(SETUP, 'X') && ischar(SETUP.X)
            imgs = SETUP.X;
            name = 'From mediation_SETUP SETUP.X';
            wh_is_image = 'X';

        elseif exist('SETUP','var') && isfield(SETUP, 'Y') && ischar(SETUP.Y)
            imgs = SETUP.Y;
            name = 'From mediation_SETUP SETUP.Y';
            wh_is_image = 'Y';

        end

        % try to find names: multi-level
        if isfield(SETUP, 'data')
            switch SETUP.cmdstring
                case 'Search for mediators'
                    imgs = SETUP.data.M;
                    name = 'Multilevel, from SETUP.data.M';
                    wh_is_image = 'M';

                case 'Search for indirect influences'
                    imgs = SETUP.data.X;
                    name = 'Multilevel, from SETUP.data.X';
                    wh_is_image = 'X';
                    
                case 'Search for mediated outcomes'
                    imgs = SETUP.data.Y;
                    name = 'Multilevel, from SETUP.data.Y';
                    wh_is_image = 'Y';
                    
                otherwise
                    error('Unknown cmdstring: "%s".', cmdstring);
            end

            if ~iscell(imgs) || ~ischar(imgs{1})
                imgs = []; % invalid data here
            end

            

        end



        if isempty(imgs)
            fprintf(1,'Could not find image list.\n');
        end

    else
        fprintf(1,'Go to valid mediation directory with SETUP.mat to use interactive plotting.\n');
    end

end

