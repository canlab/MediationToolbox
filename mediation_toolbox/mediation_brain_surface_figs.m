function all_surf_handles = mediation_brain_surface_figs(clpos, clneg, varargin)
    % all_surf_handles = mediation_brain_surface_figs(clpos, clneg, [optional])
    %
    % To just create surface figures:
    % mediation_brain_surface_figs({}, {});
    %
    % Optional:
    %  'mm', mm_to_surface = varargin{i+1};
    %  {'poscolors','pcolors'}, pcolors = varargin{i+1};
    %  {'negcolors','ncolors'}, ncolors = varargin{i+1};
    %                 
    % Tor Wager
    %
    % Examples:
    % Plot positive clusters only, in purples
    % mediation_brain_surface_figs(clpos, [], 'poscolors', {[.3 0 .7] [.15 0 .85] [0 0 1]});
    %
    % See: mediation_brain_results.m
    
    all_surf_handles = [];
    mm_to_surface = 5;
    
    pcolors = { [1 1 0] [1 .5 0] [1 .3 .3] };
    ncolors = { [0 0 1] [0 .5 1] [.3 .3 1] };

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % functional commands
                case 'mm', mm_to_surface = varargin{i+1};
                case {'poscolors','pcolors'}, pcolors = varargin{i+1};
                case {'negcolors','ncolors'}, ncolors = varargin{i+1};

                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end

    f1 = create_figure('Mediation Surface Figures');

    figure(f1);
    subplot(3, 2, 1);
    surfh = addbrain('hires');
    set(surfh, 'FaceColor', [.5 .5 .5], 'FaceAlpha', 1);
    view(270, 0)
    add_blobs
    surfhan2 = surfh;

    all_surf_handles = [all_surf_handles surfh];
    
    figure(f1);
    axh = subplot(3, 2, 2);
    copyobj(surfhan2, axh);
    view(90, 0);
    lightRestoreSingle; axis image; axis off; lighting gouraud; material dull
    %     subplot(3, 2, 2);
    %     surfh = addbrain('hires');
    %     set(surfh, 'FaceColor', [.5 .5 .5], 'FaceAlpha', 1);
    %     view(90, 0)
    %     add_blobs

    figure(f1);
    subplot(3, 2, 3);
    surfh = addbrain('hires right');
    set(surfh, 'FaceColor', [.5 .5 .5], 'FaceAlpha', 1);
    add_blobs()

    all_surf_handles = [all_surf_handles surfh];
    
    figure(f1);
    subplot(3, 2, 4);
    surfh = addbrain('hires left');
    set(surfh, 'FaceColor', [.5 .5 .5], 'FaceAlpha', 1);
    add_blobs()
    
    all_surf_handles = [all_surf_handles surfh];

    figure(f1);
    subplot(3, 2, 5);
    surfh = addbrain('limbic'); 
    surfh = [surfh addbrain('brainstem')]; 
    axis auto
    set(surfh, 'FaceColor', [.5 .5 .5]);
    set(surfh(end-1), 'FaceAlpha', .15);
    view(135, 10)
    add_blobs()

    all_surf_handles = [all_surf_handles surfh];
    
    % orbital
    figure(f1);
    axh = subplot(3, 2, 6);
    copyobj(surfhan2, axh);
    view(180, -90);
    lightRestoreSingle; axis image; axis off; lighting gouraud; material dull


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Inline functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function add_blobs()
        %set(surfh, 'FaceColor', [.5 .5 .5], 'FaceAlpha', 1);
        lightRestoreSingle(gca);
        lighting gouraud
        axis image;
        axis off
        material dull

        if(exist('clpos', 'var'))
            for i = length(clpos):-1:1
                if ~isempty(clpos{i})
                    cluster_surf(clpos{i}, pcolors(i), mm_to_surface, surfh);
                end
            end
        end

        if(exist('clneg', 'var'))
            for i = length(clneg):-1:1
                if ~isempty(clneg{i})
                    cluster_surf(clneg{i}, ncolors(i), mm_to_surface, surfh);
                end
            end
        end
    end
end