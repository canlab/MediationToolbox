% obs = mediationIVObserver(varargin)
%   Computes mediation parameters for a given point, displays independently of the InteractiveViewer, 
%   and then saves the mediation results to the workspace
%
% E.g.
% To search the brain for mediators between existing X and Y variables:
%
% load('mediation_SETUP')
% X = SETUP.data.X;
% Y = SETUP.data.Y;
% Mimgs = SETUP.data.M;
%
% for i=1:length(Mimgs)
%     args{i*2-1} = 'AssociatedDataVols';
%     args{i*2} = Mimgs{i};
% end
% iv = InteractiveViewer(args{:}, 'IdenticalSpace', 1, 'LoadDataOnDemand', 1, 'IVObserver', mediationIVObserver('X', X, 'Y', Y));

function obs = mediationIVObserver(varargin)
    X = [];
    Y = [];
    M = [];

    for i=1:length(varargin)
        if(ischar(varargin{i}))
            switch(varargin{i})
                case 'X'
                    X = varargin{i+1};
                case 'Y'
                    Y = varargin{i+1};
                case 'M'
                    M = varargin{i+1};
            end
        end
    end

    wh_empty = [isempty(X) isempty(Y) isempty(M)];
    if(sum(wh_empty) ~= 1)
        error('Exactly one of the three variables, X, Y, or M must be empty. These will come from the image sets.');
    end

    obs = IVObserver(@mediation_plot_);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inline functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function mediation_plot_(iv, mmPos, voxPos)
        numDataSets = get(iv, 'NumDataSets');

        brainData = cell(1, numDataSets);
        for i = 1:numDataSets
            brainData{i} = getDataVolProp(iv, i, 'Timeseries');
        end

        switch(find(wh_empty))
            case 1
                X = brainData;
            case 2
                Y = brainData;
            case 3
                M = brainData;
        end


        [paths, stats, wistats] = mediation(X, Y, M, 'boot', 'plots', 'verbose', 'bootsamples', 10000);
        title(sprintf('x, y, z = %3.0f, %3.0f, %3.0f', mmPos(1), mmPos(2), mmPos(3)), 'FontSize', 14);

        assignin('base', 'paths', paths);
        assignin('base', 'stats', stats);
        assignin('base', 'wistats', wistats);
        assignin('base', sprintf('paths_for_vox_%02d_%02d_%02d', round(mmPos(1)), round(mmPos(2)), round(mmPos(3))), paths);
        assignin('base', sprintf('stats_for_vox_%02d_%02d_%02d', round(mmPos(1)), round(mmPos(2)), round(mmPos(3))), stats);
        assignin('base', sprintf('wistats_for_vox_%02d_%02d_%02d', round(mmPos(1)), round(mmPos(2)), round(mmPos(3))), wistats);
    end
end