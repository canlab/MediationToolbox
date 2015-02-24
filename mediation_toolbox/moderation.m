function [paths, stats] = moderation(X, Y, M, varargin)
%
% [paths, stats] = moderation(X, Y, M, varargin)
%
% Single or multi-level moderation
% MULTILEVEL ONLY RIGHT NOW!!!
% 
% Optional inputs:
% Same as mediation.m, all optional inputs apply
%
% Y = EXPT.hr_nocovs;
% X = clpos_data2(1).timeseries;
% M = EXPT.DESIGN.speech;
%
% Moderation
% M moderates X->Y
%
% [paths, stats] = moderation(X, Y, M, 'plots', 'verbose');

N = length(Y);

for i = 1:N, XbyM{i} = scale(X{i}, 1) .* scale(Y{i}, 1); end

[paths, stats] = mediation(X, Y, XbyM, 'M', M, varargin{:});

paths = paths(:, 2);

disp('b effect is moderator effect')

%% Compute regression dependent on level of M
X2 = stats.inputOptions.X_2ndlevel;
X2 = intercept(X2, 'remove');

for i = 1:N
    my_split{i} = mediansplit(M{i});
    
    X1{i} = X{i}(my_split{i} > 0);
    Y1{i} = Y{i}(my_split{i} > 0);
    M1{i} = M{i}(my_split{i} > 0);
    X0{i} = X{i}(my_split{i} < 0);
    Y0{i} = Y{i}(my_split{i} < 0);
    M0{i} = M{i}(my_split{i} < 0);
end


wh = strmatch('plots', varargin, 'exact');
varargin(wh) = [];

stats0 = glmfit_multilevel(Y0, X0, X2, 'names', {'Intercept' 'L2 predictor1'}, varargin{:});
stats1 = glmfit_multilevel(Y1, X1, X2, 'names', {'Intercept' 'L2 predictor1'}, varargin{:});

if ~isempty(wh)
    % plot averages
    create_figure('Y as a function of X', 1, 2);
    scn_stats_helper_functions('xyplot', X0, Y0, 'samefig', 'weighted', 'colors', {'b'}, 'nostats', 'noind');

    scn_stats_helper_functions('xyplot', X1, Y1, 'samefig', 'weighted', 'colors', {'r'}, 'nostats', 'noind');
    disp('blue = Moderator low, red = Moderator high')
end



end

