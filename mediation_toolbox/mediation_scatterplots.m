function mediation_scatterplots(stats, varargin)
% mediation_scatterplots(stats, [myfontsize])
%
% Create scatterplots for output of a mediation.m stats structure
% Uses options from stats structure to create plots that correspond to
% the analysis performed.


myfontsize = 24;
for i = 1:length(varargin)
    if ~isempty(varargin{i}) && ~ischar(varargin{i})

        myfontsize = varargin{i};

    end
end


% get variables from structure
% ---------------------------------------------------------------
dorobust = strcmp(stats.inputOptions.robust, 'Yes');

fnames = {'N' 'X' 'Y' 'M' 'additionalM' 'mediation_covariates' 'vnames'};

for i = 1:length(fnames)
    eval([fnames{i} ' = stats.inputOptions.(fnames{i});']);
end

nadditional = size(additionalM, 2);

if N > 1
    covs2 = intercept(stats.inputOptions.X_2ndlevel, 'remove');

    nadditional = size(covs2, 2) - 1; % additional rows of plots

    if isempty(covs2)
        disp('No second-level covariates; omitting 2nd level scatterplots.');
        return
    end

    % We cannot plot additional 1st-level mediators here; plot 2nd-level
    % moderators of first 1st-level mediator

    n_mediators = size(stats.abetas{1}, 2);  % a paths are [int; slope] for each mediator, each S is cell

    a = cat(2, stats.abetas{:}); a = a(1, :)';
    
    for i = 1:n_mediators
        atmp(:, i) = a(i:n_mediators:end);
    end
    a = atmp(:, 1);  % just keep the first a

    b = cat(2, stats.bbetas{:}); b = b(1, :)';
    cp = cat(2, stats.cpbetas{:}); cp = cp(1, :)';
    ab = a .* b;

end

% create figure
% ---------------------------------------------------------------
switch N > 1
    case 1
        % multilevel
        f1 = create_figure('Mediation_Scatterplots', 1 + nadditional, 4);

    case 0
        % single level
        f1 = create_figure('Mediation_Scatterplots', 1 + nadditional, 3);

    otherwise, error('This should not happen.')
end

% print info
% ---------------------------------------------------------------
ynstr = {'No' 'Yes'};
fprintf('Mediation scatterplots:\n');
fprintf('Replications: %3.0f\n', N);
fprintf('Covariates controlled for in all regressions: %3.0f predictors\n', size(mediation_covariates, 2));
fprintf('Additional mediators controlled for in outcome predictions: %3.0f predictors\n', size(nadditional, 2));

% make plots
% ---------------------------------------------------------------
switch (N > 1)
    case 1
        % multilevel

        for i = 1:nadditional + 1
            subplot(1 + nadditional, 4, 4*i - 3);
            [r,str,sig,ry,rx,h] = prplot(a, covs2, i, dorobust, {'ko'});
            set(gca,'FontSize', myfontsize);
            xlabel('a effect'); ylabel(sprintf('2nd-level cov %3.0f', i));
            title('Second-level moderation')

            % re-plot with weights
            plot_with_weights(h, ry, rx, stats.w, 1)

            subplot(1 + nadditional, 4, 4*i - 2);
            [r,str,sig,ry,rx,h] = prplot(b, covs2, i, dorobust, {'ko'});
            set(gca,'FontSize', myfontsize);
            xlabel('b effect'); ylabel(sprintf('2nd-level cov %3.0f', i));
            title('  ')

            % re-plot with weights
            plot_with_weights(h, ry, rx, stats.w, 2)

            subplot(1 + nadditional, 4, 4*i - 1);
            [r,str,sig,ry,rx,h] = prplot(cp, covs2, i, dorobust, {'ko'});
            set(gca,'FontSize', myfontsize);
            xlabel('c-prime effect'); ylabel(sprintf('2nd-level cov %3.0f', i));
            title('  ')

            % re-plot with weights
            plot_with_weights(h, ry, rx, stats.w, 3)

            subplot(1 + nadditional, 4, 4*i - 0);
            [r,str,sig,ry,rx,h] = prplot(ab, covs2, i, dorobust, {'ko'});
            set(gca,'FontSize', myfontsize);
            xlabel('a*b effect'); ylabel(sprintf('2nd-level cov %3.0f', i));
            title('  ')

            % re-plot with weights
            plot_with_weights(h, ry, rx, stats.w, 5)


        end

    case 0
        % single level

        [nanvec, X, Y, M, mediation_covariates, additionalM] = nanremove(X, Y, M, mediation_covariates, additionalM);

        subplot(1 + nadditional, 3, 1);
        prplot(M, [X mediation_covariates], 1, dorobust, {'ko'});
        set(gca,'FontSize', myfontsize);
        xlabel(vnames{1}); ylabel(vnames{3});
        title('Seed predicting mediator')

        subplot(1 + nadditional, 3, 2);
        prplot(Y, [M X additionalM mediation_covariates], 1, dorobust, {'ko'});
        set(gca,'FontSize', myfontsize);
        xlabel(vnames{3}); ylabel(vnames{2});
        title('Mediator predicting outcome')

        subplot(1 + nadditional, 3, 3);
        prplot(Y, [X M additionalM mediation_covariates], 1, dorobust, {'ko'});
        set(gca,'FontSize', myfontsize);
        xlabel(vnames{1}); ylabel(vnames{2});
        title('Direct (unmediated)')

        % additional mediators
        % ---------------------------------------------------------------

        for i = 1:nadditional
            subplot(1 + nadditional, 3, (3 * i) + 1);
            % NOTE: additionalM is OK; it IS set above, in eval statement
            prplot(additionalM(:, i), [X mediation_covariates], 1, dorobust, {'ko'});
            set(gca,'FontSize', myfontsize);
            xlabel(vnames{1}); ylabel(vnames{3 + i});
            title('Seed predicting mediator')

            subplot(1 + nadditional, 3, (3 * i) + 2);
            prplot(Y, [additionalM M X  mediation_covariates], i, dorobust, {'ko'});
            set(gca,'FontSize', myfontsize);
            xlabel(vnames{3 + i}); ylabel(vnames{2});
            title('Mediator predicting outcome')

            subplot(1 + nadditional, 3, (3 * i) + 3);
            axis off
        end

    otherwise, error('This should not happen.')
end

% enlarge text
hh = findobj(f1, 'Type', 'Text'); set(hh, 'FontSize', myfontsize);

end


function plot_with_weights(h, ry, rx, w, whcol)
delete(h)
w = w(:, whcol);
N = length(ry);
w = w .* N;     % normalize so that equal weights are all 1s

basecolor = [0 0 0];
for i = 1:N
    mycolor = w(i) .* basecolor + (1-w(i)) .* ([1 1 1] - basecolor);
    mycolor(mycolor > 1) = 1;
    mycolor(mycolor < 0) = 0;
    plot(rx(i), ry(i), 'ko', 'MarkerFaceColor', mycolor);
end

end


 
