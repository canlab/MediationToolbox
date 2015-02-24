function mediation_multilevel_loess_plots(stats)


    X = stats.inputOptions.X;
    M = stats.inputOptions.M;
    Y = stats.inputOptions.Y;

    create_figure('Mediation Loess Plots', 1, 4);

    subplot(1, 4, 1); hold on;

    % X to M

    stats_plot = scn_stats_helper_functions('xycatplot', X, M, 'weighted', 'samefig');

    title('X to mediator');
    drawnow

    % M to Y

    % * note: controlling for X as categorical var not done yet...

    subplot(1, 4, 2); hold on;

    scn_stats_helper_functions('loess_partial', M, Y, X, 'samefig');

    title('Mediator to outcome, direct');
    drawnow

    % X to Y partial

    % * note cat. X not done yet...

    subplot(1, 4, 3); hold on;

    scn_stats_helper_functions('loess_partial', X, Y, M, 'samefig');

    title('X to outcome, direct');
    drawnow


    % X to Y total
    subplot(1, 4, 4); hold on;
    stats_plot = scn_stats_helper_functions('xycatplot', X, Y, 'weighted', 'samefig');

    title('X to outcome, total');
    drawnow

end