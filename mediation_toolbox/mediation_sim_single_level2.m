
% function mediation_sim2_igls(iter,varargin)
% Simulation for power and false positive rates for mediation analysis
%
% tor wager, Feb. 2007, Updated March 2007
% -------------------------------------------------------------------------
% Default: Bootstrap 1000 samples, AR(2), hierarchical weighting, no
% shift/latent
%
% case {'boot'}, bootopt = 'boot';
% case 'noboot', bootopt = 'no';
%
% case 'ar', arstropt = 'ar';
% case 'noar', arstropt = 'no';
%
% case 'latent', latentopt = 'latent';
% case 'nolatent', latentopt = 'no';
%
% case 'hierarchical', weightopt = 'hierarchical';
% case 'summarystats', weightopt = 'summarystats';
%
% case 'bootsamples', bootsamples = varargin{i+1};
%
%
% % Examples:
% -------------------------------------------------------------------------
%
% mediation_sim2_igls(10000,'bootsamples',5000); mediation_sim_output_figs(5);
%
% Generate data with different noise sigmas for each 'participant'
% N = 15; sigma = linspace(3, 5, N);
% mediation_sim2_igls(1000,'N',N,'sigma',sigma, 'noar','igls_type','i'); mediation_sim_output_figs(5);

function summary = mediation_sim_single_level2(niter, varargin)

%diary mediation_sim_output.txt
%%
% OLS vs. BOOTSTRAP

% default parameters for data generation
% -----------------------------------------------------------
N = 30;         % number of subjects
sigma = 2;
%niter = 100;

clear summary
summary.N = N;

for n = 1:length(N)

    % Generate data: Null hypothesis
    % ---------------------------------------------------------------
    randn('state',sum(100*clock))

    X = normrnd(0, sigma, N(n), niter);
    Y = normrnd(0, sigma, N(n), niter);
    M = normrnd(0, sigma, N(n), niter);


    % OLS
    % ---------------------------------------------------------------

    summary.null_p_ols = zeros(niter, 5);

    tt = clock;
    for i = 1:niter
        if mod(i, 50) == 0, fprintf('%3.0f ', i); end
        if mod(i, 1000) == 0, fprintf('\n'); end


        [paths, stats] = mediation(X(:, i), Y(:, i), M(:, i));

        summary.null_p_ols(i, :) = stats.p;

    end

    et = etime(clock, tt);

    summary.fpr_ols = sum(summary.null_p_ols < .05) ./ niter;
    summary.runtime_ols = et ./ niter;

    fprintf(1,'\nFalse-positive rates, p < .05\n-----------------------------------\n\n');
    fprintf(1,'OLS\t');
    print_matrix(summary.fpr_ols, {'a' 'b' 'c''' 'c' 'a*b'});


    % BOOT - 1000 samples
    % run separately because of min p-value across vars thing for boot samp.
    % needed, and to get time

    % ---------------------------------------------------------------
    summary.null_p_boot = zeros(niter, 5);

    tt = clock;
    for i = 1:niter
        if mod(i, 50) == 0, fprintf('%3.0f ', i); end
        if mod(i, 1000) == 0, fprintf('\n'); end


        [paths, stats] = mediation(X(:, i), Y(:, i), M(:, i), 'boot');

        summary.null_p_boot(i, :) = stats.p;

    end

    et = etime(clock, tt);

    summary.fpr_boot = sum(summary.null_p_boot < .05) ./ niter;
    summary.runtime_boot = et ./ niter;

    fprintf(1,'\nBOOT\t');
    print_matrix(summary.fpr_boot, {'a' 'b' 'c''' 'c' 'a*b'});


end % n = N


%%
fprintf(1,'\nFalse-positive rates, p < .05\n-----------------------------------\n\n');
fprintf(1,'OLS\t');
print_matrix(summary.fpr_ols, {'a' 'b' 'c''' 'c' 'a*b'});

fprintf(1,'\nBOOT\t');
print_matrix(summary.fpr_boot, {'a' 'b' 'c''' 'c' 'a*b'});


plot_output;

% end main function




% _________________________________________________________________________
%
%
%
% * Inline functions
%
%
%
%__________________________________________________________________________


    function print_iteration_output

        persistent str timestr

        if i == 1
            verbstr = 'noverbose';
            str = sprintf(' %05d',1);
            fprintf(1,['Iteration ' str]);  
        end

        if i > 1

            if mod(i,10) == 1
                erase_string(timestr);
            end

            erase_string(str);
            str = sprintf(' %05d',i);  % iteration number
            fprintf(1, str);
        end

        if mod(i,10) == 0
            timestr = sprintf('Avg. time per iteration: %05d',etime(clock,t1) ./ i); fprintf(1,timestr);
        end

    end


    function plot_output
        create_figure('a_b', 1, 2);
        plot(max(summary.null_p_boot(:, 1:2), [], 2), summary.null_p_boot(:, 5), 'ko'); refline
        axis equal
        hold on; plot([0 1], [0 1], 'r')
        xlabel('Max of a and b p-values');
        ylabel('Mediation (a*b) p-value');

        subplot(1, 2, 2);
        plot(summary.null_p_boot(:, 1), summary.null_p_boot(:, 2), 'ko');
        wh = find(summary.null_p_boot(:, 5) < .05);
        hold on;
        plot(summary.null_p_boot(wh, 1), summary.null_p_boot(wh, 2), 'ko', 'MarkerFaceColor', 'r');
        plot_horizontal_line(.05); plot_vertical_line(.05)
        set(gca, 'XLim', [0 .1], 'YLim', [0 .1]);
        xlabel('a effect p-value'); ylabel('b effect p-value');
        title('Red = a*b p < .05');
    end

end  % main function









function erase_string(str1)
fprintf(1,repmat('\b',1,length(str1))); % erase string
end

