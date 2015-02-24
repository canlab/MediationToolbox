
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

function mediation_sim_single_level1(iter,varargin)
    
    %diary mediation_sim_output.txt
%%
% OLS vs. BOOTSTRAP vs. SIGNPERM

    % default parameters for data generation
    % -----------------------------------------------------------
    N = [5:5:40];         % number of subjects
    sigma = 2;
    niter = 100;

    clear summary
    summary.N = N;
    
    for n = 1:length(N)

    % Generate data: Null hypothesis
    % ---------------------------------------------------------------
    X = ones(N(n), 1);

    randn('state',sum(100*clock))
    Y = normrnd(0, sigma, N(n), niter);

    % OLS
    % ---------------------------------------------------------------
    tt = clock; stats = glmfit_general(Y, X); et = etime(clock, tt);

    summary.fpr_ols(n) = sum(stats.p < .05) ./ niter;
    summary.runtime_ols(n) = et ./ niter;

    
    % BOOT - 1000 samples
    % run separately because of min p-value across vars thing for boot samp.
    % needed

    % ---------------------------------------------------------------
    summary.p_boot = zeros(length(N), niter);
    summary.time_boot = zeros(length(N), niter);
    
    for i = 1:niter
        if mod(i, 50) == 0, fprintf('%3.0f ', i); end
        if mod(i, 1000) == 0, fprintf('\n'); end
        
        tt = clock;
        stats = glmfit_general(Y(:, i), X, 'boot');
        et = etime(clock, tt);
        summary.p_boot(n, i) = stats.p;
        summary.time_boot(n, i) = et;
    end
    summary.fpr_boot(n) = sum(summary.p_boot(n, :) < .05) ./ niter;
    summary.runtime_boot(n) = mean(summary.time_boot(n, :));


    % SIGNPERM - 1000 samples
    % run separately because of perm setup thing

    % ---------------------------------------------------------------
    summary.p_sign = zeros(length(N), niter);
    summary.time_sign = zeros(length(N), niter);
    
    for i = 1:niter
        if mod(i, 50) == 0, fprintf('%3.0f ', i); end
        if mod(i, 1000) == 0, fprintf('\n'); end
        
        tt = clock;
        stats = glmfit_general(Y(:, i), X, 'signperm');
        et = etime(clock, tt);
        summary.p_sign(n, i) = stats.p;
        summary.time_sign(n, i) = et;
    end
    summary.fpr_sign(n) = sum(summary.p_sign(n, :) < .05) ./ niter;
    summary.runtime_sign(n) = mean(summary.time_sign(n, :));
    
    end % n = N
    
%%

    fprintf(1,'Final output:\n-----------------------------------\n\n');
    fprintf(1,'\n\nNULL-HYPOTHESIS\n\n');

    fprintf(1,'Mean path coefficients\n-----------------------------------\n\n');
    print_matrix(nullmeans)

    fprintf(1,'2nd level p-values\n-----------------------------------\n\n');
    print_matrix(nullp)

    fprintf(1,'\n\nALT-HYPOTHESIS\n\n');

    fprintf(1,'Mean path coefficients\n-----------------------------------\n\n');
    print_matrix(altmeans)

    fprintf(1,'2nd level p-values\n-----------------------------------\n\n');
    print_matrix(altp)

    save mediation_sim_output out null* alt* N sigma n phi hrf hrfparams *meth *opt

    diary off

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
            fprintf(1,['Iteration ' str]);  end

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

end  % main function






function [X, Y, M] = generate_group_data(meth,n,sigma,phi,hrf,N)
    % -----------------------------------------------------------
    % Generate data for N subjects
    %
    % If sigma is a vector, then different sigmas are used for diff.
    % subjects, to simulate case where multilevel (weighted) model is optimal
    % -----------------------------------------------------------

    % initialize random number generator
    randn('state',sum(100*clock))

    X = zeros(n,N);
    Y = zeros(n,N);
    M = zeros(n,N);

    for i = 1:N
        [X(:,i), Y(:,i), M(:,i)] = get_data_single_subject(meth,n,sigma(i),phi,hrf);
    end

end



function [ox, oy, om] = get_data_single_subject(meth,n,sigma,phi,hrf)
    % -----------------------------------------------------------
    % Generate data for one subject
    % -----------------------------------------------------------

    hrfartifactadjust = round(length(hrf) ./ 2);
                
    switch meth

        case 'mnull'
            % ===========================================================
            % null-hypothesis data: x-y relation, no path thru m
            % ===========================================================

            % sx, sy, sm are 'true, metabolic' signals in regions x, y, m
            % -----------------------------------------------------------
            % x and y are correlated .5 ./ sigma
            r = mvnrnd([0 0],[sigma .5; .5 sigma],n + hrfartifactadjust);
            sx = r(:,1);
            sy = r(:,2);

            % m is unrelated
            % (alternatively, m could be related to x but not y)
            sm = normrnd(0,sigma,n + hrfartifactadjust,1);


        case 'mtotal'
            % ===========================================================
            % alt-hypothesis data: x-y relation thru m
            % ===========================================================

            % **note: in some sense pathological, because true x-m-y is
            % caused by 2 different components.

            % sx, sy, sm are 'true, metabolic' signals in regions x, y, m
            % -----------------------------------------------------------
            % ** note: how correlated?
            r = mvnrnd([0 0],[sigma .5; .5 sigma],n + hrfartifactadjust);
            sx = r(:,1);
            sm = r(:,2);

            % y is related through m

            sy = (sm + normrnd(0,sigma,n + hrfartifactadjust,1)) ./ 2;


        case 'mtotal2'
            % ===========================================================
            % alt-hypothesis data: x-y relation thru m
            % ===========================================================

            % **note: in some sense pathological, because true x-m-y is
            % caused by 2 different components.

            % sx, sy, sm are 'true, metabolic' signals in regions x, y, m
            % -----------------------------------------------------------
            % i think x-m and m-y are correlated .5/sqrt(2) for sigma = 1
            r = mvnrnd([0 0],[sigma .5; .5 sigma],n + hrfartifactadjust);
            sx = r(:,1);
            sm = r(:,2);

            % y is related through m
            r = mvnrnd([0 0],[sigma .5; .5 sigma],n + hrfartifactadjust);
            sy = r(:,2);
            sm = (sm + r(:,1)) ./ 2;

            
             case 'niall_null1'
            % ===========================================================
            % null hypothesis data: x-y direct, no relation thru m
            % ===========================================================

            % here we generate betas for each effect, a, b, and cp
            % and add noise determined by sigma
            %
            % random effects are:   a, b, cp, x intercept, m intercept, y intercept
            % x intercept is not fit by model

            % random coefficients (all uncorrelated)
            %    B = [x0 m0 a y0 cp b]
            mu =     [2 4 0 2.5 1 0];

% %             rsigma = eye(6) ./ 4;  % sigma_between!
% % 
% % 
% %             B = mvnrnd(mu,rsigma,1);
            
            B = mu;
            
%             effect_sizes = B([3 6 5]) ./ sigma;  % for this subject, [a b c']
%             effect_sizes = mu([3 6 5]) ./ sigma;  % population, [a b c']

            % sx, sy, sm are 'true, metabolic' signals in regions x, y, m
            % -----------------------------------------------------------
            % 
            
            % signal var in sx determines effect size
            % (higher = more signal = more power)
            xsignal = normrnd(0, 1, n + hrfartifactadjust, 1);
            
            sx = B(1) + xsignal;
            
            sm = B(2) + B(3) * sx;
            
            sy = B(4) + B(6) * sm + B(5) * sx ;
            
            sx = sx + normrnd(0, sigma, n + hrfartifactadjust, 1) ;
            
            sy = sy + normrnd(0, sigma, n + hrfartifactadjust, 1) ;
            
            sm = sm + normrnd(0, sigma, n + hrfartifactadjust, 1) ;
            
            
            
            
            case 'niall_alt1'
            % ===========================================================
            % alt-hypothesis data: x-y relation thru m
            % ===========================================================

            % here we generate betas for each effect, a, b, and cp
            % and add noise determined by sigma
            %
            % random effects are:   a, b, cp, x intercept, m intercept, y intercept
            % x intercept is not fit by model
            
            
            % random coefficients (all uncorrelated)
            %    B = [x0 m0 a y0 cp b]
            mu =     [2 4 1 2.5 0 1];

% %             rsigma = eye(6) ./ 4;  % sigma_between!
% % 
% % 
% %             % % sigma = [1.00          0.00          0.00          0.00          0.00          0.00;
% %             % %          0.00          1.00          0.00          0.50          0.00          0.00;
% %             % %          0.00          0.00          1.00          0.00          0.10         -0.50;
% %             % %          0.00          0.50          0.00          1.00          0.00          0.00;
% %             % %          0.00          0.00          0.10          0.00          0.15          0.10;
% %             % %          0.00          0.00         -0.50          0.00          0.10          1.00];
% % 
% %             B = mvnrnd(mu,rsigma,1);
            
            B = mu;
            
            
%             effect_sizes = B([3 6 5]) ./ sigma;  % for this subject, [a b c']
%             effect_sizes = mu([3 6 5]) ./ sigma;  % population, [a b c']

            % sx, sy, sm are 'true, metabolic' signals in regions x, y, m
            % -----------------------------------------------------------
            % 
            
            % signal var in sx determines effect size
            % (higher = more signal = more power)
            xsignal = normrnd(0, 1, n + hrfartifactadjust, 1);
            
            sx = B(1) + xsignal;
            
            sm = B(2) + B(3) * sx;
            
            sy = B(4) + B(6) * sm + B(5) * sx ;
            
            sx = sx + normrnd(0, sigma, n + hrfartifactadjust, 1) ;
            
            sy = sy + normrnd(0, sigma, n + hrfartifactadjust, 1) ;
            
            sm = sm + normrnd(0, sigma, n + hrfartifactadjust, 1) ;
            
            
        otherwise
            error('Unknown data type in get_data_single_subject.');
    end

    % hx, hy, hm are noiseless hrf-convolved observed signals
    % -----------------------------------------------------------
    hx = fast_conv_fft(hrf,sx);
    hy = fast_conv_fft(hrf,sy);
    hm = fast_conv_fft(hrf,sm);

    % add noise coloring (increases sigma a little bit, but does not re-add
    % noise because we put error vector in as final argument)
    ox = noise_arp(n + hrfartifactadjust, phi, [], hx);
    oy = noise_arp(n + hrfartifactadjust, phi, [], hy);
    om = noise_arp(n + hrfartifactadjust, phi, [], hm);

    % truncate to remove hrf conv artifact
    ox = ox(hrfartifactadjust+1 : end);
    oy = oy(hrfartifactadjust+1 : end);
    om = om(hrfartifactadjust+1 : end);
    
end


function erase_string(str1)
    fprintf(1,repmat('\b',1,length(str1))); % erase string
end

