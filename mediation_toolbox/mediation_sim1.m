
% function mediation_sim1(iter, varargin)
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
% case {'logistic', 'logit'}, logistic = 'logit'; dologit = 1;
% case 'L2M', L2M = 'yes'; doL2M = 1;
%
%
% % Examples:
% -------------------------------------------------------------------------
%
% mediation_sim1(10000, 'bootsamples', 5000); mediation_sim_output_figs(5);
%
% Generate data with different noise sigmas for each 'participant'
% N = 15; sigma = linspace(3, 5, N);
% mediation_sim1(10000, 'N', N, 'sigma', sigma); mediation_sim_output_figs(5);
%
% logitstic (Wani)
% mediation_sim1(1000, 'bootsamples', 5000, 'noar', 'logistic'); mediation_sim_output_figs(5);
% 
% 2nd level moderation (Wani)
% mediation_sim1(1000, 'bootsamples', 5000, 'noar', 'L2M', 'sigma', .5, 'N', 50); mediation_sim_output_figs(1, [], 'L2M');


function mediation_sim1(iter, varargin)

    diary mediation_sim_output.txt

    % default parameters for data generation
    % -----------------------------------------------------------
    N = 15;         % number of subjects
    sigma = 2;
    n = 300;
    phi = [.5 .1];

    % hrf kernel
    % (later, vary across regions, etc., use real data...)
    hrfparams = [6 16 .9 1.2 10 0 30];
    hrf = spm_hrf(1, hrfparams);

    nullmeth = 'niall_null1';
    altmeth = 'niall_alt1';

    % method of generating data; null or alt. hypothesis
    % 'mnull' is null-hyp data, true x-y relationship, m is unrelated

    % default parameters for model fitting
    % -----------------------------------------------------------
    bootopt = 'boot';      % boot, or noboot, or signperm
    latentopt = 'nolatent';   % latent or otherwise for no latent
    shiftopt = 'noshift';   % shift or otherwise for no shift
    shiftvalsopt = [-3 3];  % shift vals, if opt is on
    arstropt = 'ar';           % ar or otherwise for no ar model
    aropt = 2;              % ar model order, or 0
    weightopt = 'hierarchical';  % hierarchical or summarystats
    bootsamples = 1000;     % bootstrap samples
    logistic = 'nologit'; dologit = 0;
    L2M = 'no'; doL2M = 0;

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}

                case {'boot'}, bootopt = 'boot';
                case 'noboot', bootopt = 'no';
                case 'signperm', bootopt = 'signperm';

                case 'ar', arstropt = 'ar';
                case 'noar', arstropt = 'no';

                case 'latent', latentopt = 'latent';
                case 'nolatent', latentopt = 'no';

                case 'hierarchical', weightopt = 'hierarchical';
                case 'summarystats', weightopt = 'summarystats';

                case 'bootsamples', bootsamples = varargin{i+1};
                case {'logistic', 'logit'}, logistic = 'logit'; dologit = 1;

                case 'shift', shiftopt = 'shift'; shiftvalsopt = varargin{i+1};

                case {'N'}, N = varargin{i+1};

                case 'sigma', sigma = varargin{i+1};
                case 'L2M', L2M = 'yes'; doL2M = 1;

                otherwise
                    fprintf('Unknown string input : %s\n', varargin{i});
            end
        end
    end

    if length(sigma) < N, sigma = repmat(sigma(1), 1, N); end

    % simulation process-control parameters
    % Load prior file or set up new params
    % -----------------------------------------------------------
    verbstr = 'verbose';
    t1 = clock();

    if exist(['.' filesep 'mediation_sim_output.mat'], 'file')
        load mediation_sim_output
        disp('Loaded existing mediation_sim_output');

    else
        if length(iter) > 1, error('Cannot start new sim file without entering single scalar number of iterations.'); end

        disp('Starting new simulation file.');


        nullmeans = NaN .* zeros(iter, 5);
        nullstes = NaN .* zeros(iter, 5);
        nullp = NaN .* zeros(iter, 5);

        altmeans = NaN .* zeros(iter, 5);
        altstes = NaN .* zeros(iter, 5);
        altp = NaN .* zeros(iter, 5);

        if strcmp(shiftopt, 'shift')
            nulldelays = NaN .* zeros(iter, 2);
            altdelays = NaN .* zeros(iter, 2);

            nulldelaystd = NaN .* zeros(iter, 2);
            altdelaystd = NaN .* zeros(iter, 2);

            nullconverged = NaN .* zeros(iter, 1);
            altconverged = NaN .* zeros(iter, 1);
        end

        if strcmp(shiftopt, 'latent')
            nullhrfparams = NaN .* zeros(iter, 2);
            althrfparams = NaN .* zeros(iter, 2);

            nullhrfstd = NaN .* zeros(iter, 2);
            althrfstd = NaN .* zeros(iter, 2);

            nullconverged = NaN .* zeros(iter, 1);
            altconverged = NaN .* zeros(iter, 1);
        end

    end

    % print set-up info to screen (and diary)
    % -----------------------------------------------------------

    fprintf('Set-up parameters for simulation:\n-----------------------------------\n\n');
    fprintf('Data generation:\n-----------------------------------\n');
    fprintf('N = %3.0f\n', N);
    fprintf('n = %3.0f\n', n);
    fprintf('sigma = %3.2f\n', sigma);
    fprintf('phi = %3.2f %3.2f\n', phi);
    fprintf('nullmeth = %s\n', nullmeth);
    fprintf('altmeth = %s\n', altmeth);
    fprintf('spm_hrf params = %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f\n', hrfparams);

    fprintf('\nModel fitting:\n-----------------------------------\n');
    fprintf('bootopt = %s\n', bootopt);
    fprintf('latentopt = %s\n', latentopt);
    fprintf('shiftopt = %s\n', shiftopt);
    fprintf('shiftvalsopt = %3.0f %3.0f\n', shiftvalsopt);

    fprintf('arstropt = %s\n', arstropt);
    fprintf('aropt = %3.0f\n', aropt);
    
    fprintf('logistic (Y) = %s\n', logistic);
    fprintf('2nd level moderator = %s\n', L2M);

    fprintf('weightopt = %s\n', weightopt);
    fprintf('nullmeth = %s\n', nullmeth);
    fprintf('\n-----------------------------------\n\n');

    % Set iterations
    % -----------------------------------------------------------
    if length(iter) == 1
        % entered a scalar number of iterations
        wh_iter = 1:iter;
    else
        % entered exact iteration numbers
        % Can append to saved file
        wh_iter = iter;
    end

    disp('Running iterations:')
    fprintf('%3.0f to %3.0f\n', wh_iter(1), wh_iter(end));



    for i = wh_iter

        % null-hypothesis data: x-y relation, no path thru m
        % ===========================================================

        % generate group data
        % -----------------------------------------------------------
        variables = generate_group_data(nullmeth, n, sigma, phi, hrf, N, dologit, doL2M);


        % fit model
        % -----------------------------------------------------------
        if ~doL2M
            [paths, stats2, wistats] = mediation(variables.X, variables.Y, variables.M, verbstr, bootopt, latentopt, shiftopt, shiftvalsopt, weightopt, logistic, arstropt, aropt, 'bootsamples', bootsamples);
        else
            [paths, stats2, wistats] = mediation(variables.X, variables.Y, variables.M, 'L2M', variables.L2M, verbstr, bootopt, latentopt, shiftopt, shiftvalsopt, weightopt, logistic, arstropt, aropt);%, 'bootsamples', bootsamples);
        end

        nullmeans(i,:) = stats2.mean(1,:);
        nullstes(i,:) = stats2.ste(1,:);
        nullp(i,:) = stats2.p(1,:);
        
        if doL2M
            if isequal(bootopt, 'boot'), nullmeans_L2M(i,:) = stats2.mean_L2M(2,:); 
            else nullmeans_L2M(i,:) = stats2.beta(2,:); end
            nullstes_L2M(i,:) = stats2.ste(2,:);
            nullp_L2M(i,:) = stats2.p(2,:);
        end

        if isfield(wistats, 'isconverged')
            nullconverged(i) = sum(wistats.isconverged) ./ N;
        end

        if isfield(wistats, 'delays')
            nulldelays(i,:) = mean(wistats.delays);
            nulldelaystd(i,:) = std(wistats.delays);
        end

        if isfield(wistats, 'hrfparams')
            for j = 1:N
                meanhrf(j,:) = mean(wistats.hrfparams{j}')';
            end
            nullhrfparams(i,:) = mean(meanhrf);
            nullhrfstd(i,:) = std(meanhrf);
        end


        % alt.-hypothesis data: x-y relation thru m
        % ===========================================================

        % generate group data
        % -----------------------------------------------------------
        
        variables = generate_group_data(altmeth, n, sigma, phi, hrf, N, dologit, doL2M);
        

        % fit model
        % -----------------------------------------------------------
        if ~doL2M
            [paths, stats2, wistats] = mediation(variables.X, variables.Y, variables.M, verbstr, bootopt, latentopt, shiftopt, shiftvalsopt, weightopt, logistic, arstropt, aropt, 'bootsamples', bootsamples);
        else
            [paths, stats2, wistats] = mediation(variables.X, variables.Y, variables.M, 'L2M', variables.L2M, verbstr, bootopt, latentopt, shiftopt, shiftvalsopt, weightopt, logistic, arstropt, aropt); %, 'bootsamples', bootsamples);
        end

        altmeans(i,:) = stats2.mean(1,:);
        altstes(i,:) = stats2.ste(1,:);
        altp(i,:) = stats2.p(1,:);
        
        if doL2M
            if isequal(bootopt, 'boot'), altmeans_L2M(i,:) = stats2.mean_L2M(2,:); 
            else altmeans_L2M(i,:) = stats2.beta(2,:); end
            altstes_L2M(i,:) = stats2.ste(2,:);
            altp_L2M(i,:) = stats2.p(2,:);
        end

        if isfield(wistats, 'isconverged')
            altconverged(i) = sum(wistats.isconverged) ./ N;
        end

        if isfield(wistats, 'delays')
            altdelays(i,:) = mean(wistats.delays);
            altdelaystd(i,:) = std(wistats.delays);
        end

        if isfield(wistats, 'hrfparams')
            for j = 1:N
                meanhrf(j,:) = mean(wistats.hrfparams{j}')';
            end
            althrfparams(i,:) = mean(meanhrf);
            althrfstd(i,:) = std(meanhrf);
        end

        print_iteration_output();

        if mod(i, 10) == 0
            % save
            if exist(['.' filesep 'mediation_sim_output.mat'], 'file')
                save mediation_sim_output -append null* alt*
            else
                save mediation_sim_output stats2 wistats null* alt* N sigma n phi hrf hrfparams *meth *opt
            end
        end


    end  % iterations


    fprintf('Final output:\n-----------------------------------\n\n');
    fprintf('\n\nNULL-HYPOTHESIS\n\n');

    fprintf('Mean path coefficients\n-----------------------------------\n\n');
    print_matrix(nullmeans)

    fprintf('2nd level p-values\n-----------------------------------\n\n');
    print_matrix(nullp)
    
    if doL2M
        fprintf('Mean path coefficients for 2nd level moderation\n-----------------------------------\n\n');
        print_matrix(nullmeans_L2M)
        
        fprintf('p-values for 2nd level moderation\n-----------------------------------\n\n');
        print_matrix(nullp_L2M)
    end
    
    fprintf('\n\nALT-HYPOTHESIS\n\n');

    fprintf('Mean path coefficients\n-----------------------------------\n\n');
    print_matrix(altmeans)

    fprintf('2nd level p-values\n-----------------------------------\n\n');
    print_matrix(altp)
    
    if doL2M
        fprintf('Mean path coefficients for 2nd level moderation\n-----------------------------------\n\n');
        print_matrix(altmeans_L2M)
        
        fprintf('p-values for 2nd level moderation\n-----------------------------------\n\n');
        print_matrix(altp_L2M)
    end
    
    % save
    if exist(['.' filesep 'mediation_sim_output.mat'], 'file')
        save mediation_sim_output -append null* alt*
    else
        save mediation_sim_output stats2 wistats null* alt* N sigma n phi hrf hrfparams *meth *opt
    end

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


    function print_iteration_output()
        persistent str timestr

        if i == wh_iter(1)
            verbstr = 'noverbose';
            str = sprintf(' %05d', 1);
            fprintf(['Iteration ' str]);
        else
            % later iteration; persistent vars already set

            if mod(i, 10) == 1
                erase_string(timestr);
            end

            erase_string(str);
            str = sprintf(' %05d', i);  % iteration number
            fprintf( str);
        end

        if mod(i, 10) == 0
            timestr = sprintf('Avg. time per iteration: %05d', etime(clock, t1) ./ i); fprintf(timestr);
        end

    end

end  % main function





function variables = generate_group_data(meth, n, sigma, phi, hrf, N, dologit, doL2M)
    % -----------------------------------------------------------
    % Generate data for N subjects
    %
    % If sigma is a vector, then different sigmas are used for diff.
    % subjects, to simulate case where multilevel (weighted) model is optimal
    % -----------------------------------------------------------

    % initialize random number generator
    randn('state', sum(100*clock))

    X = zeros(n, N);
    Y = zeros(n, N);
    M = zeros(n, N);
    
    if doL2M 
        sigma_L2M = 1;
        
        if isequal(meth, 'niall_null1'), b_l2m = [0 0 0]; meth_b = [0 0 1]; % 2nd level moderation effects, for [a b cp]
        elseif isequal(meth, 'niall_alt1'), b_l2m = [0 1 0]; meth_b = [1 1 0]; end
        
        M2L = normrnd(0, sigma_L2M, N, 3);
        B = repmat(b_l2m,N,1).*M2L + repmat(meth_b, N,1);
    
    else % ~doL2M
        B = NaN(N,3);
    end
    
    for i = 1:N
        [X(:,i), Y(:,i), M(:,i)] = get_data_single_subject(meth, n, sigma(i), phi, hrf, dologit, doL2M, B(i,:));
    end

    variables.X = X;
    variables.Y = Y;
    variables.M = M;
    if doL2M, variables.L2M = M2L(:,2); end
end



function [ox, oy, om] = get_data_single_subject(meth, n, sigma, phi, hrf, dologit, doL2M, B_L2M)
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
            r = mvnrnd([0 0], [sigma .5; .5 sigma], n + hrfartifactadjust);
            sx = r(:,1);
            sy = r(:,2);

            % m is unrelated
            % (alternatively, m could be related to x but not y)
            sm = normrnd(0, sigma, n + hrfartifactadjust, 1);


        case 'mtotal'
            % ===========================================================
            % alt-hypothesis data: x-y relation thru m
            % ===========================================================

            % **note: in some sense pathological, because true x-m-y is
            % caused by 2 different components.

            % sx, sy, sm are 'true, metabolic' signals in regions x, y, m
            % -----------------------------------------------------------
            % ** note: how correlated?
            r = mvnrnd([0 0], [sigma .5; .5 sigma], n + hrfartifactadjust);
            sx = r(:,1);
            sm = r(:,2);

            % y is related through m

            sy = (sm + normrnd(0, sigma, n + hrfartifactadjust, 1)) ./ 2;


        case 'mtotal2'
            % ===========================================================
            % alt-hypothesis data: x-y relation thru m
            % ===========================================================

            % **note: in some sense pathological, because true x-m-y is
            % caused by 2 different components.

            % sx, sy, sm are 'true, metabolic' signals in regions x, y, m
            % -----------------------------------------------------------
            % i think x-m and m-y are correlated .5/sqrt(2) for sigma = 1
            r = mvnrnd([0 0], [sigma .5; .5 sigma], n + hrfartifactadjust);
            sx = r(:,1);
            sm = r(:,2);

            % y is related through m
            r = mvnrnd([0 0], [sigma .5; .5 sigma], n + hrfartifactadjust);
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
            if ~doL2M
                mu =     [2 4 0 2.5 1 0];
            elseif doL2M
                mu =     [2 4 B_L2M(1) 2.5 B_L2M(3) B_L2M(2)];
            end
            

            % %             rsigma = eye(6) ./ 4;  % sigma_between!
            % %
            % %
            % %             B = mvnrnd(mu, rsigma, 1);

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

            if ~dologit
                sy = B(4) + B(6) * sm + B(5) * sx ;
                sy = sy + normrnd(0, sigma, n + hrfartifactadjust, 1) ;
            else
                sy = glmval([B(4); B(6); B(5)], [sm sx], 'logit');
                sy = sy > prctile(sy, 50);
            end
            
            sx = sx + normrnd(0, sigma, n + hrfartifactadjust, 1) ;
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
            if ~doL2M
                mu =     [2 4 1 2.5 0 1];
            elseif doL2M
                mu =     [2 4 B_L2M(1) 2.5 B_L2M(3) B_L2M(2)];
            end

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
            % %             B = mvnrnd(mu, rsigma, 1);

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

            if ~dologit
                sy = B(4) + B(6) * sm + B(5) * sx ;
                sy = sy + normrnd(0, sigma, n + hrfartifactadjust, 1) ;
            else
                sy = glmval([B(4); B(6); B(5)], [sm sx], 'logit');
                sy = sy > prctile(sy, 50);
            end
            
            sx = sx + normrnd(0, sigma, n + hrfartifactadjust, 1) ;
            sm = sm + normrnd(0, sigma, n + hrfartifactadjust, 1) ;
            
        otherwise
            error('Unknown data type in get_data_single_subject.');
    end

    % hx, hy, hm are noiseless hrf-convolved observed signals
    % -----------------------------------------------------------
    if ~dologit
        hx = fast_conv_fft(hrf, sx);
        hy = fast_conv_fft(hrf, sy);
        hm = fast_conv_fft(hrf, sm);

        % add noise coloring (increases sigma a little bit, but does not re-add
        % noise because we put error vector in as final argument)
        ox = noise_arp(n + hrfartifactadjust, phi, [], hx);
        oy = noise_arp(n + hrfartifactadjust, phi, [], hy);
        om = noise_arp(n + hrfartifactadjust, phi, [], hm);
    else
        ox = sx; om = sm; oy = sy;
    end
    % truncate to remove hrf conv artifact
    ox = ox(hrfartifactadjust+1 : end);
    oy = oy(hrfartifactadjust+1 : end);
    om = om(hrfartifactadjust+1 : end);

end


function erase_string(str1)
    fprintf(repmat('\b', 1, length(str1))); % erase string
end

