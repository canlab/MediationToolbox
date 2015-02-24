% function mediation_dcm_sim1(iter, varargin)
% Simulation for power and false positive rates for DCM analysis
%
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
% mediation_dcm_sim1(10000, 'bootsamples', 5000); mediation_sim_output_figs(5);
%
% Generate data with different noise sigmas for each 'participant'
% N = 15; sigma = linspace(3, 5, N);
% mediation_dcm_sim1(10000, 'N', N, 'sigma', sigma); mediation_sim_output_figs(5);

function mediation_dcm_sim1(iter, varargin)

    diary mediation_dcm_sim_output.txt

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

                case 'shift', shiftopt = 'shift'; shiftvalsopt = varargin{i+1};

                case 'N', N = varargin{i+1};
                case 'nscan', n = varargin{i+1};

                case 'sigma', sigma = varargin{i+1};

                otherwise error('Unknown string input.');
            end
        end
    end

    if length(sigma) < N, sigma = repmat(sigma(1), 1, N); end

    % simulation process-control parameters
    % -----------------------------------------------------------
    verbstr = 'verbose';
    time_per_it = 0;
    t1 = clock;

    nullmeans = NaN .* zeros(iter, 4);
    nullstes = NaN .* zeros(iter, 4);
    nullp = NaN .* zeros(iter, 4);

    altmeans = NaN .* zeros(iter, 4);
    altstes = NaN .* zeros(iter, 4);
    altp = NaN .* zeros(iter, 4);

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

    % print set-up info to screen (and diary)
    % -----------------------------------------------------------

    fprintf(1, 'Set-up parameters for simulation:\n-----------------------------------\n\n');
    fprintf(1, 'Data generation:\n-----------------------------------\n');
    fprintf(1, 'N = %3.0f\n', N);
    fprintf(1, 'n = %3.0f\n', n);
    fprintf(1, 'sigma = %3.2f\n', sigma);
    fprintf(1, 'phi = %3.2f %3.2f\n', phi);
    fprintf(1, 'nullmeth = %s\n', nullmeth);
    fprintf(1, 'altmeth = %s\n', altmeth);
    fprintf(1, 'spm_hrf params = %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f\n', hrfparams);

    fprintf(1, '\nModel fitting:\n-----------------------------------\n');
    fprintf(1, 'bootopt = %s\n', bootopt);
    fprintf(1, 'latentopt = %s\n', latentopt);
    fprintf(1, 'shiftopt = %s\n', shiftopt);
    fprintf(1, 'shiftvalsopt = %3.0f %3.0f\n', shiftvalsopt);

    fprintf(1, 'arstropt = %s\n', arstropt);
    fprintf(1, 'aropt = %3.0f\n', aropt);

    fprintf(1, 'weightopt = %s\n', weightopt);
    fprintf(1, 'nullmeth = %s\n', nullmeth);
    fprintf(1, '\n-----------------------------------\n\n');

    while(~exist('mat_file_name', 'var') || exist(mat_file_name, 'file'))
        mat_file_name = ['dcm-sim-results-' num2str(round(rand(1,1) * 10000))];
    end
    fprintf('Saving to %s (number suffix means nothing)\n', mat_file_name);
    
    for i = 1:iter


        % null-hypothesis data: x-y relation, no path thru m
        % ===========================================================

        % generate group data
        % -----------------------------------------------------------
        [X, Y, M] = generate_group_data(nullmeth, n, sigma, phi, hrf, N);
        psych_input = rand(n, 1);
        SPM = fakeSPM(n, psych_input);

        % fit model
        % -----------------------------------------------------------
        %         [paths, stats2, wistats] = mediation(X, Y, M, verbstr, bootopt, latentopt, shiftopt, shiftvalsopt, weightopt, arstropt, aropt, 'bootsamples', bootsamples);
        [paths, means, stes, pvals] = dcm_sim(SPM, psych_input, X, M, Y);


        nullmeans(i, :) = means;
        nullstes(i, :) = stes;
        nullp(i, :) = pvals;

        % alt.-hypothesis data: x-y relation thru m
        % ===========================================================

        % generate group data
        % -----------------------------------------------------------
        [X, Y, M] = generate_group_data(altmeth, n, sigma, phi, hrf, N);
        psych_input = rand(n, 1);
        SPM = fakeSPM(n, psych_input);


        % fit model
        % -----------------------------------------------------------
        %         [paths, stats2, wistats] = mediation(X, Y, M, verbstr, bootopt, latentopt, shiftopt, shiftvalsopt, weightopt, arstropt, aropt, 'bootsamples', bootsamples);
        [paths, means, stes, pvals] = dcm_sim(SPM, psych_input, X, M, Y);

        altmeans(i, :) = means;
        altstes(i, :) = stes;
        altp(i, :) = pvals;

        save(mat_file_name, 'null*', 'alt*', 'N', 'sigma', 'n', 'phi', 'hrf', 'hrfparams', '*meth', '*opt', 'i');
        
        print_iteration_output;
    end


    fprintf(1, 'Final output:\n-----------------------------------\n\n');
    fprintf(1, '\n\nNULL-HYPOTHESIS\n\n');

    fprintf(1, 'Mean path coefficients\n-----------------------------------\n\n');
    print_matrix(nullmeans)

    fprintf(1, '2nd level p-values\n-----------------------------------\n\n');
    print_matrix(nullp)

    fprintf(1, '\n\nALT-HYPOTHESIS\n\n');

    fprintf(1, 'Mean path coefficients\n-----------------------------------\n\n');
    print_matrix(altmeans)

    fprintf(1, '2nd level p-values\n-----------------------------------\n\n');
    print_matrix(altp)


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
            str = sprintf(' %05d', 1);
            fprintf(1, ['Iteration ' str]);  end

        if i > 1

            if mod(i, 10) == 1
                erase_string(timestr);
            end

            erase_string(str);
            str = sprintf(' %05d', i);  % iteration number
            fprintf(1, str);
        end

        if mod(i, 10) == 0
            timestr = sprintf('Avg. time per iteration: %05d', etime(clock, t1) ./ i); fprintf(1, timestr);
        end

    end

end  % main function






function [X, Y, M] = generate_group_data(meth, n, sigma, phi, hrf, N)
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

    for i = 1:N
        [X(:, i), Y(:, i), M(:, i)] = get_data_single_subject(meth, n, sigma(i), phi, hrf);
    end

end



function [ox, oy, om] = get_data_single_subject(meth, n, sigma, phi, hrf)
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
            sx = r(:, 1);
            sy = r(:, 2);

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
            sx = r(:, 1);
            sm = r(:, 2);

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
            sx = r(:, 1);
            sm = r(:, 2);

            % y is related through m
            r = mvnrnd([0 0], [sigma .5; .5 sigma], n + hrfartifactadjust);
            sy = r(:, 2);
            sm = (sm + r(:, 1)) ./ 2;


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
    hx = fast_conv_fft(hrf, sx);
    hy = fast_conv_fft(hrf, sy);
    hm = fast_conv_fft(hrf, sm);

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


function SPM = fakeSPM(n, psych_input)
    SPM.Sess.row = 1:n;
    SPM.xX.X = [psych_input ones(n, 1)];
    SPM.xX.xKXs.X = SPM.xX.X;
    SPM.xX.iB = 2;
    SPM.xX.iG = [];
    SPM.xX.K = struct('foo', 1);
    SPM.xY.RT = 2;
    SPM.Sess.U(1).dt = SPM.xY.RT / 16;
end