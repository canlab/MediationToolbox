function  bootsam = setup_boot_samples(x, bootsamples)
    % bootsam = setup_boot_samples(x, bootsamples)
    %
    % Tor Wager, Sept 2009
    %
    % Sets up a series of bootsamples bootstrap samples
    % Checks to see that there is x (predictor/input variable)
    % variance in each sample, and replaces samples with no
    % variance.  This is done because some analyses with small
    % samples and categorical predictors (e.g., contrast or effect coded variables)
    % will have bootstrap samples with no variance in
    % an input argument (x).  This could cause warnings,
    % errors, or incorrect results depending on what function
    % is subjected to bootstrapping.
    % This is used in mediation.m to make sure mediation is
    % not done on samples with no predictor (x variable)
    % variance.
    %
    % Also, the random number generator is initialized here.
    % there are problems with randomization if not done!
    % (at least in Matlab 7.5)


    % set up boot samples; make sure all are valid (otherwise
    % warnings/problems with categorical predictors + small
    % samples)
    %
    % initalize random number generator to new values; bootstrp uses this
    % there are problems with randomization if not done!
    rand('twister', sum(100*clock))
    n = size(x, 1);

    bootsam = ceil(n*rand(n, bootsamples)); % x value is the same for all boot samples

    maxiter = 500;
    
    % check and replace invalid ones with no variance in predictor
    % could do this later for mediator(s)/outcome, but may slow
    % down/pose other problems?
    wh_bad = ~any(diff(x(bootsam))); %| ~any(diff(m(bootsam)))
    icount = 1;
    
    while any(wh_bad)
        bootsam(:, wh_bad) = [];
        %ix = size(bootsam, 2) + 1;
        ntoadd = bootsamples - size(bootsam, 2);

        bootsam = [bootsam ceil(n*rand(n, ntoadd))]; % fill out rest of samples

        wh_bad = ~any(diff(x(bootsam)));
        icount = icount + 1;
        
        if icount > maxiter
            warning('MEDIATION:cannotGetBootSamples', 'Cannot get bootstrap samples; max iterations exceeded.  Sample size too small with categorical predictors?');
            bootsam(:, wh_bad) = [];
            break
        end
    end

end