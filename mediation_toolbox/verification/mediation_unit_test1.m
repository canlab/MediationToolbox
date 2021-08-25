function mediation_unit_test1

% Random single-level mediation with random variables
% Tests with/without NaNs and with/without bad data from a subject
% Varying sample sizes within-person
% some subjects will have too few trials by design
% some subjects will have NaN values by design

for N = [3 3 3 5 5 5 20 20 20]
    
    iter = 1;
    
    [X, Y, M, covs] = deal([]);
        
        X = rand(N, 1);
        
        M = rand(N, 1);
        
        Y = rand(N, 1);
        
        covs = rand(N, 2);  % up to two covariates
        
    
    if mod(iter, 3) == 0
        % no NaNs
    elseif mod(iter, 3) == 1
        % some NaNs
        
        n = numel(X);
        nnan = floor(n ./ 5);
        
        wh = randperm(n);
        wh = wh(1:nnan);
        X(wh) = NaN;
        
        wh = randperm(n);
        wh = wh(1:nnan);
        M(wh) = NaN;
        
        wh = randperm(n);
        wh = wh(1:nnan);
        Y(wh) = NaN;
        
        
                
    else
        % more NaNs
        
        n = numel(X);
        nnan = floor(n ./ 2);
        
        wh = randperm(n);
        wh = wh(1:nnan);
        X(wh) = NaN;
        
        wh = randperm(n);
        wh = wh(1:nnan);
        M(wh) = NaN;
        
        wh = randperm(n);
        wh = wh(1:nnan);
        Y(wh) = NaN;
        
    end
    
    run_models;
    
    iter = iter + 1;  % increment count for NaN purposes
    
end % for




% Random multilevel mediation with random variables
% Tests with/without NaNs and with/without bad data from a subject
% Varying sample sizes within-person
% some subjects will have too few trials by design
% some subjects will have NaN values by design

for N = [3 3 3 5 5 5 20 20 20]
    
    iter = 1;
    
    X = {};
    Y = {};
    M = cell(1, N);
    
    for i = 1:N
        
        ntrials = randperm(12);
        ntrials = ntrials(1);
        
        X{i} = rand(12, 1);
        
        M{i} = rand(12, 1);
        
        Y{i} = rand(12, 1);
        
    end
    
    
    if mod(iter, 3) == 0
        % no NaNs
    elseif mod(iter, 3) == 1
        % some NaNs
        
        add_random_nans;
        
    else
        % more NaNs
        
        for i = 1:5
            add_random_nans;
        end
        
    end
    
    run_models;
    
    iter = iter + 1;  % increment count for NaN purposes
    
end % for


% nested functions

    function add_random_nans()
        
        whsub = randperm(N);
        whsub = whsub(1);
        wh = randperm(length(X{whsub}));
        
        X{whsub}(wh(1)) = NaN;
        
        whsub = randperm(N);
        whsub = whsub(1);
        wh = randperm(length(M{whsub}));
        
        M{whsub}(wh(1)) = NaN;
        
        whsub = randperm(N);
        whsub = whsub(1);
        wh = randperm(length(X{whsub}));
        
        Y{whsub}(wh(1)) = NaN;
        
    end


    function run_models()
        
        [paths, stats] = mediation(X, Y, M, 'plots', 'verbose', 'boot'); %, 'summarystats');
        
        [paths, stats] = mediation(X, Y, M, 'plots', 'verbose', 'summarystats');
        
        [paths, stats] = mediation(X, Y, M, 'plots', 'verbose', 'boot', 'summarystats');
        
        [paths, stats] = mediation(X, Y, M, 'plots', 'verbose', 'boot', 'bootsamples', 5000); %, 'summarystats');
        
        [paths, stats] = mediation(X, Y, M);
        
        [paths, stats] = mediation(X, Y, M, 'names', {'x' 'y' 'm'});
        
        [paths, stats] = mediation(X, Y, M, 'names', {'x' 'y' 'm'}, 'signperm', 'verbose');
        
        [paths, stats] = mediation(X, Y, M, 'names', {'x' 'y' 'm'}, 'signperm', 'verbose', 'summarystats');
        
    end


end % Main function
