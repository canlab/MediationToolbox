function [paths, abetas, bbetas, cpbetas, cbetas, sterrs, intcpt, n, residm, residy, residy2, totalsse, d, isconverged] = ...
        mediation_shift(x,y,m,drange,meth,domultilev,dorobust,boot1)
    %
    % [paths, abetas, bbetas, cpbetas, cbetas, sterrs, n, residm, residy, residy2, totalsse, d, isconverged] = ...
    % = mediation_shift(x, y, m, drange, meth, domultilev, dorobust, boot1)
    %
    % OLD: [paths, totalsse,d,isconverged] = mediation_shift(x,y,m,intcpt,drange,meth,domultilev,dorobust,boot1)
    %
    % inputs: see mediation_path_coefficients and mediation.m
    % ---------------------------------------------------------------------
    %
    % meth = 'ga' or 'exhaustive' or 'combo'
    % 
    % outputs:
    % ---------------------------------------------------------------------
    % d = delay values
    %
    % examples:
    % ---------------------------------------------------------------------
    % [paths, totalsse,d,isconverged] = mediation_shift(x,y,m,[-6 6]); d
    %
    % tor wager, feb 2007

    global mediation_covariates
    
    intcpt = ones(size(x));
    
    % Remove covariates, if specified
    if ~isempty(mediation_covariates)
        covs = [mediation_covariates intcpt];
        x_pinv_covs = covs * pinv(covs);
        x = x - x_pinv_covs * x;
        y = y - x_pinv_covs * y;
        m = m - x_pinv_covs * m;
    end
        
    objfun = @(d)mediation_shift_sse(d, x, y, m, intcpt);

    % defaults
    if nargin < 6, meth = 'ga'; end
    if nargin < 7, domultilev = 1; end
    if nargin < 8, dorobust = 0; end
    if nargin < 9, boot1 = 0; end
    
    isconverged = 0;

    
    
    % range of delays possible
    %drange = [-20 20];

    % ----------------------------------------------------------------------
    % optimize (find solution for best delays)
    % ----------------------------------------------------------------------

    switch lower(meth)

        case 'ga'

            % with sobol start, this always got the right solution for me in the first
            % iteration
            
            % sobol start state bounds
            start = [drange(1) drange(1)];
            start(:,:,2) = [drange(2) drange(2)];

            objfun_ga = @(d)-log(mediation_shift_sse(d, x, y, m, intcpt));
            
            % start state should have resolution of TR/3
            % we need vals ^ params organisms to cover each point of a grid
            % at TR/3 spacing, so we need with 2 params:
            startorgs = (range(drange) .* 3) .^ 2;
            
            % 576 is 24^2, so start state will span space with 24 vals
            % within range.
            [best_params,fit,beff,in,isconverged] = tor_ga(startorgs,30,{start},objfun_ga,'genconverge',3,'noverbose');
            d = best_params{1};
            totalsse = exp(-max(beff));

        case 'exhaustive'


            % make exhaustive map
            stepsize = .5;

            [d, totalsse] = exhaustive_map(stepsize,drange,objfun);

            isconverged = 1;


        case 'combo'
            
            % first low-res exhaustive
            stepsize = 2;
            [d, totalsse] = exhaustive_map(stepsize,drange,objfun);
            
            % then GA
            newrange = [min(drange(1),min(d) - 3)  max(drange(2),max(d) + 3)];
            
                        % sobol start state bounds
            start = [newrange(1) newrange(1)];
            start(:,:,2) = [newrange(2) newrange(2)];
            
            objfun_ga = @(d)-log(mediation_shift_sse(d, x, y, m, intcpt));
            [best_params,fit,beff,in,isconverged] = tor_ga(200,20,{[newrange(1); newrange(2)]},objfun_ga,'genconverge',5,'noverbose');
            d = best_params{1};
            totalsse = exp(-max(beff));
            
            

        otherwise
            error('Unknown optimization method.')
    end


    %[totalsse paths] = objfun(d);

    d = d';
    
    % ----------------------------------------------------------------------
    % now do full model with AR, stes, etc. for the best delay
    % ----------------------------------------------------------------------
    
    % shift variables x and m, based on delays d relative to y
    x = shift_variable(x, d(1));
    m = shift_variable(m, d(2));
    
    [paths, abetas, bbetas, cpbetas, cbetas, sterrs, intcpt, n, residm, residy, residy2] = ...
        mediation_path_coefficients(x, y, m, domultilev, dorobust, boot1);
    
end



function [d, totalsse, ssemap] = exhaustive_map(stepsize,drange,objfun)


    d1range = [drange(1):stepsize:drange(2)]; d2range = [drange(1):stepsize:drange(2)];
    [d1,d2] = meshgrid(d1range,d2range);

    ssemap = NaN .* zeros(size(d1));

    [m,n] = size(d1);
    for j = 1:m
        for k = 1:n
            ssemap(j,k) = objfun([d1(j,k) d2(j,k)]);
        end
    end
    [w1,w2] = find(ssemap == min(ssemap(:)));
    totalsse = ssemap(w1,w2);
    best_d1 = d1(w1,w2);
    best_d2 = d2(w1,w2);
    d = [best_d1 best_d2];

end




function vec = shift_variable(vec, shift, interp_method)
    if(~exist('interp_method', 'var') || isempty(interp_method))
        interp_method = 'linear';
    end

    if ~iscol(vec), vec = vec'; end

    if shift == 0
        % we're done
        return

    elseif (mod(shift, 1) ~= 0)
        % non-integer shift value
        
        % define xi, new x values; works for both pos and neg shift
        xi = ((-shift) + 1) : (length(vec) - shift);
        shifted_vec = interp1(vec, xi, interp_method)';
        
        whnan = isnan(shifted_vec);
        vec = shifted_vec;
        vec(whnan) = mean(vec(~whnan));

    else
        % integer shift value
        n = length(vec);
        if shift > 0
            % shift forward, imputing mean to missing initial elements
            vec = [repmat(mean(vec), shift, 1); vec];
            vec = vec(1:n);

        elseif shift < 0
            % shift back, imputing mean to missing final elements
            vec = [vec((-shift)+1:end); repmat(mean(vec), -shift, 1)];
        end

    end


end