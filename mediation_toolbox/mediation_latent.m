function [paths, abetas, bbetas, cpbetas, cbetas, sterrs, intcpt, n, residm, residy, residy2, totalsse, hrfparams, hrf_xmy, isconverged] = ...
        mediation_latent(x,y,m,meth,domultilev,dorobust,boot1)
    %
    % [paths, abetas, bbetas, cpbetas, cbetas, sterrs, intcpt, n, residm, residy, residy2, totalsse, hrfparams, hrf_xmy, isconverged] = ...
    %   mediation_latent(x, y, m, 'ga', domultilev, dorobust, boot1)
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
    % [paths, totalsse,d,isconverged] = mediation_latent(x,y,m,[-6 6]); d
    %
    % tor wager, feb 2007

    intcpt = ones(size(x));
    
    % ----------------------------------------------------------------------
    % HRF basis set stuff
    % ----------------------------------------------------------------------    
    % define basis set here
    bforder = 2;
    bfstruct = spm_get_bf(struct('name','Gamma functions','dt',2,'length',30,'order',bforder));
    bf = bfstruct.bf;
    
    % starting values for HRFs x m y
    % range of -1 to 1 for each param
    % sobol start state bounds
    hrfparams = .0001 * ones(bforder,3);  %[1 1; 1 1; 1 1]';
    hrfparams(:,:,2) = ones(bforder,3);   %[1 1; 1 1; 1 1]';
    nbf = size(hrfparams,1);
    

    % ----------------------------------------------------------------------
    % options and defaults
    % ----------------------------------------------------------------------
    % defaults
    if nargin < 6, meth = 'ga'; end
    if nargin < 7, domultilev = 1; end
    if nargin < 8, dorobust = 0; end
    if nargin < 9, boot1 = 0; end
    
    isconverged = 0;

    objfun = @(d)mediation_latent_sse(d, x, y, m, intcpt,bf);

    
    % ----------------------------------------------------------------------
    % optimize (find solution for best delays)
    % ----------------------------------------------------------------------

    switch lower(meth)

        case 'ga'

            objfun_ga = @(hrfparams)-log(mediation_latent_sse(hrfparams, x, y, m, intcpt,bf));
            % need startlevels w/i each param ^ nparams gensize.
            [best_params,fit,beff,in,isconverged] = tor_ga(3^(3*nbf),20,{hrfparams},objfun_ga,'genconverge',3,'noverbose');
            hrfparams = best_params{1};
            totalsse = exp(-max(beff));

        case 'exhaustive'

            % NOT DEFINED YET FOR THIS FUNCTION
            
            % make exhaustive map
            stepsize = .5;

            [d, totalsse] = exhaustive_map(stepsize,drange,objfun);

            isconverged = 1;


        case 'combo'
            % NOT DEFINED YET FOR THIS FUNCTION

            % first low-res exhaustive
            stepsize = 2;
            [d, totalsse] = exhaustive_map(stepsize,drange,objfun);
            
            % then GA
            newrange = [min(drange(1),min(d) - 3)  max(drange(2),max(d) + 3)];
            
                        % sobol start state bounds
            start = [newrange(1) newrange(1)];
            start(:,:,2) = [newrange(2) newrange(2)];
            
            objfun_ga = @(d)-log(mediation_latent_sse(d, x, y, m, intcpt));
            [best_params,fit,beff,in,isconverged] = tor_ga(3^(3*nbf),20,{[newrange(1); newrange(2)]},objfun_ga,'genconverge',5,'noverbose');
            d = best_params{1};
            totalsse = exp(-max(beff));
            
            

        otherwise
            error('Unknown optimization method.')
    end


    
    % ----------------------------------------------------------------------
    % now do full model with AR, stes, etc. for the best delay
    % ----------------------------------------------------------------------
    
    % get hrfs for x, y, m from params
    [hrf_xmy x m y] = get_latent(hrfparams, x, m, y, bf);
    
    [paths, abetas, bbetas, cpbetas, cbetas, sterrs, intcpt, n, residm, residy, residy2] = ...
        mediation_path_coefficients(x, y, m, domultilev, dorobust, boot1);
    
end






    function [hrf_xmy x m y] = get_latent(hrfparams, x, m, y, bf)
    % get hrfs for x, y, m from params
    % constrain params to be all positive
    hrfparams(hrfparams <= 0) = .0001;
  
    hrf_xmy = bf * hrfparams;

    % norm so that equal area under curve for each HRF.
    % note: this means there is no unique param solution for hrfparams, but
    % rather relative param values.  also, diff hrf param values affect
    % variance of latent x/ y / m, so this introduces some uncertainty in
    % scaling.  maybe standardizing x/y/m, below, is better?
    hrf_xmy = hrf_xmy ./ repmat(sum(hrf_xmy),size(hrf_xmy,1),1);

    % deconvolve all vars to get latent
    x = fast_conv_fft(hrf_xmy(:,1),x, 'deconv');
    m = fast_conv_fft(hrf_xmy(:,2),m, 'deconv');
    y = fast_conv_fft(hrf_xmy(:,3),y, 'deconv');

    % scale vars
    x = scale(x);
    m = scale(m);
    y = scale(y);

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