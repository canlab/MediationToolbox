function [totalsse paths hrf_xmy] = mediation_latent_sse(hrfparams,x, y, m, intcpt,bf)
    % [totalsse paths hrf_xmy] = mediation_latent_sse(hrfparams,x, y, m,intcpt,bf)
    %
    % hrfparams = [x1 x2; m1 m2; y1 y2]';
    %
    %%% **** action item: make subfunction; persistent px pmx
    % if search y, save px, pmx, a, c
    % if search x, save nothing
    % if search m, save px, c
    % Check whether grouping px * [m y] is faster
    % check how regress, glmfit calculate betas, sterrs --
    %
    % bfstruct = spm_get_bf(struct('name','Gamma functions','dt',2,'length',20,'order',2));
    % bf = bfstruct.bf;

    n = size(x,1);
    if n ~= size(y,1) || n ~= size(m,1), error('Data vectors are not equal length.'); end
    if n ~= size(intcpt,1), error('Intercept vector is wrong length.'); end

    % get hrfs for x, y, m from params
    [hrf_xmy x m y] = get_latent(hrfparams, x, m, y, bf);

    xx = [x intcpt];
    px = (xx'*xx)\xx';	                    % X beta-forming matrix
    mx = [m xx];
    pmx = (mx'*mx)\mx';                     % M+X beta-forming matrix

    % Is X related to the mediator?
    % Eq. 2, a = M ~ X
    betas = px * m;
    r = y - xx * betas;
    sse(1) = r' * r;
    a = betas(1);

    % Are M (b path) and X (c' path) independently related to Y?
    % Eq. 3, [b, c'] = Y ~ X + M
    betas = pmx * y;
    r = y - mx * betas;
    sse(2) = r' * r;
    b = betas(1);
    cp = betas(2);

    % is X related to Y without mediator?
    % Eq. 1, Y ~ X
    betas = px * y;
    r = y - xx * betas;
    sse(3) = r' * r;
    c = betas(1);

    ab = a .* b;

    totalsse = sum(sse);
    paths = [a b cp c ab];
end




function [hrf_xmy x m y] = get_latent(hrfparams, x, m, y, bf)
    % get hrfs for x, y, m from params
    % constrain params to be all positive
    hrfparams(hrfparams <= 0) = .0001;

    if size(bf,2) ~= size(hrfparams,1)
        error('Number of basis functions and number of parameters do not match!');
    end
    
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

