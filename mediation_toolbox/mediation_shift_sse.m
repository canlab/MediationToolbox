function [totalsse paths] = mediation_shift_sse(d, x, y, m, intcpt)

    %%% **** action item: make subfunction; persistent px pmx
    % if search y, save px, pmx, a, c
    % if search x, save nothing
    % if search m, save px, c
    % Check whether grouping px * [m y] is faster
    % check how regress, glmfit calculate betas, sterrs --

    n = size(x,1);
    if n ~= size(y,1) || n ~= size(m,1), error('Data vectors are not equal length.'); end
    if n ~= size(intcpt,1), error('Intercept vector is wrong length.'); end
    
    % shift variables x and m, based on delays d relative to y
    x = shift_variable(x, d(1));
    m = shift_variable(m, d(2));


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

