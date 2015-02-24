% [a, b, c1, c, ab, aste, bste, c1ste, cste, abste, ap, bp, c1p, cp, abp, ...
% aind, bind, c1ind, cind, abind, aiste, biste, c1iste, ciste, abiste] = ...
% glmfit_brain_multilev_wrapper(dmpfc, hr, pag, 'boot');
%
% Wrapper function for mediation_brain_multilev
% Returns separate outputs for each variable that deserves a separate image name
%
% Tor Wager, Jan 2010
%
% Outputs:
% % names = {'intercept_b.img' 'slope_b.img' 'intercept_rfxvar_b.img' 'slope_rfxvar_b.img' ...
% %     'intercept_indiv.img' 'slope_indiv.img' ...
% %     'sigma.img' 'isconverged.img' 'intercept_t.img' 'slope_t.img'  ...
% %     'intercept_rfxvar_t.img' 'slope_rfxvar_t.img' 'intercept_p.img' 'slope_p.img'  ...
% %     'intercept_rfxvar_p.img' 'slope_rfxvar_p.img' };
%
% All this function does is run igls.m on whatever you put in, and then
% separates the output from a structure into separate output variables,
% which will then be written as separate image files.

function [intercept_b, slope_b, intercept_rfxvar_b, slope_rfxvar_b, intercept_indiv, slope_indiv  ...
        sigma, isconverged, intercept_t, slope_t, intercept_rfxvar_t, slope_rfxvar_t, ...
        intercept_p, slope_p, intercept_rfxvar_p, slope_rfxvar_p, intercept_LRTrfxvar_p, slope_LRTrfxvar_p, ...
        cov_b, cov_t, cov_p] = ...
        igls_brain_multilev_wrapper(X, Y, varargin)

    % main function call
    out = glmfit_multilevel(Y, X1, X2, varargin{:});
    
    % parse outputs into separate variables, for writing to separate images
    intercept_b = out.beta(1);
    slope_b = out.beta(2);
    
    intercept_rfxvar_b = out.betastar(1);
    slope_rfxvar_b = out.betastar(2);
    
    intercept_indiv = out.beta_indiv(1, :)';
    slope_indiv = out.beta_indiv(2, :)';
    
    sigma = out.Sigma;
    
    isconverged = out.isconverged;
    
    intercept_t = out.t(1);
    slope_t = out.t(2);
    
    % NOTE: the t-values / p-values for rfx are  different from the likelihood ratio test (LRT)!
    % the LRT p is the preferred one
    intercept_rfxvar_t = out.t_randvariance(1);
    slope_rfxvar_t = out.t_randvariance(2);
    
    intercept_p = out.p(1);
    slope_p = out.p(2);
    
    intercept_rfxvar_p = out.p_randvariance(1);
    slope_rfxvar_p = out.p_randvariance(2);
    
    intercept_LRTrfxvar_p = out.pLRT_randvariance(1);
    slope_LRTrfxvar_p = out.pLRT_randvariance(2);

    % Outputs for the btwn-subjects covariate
    % igls Works now for only ONE covariate!
    if length(out.beta) > 2
        cov_b = out.beta(3:end);
        cov_t = out.t(3:end);
        cov_p = out.p(3:end);

    else
        cov_b = 0;
        cov_t = 0;
        cov_p = 0;
    end
    
end




       