% MC_FWE = mediation_permutation_svc_fwe(X, Y, dat, nperms, varargin)
% % Takes all standard options that mediation.m does
%
% Single-level mediation only!
%
% Permutation test to assess familywise error rate corrected p-values for
% mediation test in a region of interest.
% 
% Permutes rows of data matrix dat and conducts mediation test on permuted
% rows.  Columns of dat are variables (e.g. voxels in the brain).
%
% Entering a single column of dat the single-mediation p-values.
% Entering multiple columns provides corrected p-values across tests (each
% col is a test) accounting for the correlations across tests.
%
% EXAMPLES:
% ----------------------------------------------------------------------
% cl = mask2clusters('mask_spm2_amy.img')
% cl = tor_extract_rois(SETUP.M, cl);
% 
% dat = cl(1).all_data;
% MC_FWE = mediation_permutation_svc_fwe(X, Y, dat, 100); % no bootstrap
% save FWEsim_L_amy_ROI_noboot MC_FWE
%
% Example : single voxel test (no MC)
% mediation(SETUP.X, SETUP.Y, dat, 'boot', 'plots');
% MC_FWE = mediation_permutation_svc_fwe(X, Y, dat, 1000, 'boot'); % bootstrap

function MC_FWE = mediation_permutation_svc_fwe(X, Y, dat, nperms, varargin)

% Takes all standard options that mediation.m does
options = varargin;

%fhan = @(m) mediation_wrapper(SETUP.X, SETUP.Y, m, 'boot');
fhan = @(m) mediation_wrapper(X, Y, m, options{:});   

permindx = permute_setupperms(size(dat, 1), nperms);

tmp = triu(corrcoef(dat)); tmp = tmp - eye(size(tmp)); tmp = tmp(:); tmp(tmp == 0) = [];
fprintf('Mean Correlation among voxels in ROI is: %3.2f\n', mean(tmp)); 

[actualsig, actualb, actualp] = matrix_eval_function(dat, fhan);

primary_thresholds = {.001 .002 .005 .05};
    
% permute data
% ----------------------------------------------------------------------
tic
fprintf('%05d', 0);

for i = 1:nperms

    fprintf('\b\b\b\b\b%05d', i);
    
    permdat = dat(permindx(i, :)', :);


    [sig, b, p] = matrix_eval_function(permdat, fhan, 'noverb');


    minpdist(i, :) = min(p, [], 1);

    % for global null: highest p-value under null across mediation test ; actual must be
    % lower than 95th% (e.g.) of highest under null to be significant in
    % ANY TEST.  (global null: NO effects)
    lowest_maxp(i, :) = min( max(p(:, [1 2 5]), [], 2)  );  
    
    
    % for conjunction null: lowest p-value under null ; actual must be
    % lower than 95th% (e.g.) of lowest under null to be significant in ALL
    % TESTS
    minpdist_a_b_ab_union(i, :) = min(minpdist(i, [1 2 5]));

    for j = 1:length(primary_thresholds)
        expected_fpr{j}(i, :) = sum(p < primary_thresholds{j}, 1);

        expected_fpr_conj{j}(i, :) = sum(expected_fpr{j}(i, [1 2 5]), 2);

    end


end

toc

% ----------------------------------------------------------------------
% summarize
% ----------------------------------------------------------------------

p_thresh_fwe_corr = prctile(minpdist, 5);

p_conj_fwe_corr = prctile(minpdist_a_b_ab_union, 5);

p_thresh_fwe_10 = prctile(minpdist, 10);

p_conj_fwe_10 = prctile(minpdist_a_b_ab_union, 10);

p_globalconj_fwe_corr = prctile(lowest_maxp, 5);
p_globalconj_fwe_10 = prctile(lowest_maxp, 10);

for j = 1:length(primary_thresholds)
    % how many times do we get a false positive somewhere in the region?
    corrp_at_primary{j} = max(sum(expected_fpr{j} > 0) ./ nperms, 1 ./ nperms);

    corrp_at_primary_conj{j} = max(sum(expected_fpr_conj{j} > 0) ./ nperms, 1 ./ nperms);
    
    fwe_a_b_ab_global_conj{j} = sum(all(minpdist(:, [1 2 5]) < primary_thresholds{j}, 2)) ./ nperms;
    fwe_a_b_global_conj{j} = sum(all(minpdist(:, [1 2]) < primary_thresholds{j}, 2)) ./ nperms;
    fwe_a_b_conj{j} = sum(any(minpdist(:, [1 2]) < primary_thresholds{j}, 2)) ./ nperms;


end

fprintf('------------------------------------\n');
fprintf('FWE Mult Correction Report\n');
fprintf('Two-tailed p-values\n');
fprintf('------------------------------------\n');

fprintf('Voxels: %3.0f, Avg. r = %3.2f, Permutations: %3.0f\n', size(dat, 2), mean(tmp), nperms);

fprintf('\nPrimary p-value thresholds needed for p < .1 corrected\n');
fprintf('Effect:\ta\tb\tc\tc''\tab\ta+b+ab conj\n');
fprintf('Primary:\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t\n', p_thresh_fwe_10, p_conj_fwe_10, p_globalconj_fwe_10);

fprintf('\nPrimary p-value thresholds needed for p < .05 corrected\n');
fprintf('Effect:\ta\tb\tc\tc''\tab\ta+b+ab conj\n');
fprintf('Primary:\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t\n', p_thresh_fwe_corr, p_conj_fwe_corr, p_globalconj_fwe_corr);

for j = 1:length(primary_thresholds)
    fprintf('\nFWE Corrected p-values at primary p < %3.6f\n', primary_thresholds{j});
    fprintf('Effect:\ta\tb\tc\tc''\tab\ta+b+ab conj\ta+b conj\ta+b+ab global null conj\ta+b global null conj\n');
    fprintf('Corr. P-val:\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t\n', corrp_at_primary{j}, corrp_at_primary_conj{j}, ...
        fwe_a_b_conj{j}, fwe_a_b_ab_global_conj{j}, fwe_a_b_global_conj{j});

    fprintf('Expected FP vox:\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t\n', mean(expected_fpr{j}), mean(expected_fpr_conj{j}));
end

permstats = struct('permindx', permindx, 'minpdist', minpdist, 'minpdist_a_b_ab_union', minpdist_a_b_ab_union);
permstats.expected_fpr = expected_fpr;
permstats.expected_fpr_conj = expected_fpr_conj;
permstats.lowest_maxp = lowest_maxp;

MC_FWE = struct('actualb', actualb, 'actualp', actualp, 'p_thresh_fwe_corr', p_thresh_fwe_corr, 'p_conj_fwe_corr', p_conj_fwe_corr, ...
    'stats_by_permutations', permstats, 'nperms', nperms);
MC_FWE.corrp_at_primary = corrp_at_primary;
MC_FWE.corrp_at_primary_conj = corrp_at_primary_conj;
MC_FWE.p_thresh_fwe_10 = p_thresh_fwe_10; MC_FWE.p_conj_fwe_10 = p_conj_fwe_10;
MC_FWE.p_globalconj_fwe_corr = p_globalconj_fwe_corr;
MC_FWE.p_globalconj_fwe_10 = p_globalconj_fwe_10;
MC_FWE.fwe_a_b_ab_global_conj = fwe_a_b_ab_global_conj;
MC_FWE.fwe_a_b_global_conj = fwe_a_b_global_conj;
MC_FWE.fwe_a_b_conj = fwe_a_b_conj;

% Get sig values
% ----------------------------------------------------------------------

sig = MC_FWE.actualp < MC_FWE.p_conj_fwe_10; 
MC_FWE.sig.conjsig10 = all(sig(:, [1 2 5]), 2);

fprintf('\nSig. conjunction voxels at p < .10 corrected\n');
fprintf('Conj null: \t%3.0f voxels\tUniversal threshold applied, ALL tests significant\n', sum(MC_FWE.sig.conjsig10));

% Threshold each effect at its own corrected threshold
sig = MC_FWE.actualp < MC_FWE.p_conj_fwe_10;
MC_FWE.sig.asig = sig(:,1);
MC_FWE.sig.bsig = sig(:,2);
MC_FWE.sig.absig = sig(:, 5);
MC_FWE.sig.ind_conjsig10 = MC_FWE.sig.asig & MC_FWE.sig.bsig & MC_FWE.sig.absig;

fprintf('Conj null: \t%3.0f voxels\tCorrected threshold applied for each test, ALL tests significant\n', sum(MC_FWE.sig.ind_conjsig10));

sig = MC_FWE.actualp < MC_FWE.p_globalconj_fwe_10;
MC_FWE.sig.global_conjsig10 = all(sig(:, [1 2 5]), 2);

fprintf('Conj global null: \t%3.0f voxels\tUniversal threshold applied, AT LEAST ONE test significant\n', sum(MC_FWE.sig.global_conjsig10));


end



function [sig, b, p] = mediation_wrapper(X, Y, M, varargin)
    
    
    
    [paths, stats] = mediation(X, Y, M, varargin{:});
    
    b = nanmean(stats.paths, 1);
    p = nanmean(stats.p, 1);
    sig = p < .05;
    
    
end
