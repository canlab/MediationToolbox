% [a, b, c1, c, ab, aste, bste, c1ste, cste, abste, ap, bp, c1p, cp, abp, ...
% aind, bind, c1ind, cind, abind, aiste, biste, c1iste, ciste, abiste] = ...
% mediation_brain_multilev_wrapper(dmpfc, hr, pag, 'boot');
%
% Wrapper function for mediation_brain_multilev
% Returns separate outputs for each variable that deserves a separate image name
%
% Tor Wager, Oct 2007
%
% Outputs:
% % names = {'X-M_effect.img' 'M-Y_effect.img' 'X-Y_direct_effect.img' 'X-Y_total_effect.img' 'X-M-Y_effect.img' ...
% %     'X-M_ste.img' 'M-Y_ste.img' 'X-Y_direct_ste.img' 'X-Y_total_ste.img' 'X-M-Y_ste.img' ...
% %     'X-M_pvals.img' 'M-Y_pvals.img' 'X-Y_direct_pvals.img' 'X-Y_total_pvals.img' 'X-M-Y_pvals.img' ...
% %     'X-M_indiv_effect.img' 'M-Y_indiv_effect.img' 'X-Y_direct_indiv_effect.img' 'X-Y_total_indiv_effect.img' 'X-M-Y_indiv_effect.img' ...
% %     'X-M_indiv_ste.img' 'M-Y_indiv_ste.img' 'X-Y_direct_indiv_ste.img' 'X-Y_total_indiv_ste.img' 'X-M-Y_indiv_ste.img'};
%
% Additional names for level-2 moderators:
% 'a_L2mod' 'b_L2mod' c1_L2mod' 'c_L2mod' 'ab_L2mod' 'ap_L2mod' 'bp_L2mod' c1p_L2mod' 'cp_L2mod' 'abp_L2mod'
%
% Additional mediators are added as additional volumes in a, b, ab, and ap,
% bp, abp images.

function [a, b, c1, c, ab, aste, bste, c1ste, cste, abste, ap, bp, c1p, cp, abp, ...
        aind, bind, c1ind, cind, abind, aiste, biste, c1iste, ciste, abiste, ...
        a_L2mod, b_L2mod, c1_L2mod, c_L2mod, ab_L2mod, ap_L2mod, bp_L2mod, c1p_L2mod, cp_L2mod, abp_L2mod] = ...
        mediation_brain_multilev_wrapper(X, Y, M, varargin)

    [indiv_paths, toplevelstats, firstlevelstats] = mediation(X, Y, M, varargin{:});
    aind = indiv_paths(:,1);
    bind = indiv_paths(:,2);
    c1ind = indiv_paths(:,3);
    cind = indiv_paths(:,4);
    abind = indiv_paths(:,5);

    meanpath = toplevelstats.mean;
    a = meanpath(1);
    b = meanpath(2);
    c1 = meanpath(3);
    c = meanpath(4);
    ab = meanpath(5);

    stepath = toplevelstats.ste(1, :);  % for intercept in 2nd level model: group within-subjects effect
    aste = stepath(1);
    bste = stepath(2);
    c1ste = stepath(3);
    cste = stepath(4);
    abste = stepath(5);

    p = toplevelstats.p(1, :);  % for intercept in 2nd level model: group within-subjects effect;
    ap = p(1);
    bp = p(2);
    c1p = p(3);
    cp = p(4);
    abp = p(5);

    indiv_ste = firstlevelstats.ste;
    aiste = indiv_ste(1);
    biste = indiv_ste(2);
    c1iste = indiv_ste(3);
    ciste = indiv_ste(4);
    abiste = indiv_ste(5);


    % Additional mediators, if any
    naddm = 0;
    if iscell(toplevelstats.inputOptions.additionalM)
        naddm = size(toplevelstats.inputOptions.additionalM{1}, 2);
    end

    for i = 1:naddm

        newpaths = toplevelstats.beta(1, 5 + 3*i - 2 : 5 + 3*i);  % these are the a, b, and ab paths for this mediator
        [a(i + 1, 1) b(i + 1, 1) ab(i + 1, 1)] = deal(newpaths(1), newpaths(2), newpaths(3));

        newp = toplevelstats.p(1, 5 + 3*i - 2 : 5 + 3*i);  % these are the a, b, and ab paths for this mediator
        [ap(i + 1, 1) bp(i + 1, 1) abp(i + 1, 1)] = deal(newp(1), newp(2), newp(3));

    end

    [a_L2mod, b_L2mod, c1_L2mod, c_L2mod, ab_L2mod, ap_L2mod, bp_L2mod, c1p_L2mod, cp_L2mod, abp_L2mod] = deal(NaN);

    % 2nd level moderators, if any
    nmod = size(toplevelstats.p, 1) - 1;
    if nmod
        meanpath = toplevelstats.beta(2:end, :);
        a_L2mod = meanpath(:, 1);
        b_L2mod = meanpath(:, 2);
        c1_L2mod = meanpath(:, 3);
        c_L2mod = meanpath(:, 4);
        ab_L2mod = meanpath(:, 5);

        p = toplevelstats.p(2:end, :);   % 2nd-level moderators
        ap_L2mod = p(:, 1);
        bp_L2mod = p(:, 2);
        c1p_L2mod = p(:, 3);
        cp_L2mod = p(:, 4);
        abp_L2mod = p(:, 5);

        % moderation of additional mediators
        % not saved; would require reformatting outputs, which is a pain
        % %            for i = 1:naddm
        % %
        % %         newpaths = toplevelstats.beta(2:end, 5 + 3*i - 2 : 5 + 3*i);  % these are the a, b, and ab paths for this mediator
        % %         [a_L2mod(i + 1, 1) b_L2mod(i + 1, 1) ab_L2mod(i + 1, 1)] = deal(newpaths(:, 1), newpaths(:, 2), newpaths(:, 3));
        % %
        % %            end

    end


    %a = cat(2, toplevelstats.bbetas{:})';


    % % names = {'X-M_effect.img' 'M-Y_effect.img' 'X-Y_direct_effect.img' 'X-Y_total_effect.img' 'X-M-Y_effect.img' ...
    % %     'X-M_ste.img' 'M-Y_ste.img' 'X-Y_direct_ste.img' 'X-Y_total_ste.img' 'X-M-Y_ste.img' ...
    % %     'X-M_pvals.img' 'M-Y_pvals.img' 'X-Y_direct_pvals.img' 'X-Y_total_pvals.img' 'X-M-Y_pvals.img' ...
    % %     'X-M_indiv_effect.img' 'M-Y_indiv_effect.img' 'X-Y_direct_indiv_effect.img' 'X-Y_total_indiv_effect.img' 'X-M-Y_indiv_effect.img' ...
    % %     'X-M_indiv_ste.img' 'M-Y_indiv_ste.img' 'X-Y_direct_indiv_ste.img' 'X-Y_total_indiv_ste.img' 'X-M-Y_indiv_ste.img'};


end
