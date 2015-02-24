function mediation_threepaths_sim_output_figs(wh_col, file_to_load)
    % mediation_sim_output_figs(wh_col)
    %
    % goes with output from mediation_sim1
    % just go to directory and run

    if(~exist('file_to_load', 'var') || isempty(file_to_load))
        file_to_load = 'mediation_threepaths_sim_output';
    end
    
    load(file_to_load);

    %% setup

    nsubplots = 3;

    %%% wh_col = 5;
    names = {'b1'  'b2'  'b3' 'c'''  'c'  'b1b2b3'};


    pthr = [.0001:.0001:.99];
    npvals = length(pthr);
    n = size(nullp,1);

    hr = NaN .* zeros(npvals,1);
    far = NaN .* zeros(npvals,1);
    phit =  NaN .* zeros(npvals,1); % overall p(yes)

    % ========================================================
    %% run
    % ========================================================

    for i = 1:length(pthr)
        whnull = nullp(:, wh_col) < pthr(i);
        whalt = altp(:, wh_col) < pthr(i);

        hr(i,1) = sum(whalt) ./ n;
        far(i,1) = sum(whnull) ./ n;

        phit(i,1) = sum([whalt;whnull]) ./ (2*n);

        % %     roc(i,1) = sum(whnull) ./ n;
        % %     roc(i,2) = sum(whalt) ./ n;

    end

    % ========================================================
    %% ROC plot
    % ========================================================

    dprime = norminv(hr) - norminv(far);


    p_act_given_hit = hr .* .5 ./ phit;
    p_null_given_no = (1 - far) .* .5 ./ (1 - phit);

    % Probability of correct classification (signif vs. not, yes vs no) given
    % true activation
    % i.e., p that sig. matches true state
    % ------------------------------------------------------------------------
    % p(correct classification) = p(yes|activ)p(activ) + p(no|inactiv)p(inactiv)
    % if p(activ) = p(inactiv), p(correct) = .5HR + .5CR
    pcorrect = .5 * (hr + (1 - far));

    % Probability of correct true state  given
    % observed classification (signif vs. not, yes vs no)
    % i.e., p that true state is what I said it was
    % ------------------------------------------------------------------------
    pcorrectclass = .5 * (p_act_given_hit + p_null_given_no);


    rocfig = create_figure('Sim Output', 1, nsubplots); plot(far,hr); xlabel('False Positive Rate'); ylabel('True Positive Rate');
    title(['ROC plot for ' names{wh_col}]);

    % ========================================================
    %% correct classification rate plot
    % ========================================================
    subplot(1,nsubplots,2); plot(pthr,pcorrect); xlabel('p-value cutoff'); ylabel('p(correct classification)');

    % maybe i did something wrong here?  seems like they should be the same...
    hold on; plot(pthr,pcorrectclass,'r'); han = legend({'p(sig/not given true)' 'p(true given sig/not)'});
    set(han,'Location','SouthWest');
    
    title('Correct classification');

    % ========================================================
    %% put points for optimal classification on ROC
    % ========================================================

    % stats at optimal threshold for correct classification
    % ---------------------------
    [best_classrate,wh] = max(pcorrect);
    best_classrate = best_classrate(1);
    wh = wh(1);

    best_pthr = pthr(wh);

    figure(rocfig); subplot(1,nsubplots,1);

    plot(far(wh), hr(wh), 'o','Color',[.2 .2 .2],'MarkerFaceColor','y','MarkerSize',8,'LineWidth',2);
    text(far(wh) + .05, hr(wh), sprintf('p threshold: %3.4f',best_pthr),'FontSize',16);
    text(far(wh) + .05, hr(wh) - .05, sprintf('p(correct): %3.0f%%',100*best_classrate),'FontSize',16);

    % stats at .05
    % ---------------------------
    wh = find(pthr == .05);
    plot(far(wh), hr(wh), 'o','Color',[.2 .2 .2],'MarkerFaceColor','r','MarkerSize',8,'LineWidth',2);
    text(far(wh) + .05, hr(wh), sprintf('p threshold: %3.4f',.05),'FontSize',16);
    text(far(wh) + .05, hr(wh) - .05, sprintf('p(correct): %3.0f%%',100*pcorrect(wh)),'FontSize',16);
    text(far(wh) + .05, hr(wh) - .1, sprintf('Obs. FPR: %3.4f',far(wh)),'FontSize',16);

    % ========================================================
    %% plot obs. fpr (alpha) vs. nominal alpha
    % ========================================================
    subplot(1,nsubplots,3);

    plot(pthr,far,'b')
    hold on; plot([0 1],[0 1],'k--');

    title('Bias: Observed alpha vs. nominal alpha');
    xlabel('Nominal alpha'); ylabel('Observed alpha');


    set(gcf,'Position',[144        1029        1089         441]);
    scn_export_papersetup(400);
    saveas(gcf,['rocplot_' names{wh_col}],'png');
    disp('Saved rocplot.png');


end