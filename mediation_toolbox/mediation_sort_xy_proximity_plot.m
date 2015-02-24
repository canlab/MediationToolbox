function out = mediation_sort_xy_proximity_plot(clpos_data, xdata, ydata)
    % out = mediation_sort_xy_proximity_plot(clpos_data, xdata, ydata)
    %
    % % load mediation_SETUP
    % % load clusters_005_01_05_prune_graymask_withdata
    % %     xdata = SETUP.data.X;
    % %     ydata = SETUP.data.Y;

    kval = .5;
    N = length(clpos_data);
    nregions = length(clpos_data{1});
    
    tthr = tinv(1 - (.05) ./ 2, N - 1);

    for i = 1:length(clpos_data), mdata{i} = cat(2, clpos_data{i}.timeseries); end
    
        % Intra-brain relationships, standardized betas
        % Each connection with each other one, controlling for all other
        % regions
    % --------------------------------------------------------
    disp('Getting intra-brain connectivity');
    out = matrix_direct_effects_ridge(mdata);
    
        
    % Each brain time series with x and y controlling for the other -- for non-ridge plot
    % --------------------------------------------------------
    disp('Getting XY connectivity');
    
    for i = 1:N
        X{i} = [xdata{i} ydata{i}];  % design mtx, no intercept
    end
    
    for i = 1:nregions
        % For each region, calculate analysis of X and Y predictors:
        for j = 1:N, Y{j} = mdata{j}(:, i); end
        
        stats = glmfit_multilevel(Y, X, ones(N, 1), 'beta_names', {'Intercept' 'X' 'Y'}, 'weighted', 'noverbose');
        
        bxy(i, :) = stats.beta(2:3);  % could use b_star ***?
        txy(i, :) = stats.t(2:3); 
        pxy(i, :) = stats.p(2:3);
        sigxy(i, :) = pxy(i, :) < .05;
    end
    
    %[bxy, txy, sigxy] = get_xy_effects(xdata, ydata, mdata, tthr);
    
    % PLOT: Each brain region, x controlling for Y and vice versa
    % --------------------------------------------------------
    % betas, etc. for each of x and y controlling for the other; average over
    % subjects (t based on ste over subjects); one cell per region
    %bxy = cat(1, bxy{:});
    % Below is new x-coords, 0 for all X, 1 for all Y
    xycomb = txy(:, 2) .^ 2 ./ (txy(:, 1) .^ 2 + txy(:, 2) .^ 2);

    
    whx = find(sigxy(:, 1) & ~sigxy(:, 2) & bxy(:, 1) > 0);
montage_clusters([], clpos_data{1}(whx), {[0 0 1]});
saveas(gcf,'Xonly_pos', 'png');
whx = find(sigxy(:, 1) & ~sigxy(:, 2) & bxy(:, 1) < 0);
montage_clusters([], clpos_data{1}(whx), {[0 0 1]});
saveas(gcf,'Xonly_neg', 'png');
why = find(sigxy(:, 2) & ~sigxy(:, 1) & bxy(:, 2) > 0);
montage_clusters([], clpos_data{1}(why), {[1 0 0]});
saveas(gcf,'Yonly_pos', 'png');
why = find(sigxy(:, 2) & ~sigxy(:, 1) & bxy(:, 2) < 0);
why
montage_clusters([], clpos_data{1}(why), {[1 0 0]});
saveas(gcf,'Yonly_neg', 'png');

    
    
    whx = find(sigxy(:, 1) & ~sigxy(:, 2));  % x not y
    why = find(sigxy(:, 2) & ~sigxy(:, 1));  % y not x
    whb = find(sigxy(:, 2) & sigxy(:, 1));  % y and x
    
    whn = find(~sigxy(:, 1) & ~sigxy(:, 2)); % neither
    
    montage_clusters([], clpos_data{1}(whx), {[0 0 1] [1 0 0] [1 1 0]}, clpos_data{1}(why), clpos_data{1}(whb));
    
    scn_export_papersetup(800)
    saveas(gcf, 'Sort_XY_stim_resp_only', 'png');


    montage_clusters([], clpos_data{1}(whn), {[1 0 1]});
    scn_export_papersetup(800)
    saveas(gcf, 'Sort_XY_cannottell', 'png');
    
    xystrength = zeros(nregions, 1);
    xystrength(whx) = txy(whx, 1); %
    xystrength(why) = txy(why, 2);
    
    xystrength(xystrength == 0) = sign(txy((xystrength == 0), 1) .* txy((xystrength == 0), 2)) .* sum(txy(xystrength == 0, :) .^ 2, 2) .^ .5;
    
% %     ( bxy(:, 1) .* bxy(:, 2) );
% %     xystrength = sign(xystrength) .* (abs(xystrength) .^ .5);
    
    % For each subject, regress mdata for each region on xdata and ydata
    % (All vars are z-scored first) scale(xx) 
    % save betay, betax for each region, and also beta^2y(sum(beta^x,beta^2y)

    xycoords(:, 1) = xycomb;
    xycoords(:, 2) = xystrength;
    create_figure('X-Y Prox plot', 1, 1); set(gca, 'XLim', [-.1 1.1]);
    
    % Position 0: X
    plot(0, 0, 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 12)

    % Position 8: Y
    plot(1, 0, 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 12)
    
    % plot points
    for i = 1:size(xycoords, 1)
        mfcolor = xycoords(i, 1) .* [1 0 0] + (1 - xycoords(i, 1)) .* [0 0 1]; 
        plot(xycoords(i, 1), xycoords(i, 2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', mfcolor);
    end
    
    set(gca, 'FontSize', 24)
    xlabel('Relationship with stim. vs. response');
    set(gca, 'XTick', [0 1], 'XTickLabel', {'Stim only' 'Resp only'});
    ylabel('Mediation Strength');
    %set(gca, 'YTick', [-1 min(xystrength) 0 max(xystrength) 1], 'YTickLabel', {'Max Poss Neg' 'Negative' 'Zero' 'Positive' 'Max Poss Pos'});
    
    sigmat = out.fdrsig;
    %handles = nmdsfig_tools('drawlines', xycoords, sigmat, [0 0 0;1 0 0],{'-',':'}, .01);

    % Put names on
    if isfield(clpos_data{1}, 'shorttitle')
        for i = 1:length(clpos_data{1})
            if xycoords(i, 1) ~= 0 && xycoords(i, 2) ~= 0
                text(xycoords(i, 1) + .1, xycoords(i, 2), clpos_data{1}(i).shorttitle, 'Color', 'b', 'FontSize', 16, 'FontWeight', 'b');
            end
        end
    end

    lineh = plot_horizontal_line(0); set(lineh, 'Color', 'g', 'LineStyle', '--');
    drawnow
    
    scn_export_papersetup(500); saveas(gcf,'Mediation_sort_XY_plot', 'png');
    % --------------------------------------------------------
    
    
    



    % STIMULUS-RELATED, standardized betas
    % --------------------------------------------------------
    disp('Getting X connectivity');
    [out.bx, out.tx, out.sigx] = get_xeffect_ridge(xdata, mdata, kval, tthr);

    % RESPONSE-RELATED, standardized betas
    % --------------------------------------------------------
    disp('Getting Y connectivity');
    [out.by, out.ty, out.sigy] = get_yeffect_ridge(ydata, mdata, kval, tthr);



    % Define x-values and region indices
    % ----------------------------------------------------
    [wh_regions, grpcoords, xpositions] = get_region_indices_by_xcoords(out);

    % Define y-values and region indices
    % ----------------------------------------------------
    dosort = 0;  % if sorting, will return index values for axis plots rather than actual values
    [yvals, nrows, offset] = get_y_values(wh_regions, dosort, out, xpositions);



    
    
    
    %% Figure
    create_figure('Parallel Sort Plot', 1, 1); set(gca, 'XLim', [-.5 8.5]);

    set(gca,'XTick', [0 xpositions 8], 'XTickLabels', {'X' 'X only' 'X_2' 'X and Y' 'XY_2' 'Y_2' 'Y only' 'Y'});


    % X coords in plot define relation to X vs Y
    % Position 0: X
    % Position 8: Y
    % Position 1: X but not Y
    % Position 4: X and Y
    % Position 7: Y but not X

    % Position 0: X
    plot(0, 0, 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 12)

    % Position 8: Y
    plot(8, 0, 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 12)

    colors = {[0 0 .8] [.3 0 .8] [.8 0 .8] [.8 0 .5] [.8 0 .3] [.8 0 0]};

    for i = 1:length(xpositions)
        plot_layer(wh_regions{i}, xpositions(i), yvals{i}, colors{i}, out, offset);

        if ~isempty(wh_regions{i}), grpcoords(wh_regions{i}, 2) = yvals{i}; end

    end


    whzero = find(grpcoords(:, 1) == 0);
    fprintf('Remaining regions (omitted): %3.0f\n', length(whzero))

    sigmat = out.fdrsig;
    grpcoords_nonzero = grpcoords;
    grpcoords_nonzero(whzero, :) = [];
    sigmat(whzero, :) = [];
    sigmat(:, whzero) = [];

    handles = nmdsfig_tools('drawlines',grpcoords_nonzero, sigmat, [0 0 0;1 0 0],{'-',':'}, .01);
    %set([handles.hhp handles.hhn], 'LineWidth', 2);

    % Put names on
    if isfield(clpos_data{1}, 'shorttitle')
        for i = 1:length(clpos_data{1})
            if grpcoords(i, 1) ~= 0 && grpcoords(i, 2) ~= 0
                text(grpcoords(i, 1) + .1, grpcoords(i, 2), clpos_data{1}(i).shorttitle, 'Color', 'b', 'FontSize', 16, 'FontWeight', 'b');
            end
        end
    end

    lineh = plot_horizontal_line(0); set(lineh, 'Color', 'g', 'LineStyle', '--');

    yl = get(gca, 'YLim');
    text(0, yl(2), 'Positive net effect', 'Color', 'g', 'FontSize', 16, 'FontWeight', 'b');
    text(0, yl(1) + .1, 'Negative net effect', 'Color', 'g', 'FontSize', 16, 'FontWeight', 'b');

    

    plot_slice_figure();


    create_figure('XYplot'); plot(out.bx, out.by, 'ko', 'MarkerFaceColor', 'r');
    xlabel('Relationship with temperature'); ylabel('Relationship with pain report');
    plot_vertical_line(0); plot_horizontal_line(0);


    % inline

    function plot_slice_figure

        %% Figure

        
        % Define x-values and region indices
        % ----------------------------------------------------
        [wh_regions, grpcoords, xpositions] = get_region_indices_by_xcoords(out);

        % Define y-values and region indices
        % ----------------------------------------------------
        dosort = 1;  % if sorting, will return index values for axis plots rather than actual values
        [yvals, nrows, offset] = get_y_values(wh_regions, dosort, out, xpositions);
        

        create_figure('Parallel Sort Plot Slices', 1, 1); set(gca, 'XLim', [-.5 8.5]);

        set(gca,'XTick', [0 xpositions 8], 'XTickLabels', {'X' 'X only' 'X_2' 'X and Y' 'XY_2' 'Y_2' 'Y only' 'Y'});


        % Position 0: X
        plot(0, offset, 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 12)

        % Position 8: Y
        plot(8, offset, 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 12)

        for i = 1:length(xpositions)
            plot_layer(wh_regions{i}, xpositions(i), yvals{i}, [1 0 0], out, offset);

            if ~isempty(wh_regions{i}), grpcoords(wh_regions{i}, 2) = yvals{i}; end

        end

        %%

        whzero = find(grpcoords(:, 1) == 0);
        fprintf('Remaining regions (omitted): %3.0f\n', length(whzero))

        sigmat = out.fdrsig;
        grpcoords(whzero, :) = [];
        sigmat(whzero, :) = [];
        sigmat(:, whzero) = [];

        handles = nmdsfig_tools('drawlines',grpcoords, sigmat, [0 0 0;0 0 1],{'-',':'}, .01);


        %%
        axh = gca;
        bigaxpos = get(axh, 'Position');
        bigxrange = bigaxpos(3); % range([bigaxpos(1) bigaxpos(1) + bigaxpos(3)]);
        bigyrange = bigaxpos(4); %range([bigaxpos(2) bigaxpos(2) + bigaxpos(4)]);
        xl = get(axh, 'XLim');
        yl = get(axh, 'YLim');
        xrange = range(xl);
        yrange = range(yl);

        for i = 1:length(xpositions)

            % slices
            for j = 1:length(wh_regions{i})

                % no, need to get location on axis for position....
                pos = [xpositions(i) - .4 yvals{i}(j) - .4 .8 .8]; % relative to axis
                % scale relative to zero/1 on big axis
                pos(1) = (pos(1) - xl(1)) ./ xrange;
                pos(3) = pos(3) ./ xrange;
                pos(2) = (pos(2) - yl(1)) ./ yrange;
                pos(4) = pos(4) ./ yrange;

                % adjust for fact that big axis doesn't cover whole area from 0 - 1
                pos(1) = pos(1) .* bigxrange + bigaxpos(1);
                pos(3) = pos(3) .* bigxrange;
                pos(2) = pos(2) .* bigyrange + bigaxpos(2);
                pos(4) = pos(4) .* bigyrange;

                han(j) = axes('Position', pos);
                axes(han(j))
                montage_clusters_maxslice([], clpos_data{1}(wh_regions{i}(j)), colors(i));

            end
        end

    end
    %%

end % main function

function plot_layer(wh, xaxispos, yvals, color, out, offset)
    % wh, which objects
    % xaxispos, x-axis position.  yvals, y-axis position for each object

    if nargin < 6, offset = 0; end
    
    nobj = length(wh);

    % dots
    plot(xaxispos * ones(1, nobj), yvals, 'ko', 'MarkerFaceColor', color, 'MarkerSize', 12)

    % xlines
    whx = find(out.sigx(wh));
    xsigns = sign(out.bx(wh));
    
    nobjx = length(whx);
    if nobjx > 0
        for j = 1:nobjx
            switch xsigns(j)
                case 1, lcolor = 'k';
                case -1, lcolor = 'r--';
                otherwise, error('Uh-OH!');
            end
                    
            plot([0 xaxispos]', [offset yvals(whx(j))]', lcolor, 'LineWidth', 1);
        end
        %plot(repmat([0 xaxispos], nobjx, 1)', [repmat(offset, nobjx, 1) yvals(whx)]', 'k');
    end

    % ylines
    why = find(out.sigy(wh));
    ysigns = sign(out.by(wh));
    
    nobjy = length(why);
    if nobjy > 0
        for j = 1:nobjy
            switch ysigns(j)
                case 1, lcolor = 'k';
                case -1, lcolor = 'r--';
                otherwise, error('Uh-OH!');
            end

            plot([8 xaxispos]', [offset yvals(why(j))]', lcolor, 'LineWidth', 1);
        end
        
        %plot(repmat([8 xaxispos], nobjy, 1)', [repmat(offset, nobjy, 1) yvals(why)]', 'k');
    end
    
    plot(xaxispos * ones(1, nobj), yvals, 'ko', 'MarkerFaceColor', color, 'MarkerSize', 12)

end


function [bx, tx, sigx] = get_xeffect_ridge(xdata, mdata, kval, tthr)
    N = length(xdata);
    v = size(mdata{1}, 2);

    b = zeros(N, v);
    for i = 1:N

        for j = 1:v

            y = mdata{i}(:, j);  % outcome is brain data in mediation

            X = [xdata{i} mdata{i}];  % predictors are X and other brain data in mediation
            X(:, j + 1) = [];         % remove predictor of interest
            X(:,end+1) = 1;

            [nanvec, X] = nanremove(X'); X = X';

            y = scale(y);
            X = scale(X);
            X(:,end) = 1;

            bsubi =  ridge(y, X, kval);
            bsubi = naninsert(nanvec, bsubi)';

            b(i, j) = bsubi(1);     % take relation to X only
        end

    end

    bx = nanmean(b)';
    tx = bx ./ ste(b)';
    sigx = abs(tx) > tthr;

end



function [by, ty, sigy] = get_yeffect_ridge(ydata, mdata, kval, tthr)
    N = length(ydata);
    v = size(mdata{1}, 2);

    b = zeros(N, v + 1);

    for i = 1:N
        y = ydata{i};  % outcome is x or y var in mediation
        X = mdata{i};  % predictors are brain data in mediation
        X(:,end+1) = 1;

        [nanvec, X] = nanremove(X'); X = X';

        y = scale(y);
        X = scale(X);
        X(:,end) = 1;

        bsubi =  ridge(y, X, kval);
        b(i,:) = naninsert(nanvec, bsubi)';

    end

    b(:, end) = []; % remove intercept parameters

    by = nanmean(b)';
    ty = by ./ ste(b)';
    sigy = abs(ty) > tthr;

end


function [bxy, txy, sigxy] = get_xy_effects(xdata, ydata, mdata, tthr)
% betas, etc. for each of x and y controlling for the other; average over
% subjects (t based on ste over subjects); one cell per region

N = length(ydata);
v = size(mdata{1}, 2); % num regions

b = cell(1, v); % betas

for i = 1:N
    y = ydata{i};  % outcome is x or y var in mediation
    x = xdata{i};  % predictors are brain data in mediation

    for r = 1:v
        m = mdata{i}(:, r); % brain region to predict

        [nanvec, xx, yy, m] = nanremove(x, y, m);

        yy = scale(yy);
        xx = scale(xx);
        m = scale(m);

        betas = pinv([yy xx ones(size(yy,1), 1)]) * m;  % x and y predicting m

        b{r}(i,:) = betas(1:2)';

    end

end

for r = 1:v
    bxy{r} = nanmean(b{r});
    txy{r} = bxy{r} ./ ste(bxy{r});
    sigxy{r} = abs(bxy{r}) > tthr;

end

end


function [wh_regions, grpcoords, xpositions] = get_region_indices_by_xcoords(out) %, xpositions)
    grpcoords = zeros(length(out.sigx), 2);

    % X coords in plot define relation to X vs Y
    % Position 0: X
    % Position 8: Y
    % Position 1: X but not Y
    % Position 4: X and Y
    % Position 7: Y but not X

    xpositions = [1 2.5 4 5 6 7];

    % Define list of regions for each x position on graph
    % ----------------------------------------------------
    for i = [1 4 7 2.5 5 6]  %1:length(xpositions)
        
        myxpos = find(xpositions == i);
        
        switch i %xpositions(i)
            case 1
                % Position 1: X but not Y
                wh_regions{myxpos} = find(out.sigx & ~out.sigy);
            case 4
                % Position 4: X and Y
                wh_regions{myxpos} = find(out.sigx & out.sigy);
            case 7
                % Position 7: Y but not X
                wh_regions{myxpos} = find(~out.sigx & out.sigy);

            case 2.5
                % Position 2.5: Connected to Position 1 vars (2nd layer)
                % (X-connected but not y-connected), and not yet on plot
                pos1vars = find(grpcoords(:, 1) == 1);
                connectpos1 = any(out.fdrsig(pos1vars, :))';  % sig is rows predicting columns, so this is pos1vars to connected vars localized

                % eliminate vars already on plot, if any
                connectpos1(grpcoords(:, 1) ~= 0) = 0;

                wh_regions{myxpos} = find(connectpos1);

            case {4.5, 5} % XY-2
                % Connected to XY mediators...
                pos1vars = find(grpcoords(:, 1) == 4);
                connectpos1 = any(out.fdrsig(pos1vars, :))';  % sig is rows predicting columns, so this is pos1vars to connected vars localized

                % eliminate vars already on plot, if any
                connectpos1(grpcoords(:, 1) ~= 0) = 0;

                wh_regions{myxpos} = find(connectpos1);

            case 6
                % Position 6: Connected to Position 1 vars (2nd layer)
                % (X-connected but not y-connected), and not yet on plot
                pos1vars = find(grpcoords(:, 1) == 7);
                connectpos1 = any(out.fdrsig(pos1vars, :))';  % sig is rows predicting columns, so this is pos1vars to connected vars localized

                % eliminate vars already on plot, if any
                connectpos1(grpcoords(:, 1) ~= 0) = 0;
                wh_regions{myxpos} = find(connectpos1);

            case 4.5  % XY-2
% %                 % Position 4.5: Connected to Position 4 vars (2nd layer) (X-connected but not y-connected), and not yet on plot
% %                 pos1vars = find(grpcoords(:, 1) == 4);
% %                 connectpos1 = any(out.fdrsig(pos1vars, :))';  % sig is rows predicting columns, so this is pos1vars to connected vars localized
% % 
% %                 % eliminate vars already on plot, if any
% %                 connectpos1(grpcoords(:, 1) ~= 0) = 0;
% %                 wh_regions{i} = find(connectpos1);

        end

        grpcoords(wh_regions{myxpos}, 1) = i; %xpositions(i);
    end


end

function [yvals, nrows, offset] = get_y_values(wh_regions, dosort, out, xpositions)
    % Define y-values
    % ----------------------------------------------------
    clear yvals

    for i = 1:length(wh_regions), nreg(i) = length(wh_regions{i}); end
    nrows = max(nreg);
    if dosort, offset = mean(1:nrows); else, offset = 0; end

    for i = 1:length(xpositions)
        switch xpositions(i)
            case 1
                % Position 1: X but not Y
                yvals{i} = out.bx(wh_regions{i});

            case 7
                % Position 1: X but not Y
                yvals{i} = out.by(wh_regions{i});
                
            case {4, 2.5, 6, 4.5, 5}
                % Position 4: X and Y
                prodvals = out.bx(wh_regions{i}) .* out.by(wh_regions{i});
                yvals{i} = sign(prodvals) .* sqrt(abs(prodvals));
        end

        if dosort
            yvals{i} = rankdata(yvals{i});
            if ~isempty(yvals{i}), yvals{i} = yvals{i} - mean(yvals{i}) + offset; end

        end
    end

end


  
% %
% % %%
% % overlay = '/Users/tor/Documents/Tor_Documents/CurrentExperiments/Mediation_Grant/NSF_pain_map_study/nsf_mean_spm5wspgr.img'
% %
% % montage_clusters(overlay, clpos_data{1}(wh), clpos_data{1}(whr), {'r' 'b'});
% %
% % create_figure('scatter'); plot(bpop, bpopr, 'ko', 'MarkerFaceColor', 'b')
% % plot(bpop(wh), bpopr(wh), 'ro', 'MarkerSize', 12, 'LineWidth', 2);
% % plot(bpop(whr), bpopr(whr), 'bo', 'MarkerSize', 16 , 'LineWidth', 2);
% %
% % plot(bpop(whn), bpopr(whn), 'o', 'Color', [1 .3 .3], 'MarkerSize', 12, 'LineWidth', 2);
% % plot(bpop(whrn), bpopr(whrn), 'co', 'MarkerSize', 16 , 'LineWidth', 2);
% % xlabel('Stimulus relation (standardized beta)');
% % ylabel('Response relation (standardized beta)');
% %
% % input('Press return...');
% % montage_clusters([], clpos_data{1}(whn), clpos_data{1}(whrn), {'m' 'c'});
% %
% % %%
% %
% % classes = zeros(out.k, 1);
% % classes(wh) = 1;  % pos related to X
% % classes(whr) = 2;  % pos related to Y
% % classes(whn) = 3;  % neg related to X
% % classes(whrn) = 4;  % neg related to Y
% %
% %
% % colors = {'ko' 'ro' 'bo' 'mo' 'co'};
% % nmdsfig(out.GroupSpace,'classes',classes + 1, 'names',[],'sig',out.fdrsig, 'legend',{'Pos' 'Neg'},'sizescale',[4 16], 'colors', colors);
% %
% % %



% %     % Position 1: X but not Y
% %     wh = find(out.sigx & ~out.sigy);
% %     xaxispos = 1;
% %     yvals = out.bx(wh);
% %     plot_layer(wh, xaxispos, yvals, [0 0 .8], out);
% %     grpcoords(wh, 1) = xaxispos; grpcoords(wh, 2) = yvals;
% % 
% %     % Position 4: X and Y
% %     wh = find(out.sigx & out.sigy);
% %     xaxispos = 4;
% %     yvals = mean([out.bx(wh) out.by(wh)], 2);
% %     plot_layer(wh, xaxispos, yvals, [.8 0 .8], out);
% %     grpcoords(wh, 1) = xaxispos; grpcoords(wh, 2) = yvals;
% % 
% %     % Position 7: Y but not X
% %     wh = find(~out.sigx & out.sigy);
% %     xaxispos = 7;
% %     yvals = mean([out.bx(wh) out.by(wh)], 2);
% %     plot_layer(wh, xaxispos, yvals, [.8 0 0], out);
% %     grpcoords(wh, 1) = xaxispos; grpcoords(wh, 2) = yvals;
% % 
% %     % Position 2.5: Connected to Position 1 vars (2nd layer) (X-connected but not y-connected), and not yet on plot
% %     pos1vars = find(grpcoords(:, 1) == 1);
% %     connectpos1 = any(out.fdrsig(pos1vars, :))';  % sig is rows predicting columns, so this is pos1vars to connected vars localized
% % 
% %     % eliminate vars already on plot, if any
% %     connectpos1(grpcoords(:, 1) ~= 0) = 0;
% %     wh = find(connectpos1);
% % 
% %     xaxispos = 2.5;
% %     yvals = mean([out.bx(wh) out.by(wh)], 2);
% %     plot_layer(wh, xaxispos, yvals, [.3 0 .3], out);
% %     grpcoords(wh, 1) = xaxispos; grpcoords(wh, 2) = yvals;
% % 
% %     % Position 6: Connected to Position 1 vars (2nd layer) (X-connected but not y-connected), and not yet on plot
% %     pos1vars = find(grpcoords(:, 1) == 7);
% %     connectpos1 = any(out.fdrsig(pos1vars, :))';  % sig is rows predicting columns, so this is pos1vars to connected vars localized
% % 
% %     % eliminate vars already on plot, if any
% %     connectpos1(grpcoords(:, 1) ~= 0) = 0;
% %     wh = find(connectpos1);
% % 
% %     xaxispos = 6;
% %     yvals = mean([out.bx(wh) out.by(wh)], 2);
% %     plot_layer(wh, xaxispos, yvals, [.3 0 .3], out);
% %     grpcoords(wh, 1) = xaxispos; grpcoords(wh, 2) = yvals;
% % 
% % 
% %     % Position 4.5: Connected to Position 4 vars (2nd layer) (X-connected but not y-connected), and not yet on plot
% %     pos1vars = find(grpcoords(:, 1) == 4);
% %     connectpos1 = any(out.fdrsig(pos1vars, :))';  % sig is rows predicting columns, so this is pos1vars to connected vars localized
% % 
% %     % eliminate vars already on plot, if any
% %     connectpos1(grpcoords(:, 1) ~= 0) = 0;
% %     wh = find(connectpos1);
% % 
% %     xaxispos = 4.5;
% %     yvals = mean([out.bx(wh) out.by(wh)], 2);
% %     plot_layer(wh, xaxispos, yvals, [.8 0 .8], out);
% %     grpcoords(wh, 1) = xaxispos; grpcoords(wh, 2) = yvals;


    % NEXT:: Names
    % Put pictures below graph: Axes of each cluster
    %
    % Y positions?

    % X coords in plot define relation to X vs Y
    % Position 0: X
    % Position 8: Y
    % Position 1: X but not Y
    % Position 4: X and Y
    % Position 7: Y but not X
% % %     subplot(2, 1, 2); axis off
% % % 
% % %     subplot(2, 1, 1); set(gca, 'Position', [0.1300    0.6838    0.7750    0.2712]);
% % %     %%
% % % 
% % %     yposvalues = .55:-.1:.05; %[.45 .35 .25 .15 .05];
% % %     xposvalues = [.23:.1:.73];
% % %     indx =1 ;
% % % 
% % %     for x = 1:length(xposvalues)
% % %         for y = 1:length(yposvalues)
% % % 
% % %             ah2(x, y) = axes('position', [xposvalues(x) yposvalues(y) .09 .09]);
% % %             axis off
% % %         end
% % %     end
% % % 
% % %     for x = 1:length(xpositions)
% % %         wh = find(grpcoords(:, 1) == xpositions(x));
% % %         % re-sort, highest to lowest y
% % %         [yvals, ii] = sort(grpcoords(wh, 2), 1, 'descend');
% % %         wh = wh(ii);
% % %         indx = 1;p
% % %         for i = 1:length(wh)
% % % 
% % %             if indx <= size(ah2, 2)
% % %                 axes(ah2(x, indx));
% % %                 montage_clusters_maxslice([], clpos_data{1}(wh(i)), {[.4 0 .4]});
% % %             else
% % %                 disp('Cannot display cluster: Too many to plot at this x-value.');
% % %             end
% % %             indx = indx + 1;
% % %         end
% % % 
% % %     end

