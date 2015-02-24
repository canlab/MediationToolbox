function [cl_apos_bpos, cl_aneg_bneg, cl_apos_bneg, cl_aneg_bpos] = mediation_brain_results_a_b_overlap(varargin)
% [cl_apos_bpos, cl_aneg_bneg, cl_apos_bneg, cl_aneg_bpos] = mediation_brain_results_a_b_overlap(varargin)
%
% Compute images and clusters for the overlap of a and b path images
% 
% Optional inputs:
% case 'p', p = varargin{i+1};
% case 'k', k = varargin{i+1};
% case 'overlay', overlay = varargin{i+1};
% case 'mask', mask = varargin{i+1};
% case 'fdr', dofdr = 1;        % use across-image FDR correction
% case 'slices', doslices = 1;  % make montages of slices
% case 'save', dosave = 1       % save images and clusters
%
% Tor Wager, Jan 09
% NOTE: DOCUMENTATION AND BUG TESTING NOT COMPLETE
%
% Examples:
% [cl_apos_bpos, cl_aneg_bneg] = mediation_brain_results_a_b_overlap('fdr', 'k', 3, 'overlay', EXPT.overlay, 'save');

% Programmer's notes:
% 9/25/09 : Tor : fixed incompatibility with single-level model results in
% terms of data extraction

overlay = [];
mask = fullfile(pwd, 'mask.img');

k = [];
p = [];
dofdr = 0;
doslices = 0;
dosave = 0;

fname = fullfile(pwd, 'mediation_SETUP.mat');
if exist(fname, 'file'), load(fname); else error('No mediation_SETUP.mat. Go to valid mediation directory first.'); end

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case 'p', p = varargin{i+1};
            case 'k', k = varargin{i+1};
                
            case 'overlay', overlay = varargin{i+1}; varargin{i+1} = [];
            case 'mask', mask = varargin{i+1}; varargin{i+1} = [];
    
            case 'fdr', dofdr = 1;
            case 'slices', doslices = 1;
            case 'save', dosave = 1;

            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if dofdr
    if ~isfield(SETUP, 'fdr_p_thresh')
        SETUP = mediation_brain_corrected_threshold('fdr', 'mask', mask);
    end

    p = SETUP.fdr_p_thresh;
end

% get output string
fmts = ['%0.' num2str(ceil(log10(1/p))) 'f'];
mystr = sprintf(fmts, p);
wh = find(mystr == '.');
if ~isempty(wh) && wh(1) ~= 1, mystr(wh(1)-1:wh(1)) = []; end % remove zero_dot

% output image names
posimagename = ['intersect_apos_bpos_' mystr '.img'];
negimagename = ['intersect_aneg_bneg_' mystr '.img'];

% input images
imgs = char('X-M_effect.img', 'X-M_pvals.img', 'M-Y_effect.img', 'M-Y_pvals.img');

volInfo = iimg_read_img(mask, 2);
dat = iimg_get_data(volInfo, imgs);  % rows are images

% compute it
pdat = all(dat([2 4], :) < p);
signdat = sign(dat([1 3], :));

% valid_vox = double(pdat & all(signdat > 0))';
% [cl_apos_bpos, datpospos] = iimg_indx2clusters(valid_vox, volInfo, .5, k);
% 
% valid_vox = double(pdat & all(signdat > 0))';
% [cl_aneg_bneg, datnegneg] = iimg_indx2clusters(valid_vox, volInfo, .5, k);


str = ['i1>0 & i2< ' num2str(p) ' & i3>0 & i4<' num2str(p) ' & ~isnan(i1) & ~isnan(i2) & i5 > .5'];
cl_apos_bpos = mask_intersection2(k, posimagename, char(imgs, mask), str);

str = ['i1<0 & i2< ' num2str(p) ' & i3<0 & i4<' num2str(p) ' & ~isnan(i1) & ~isnan(i2) & i5 > .5'];
cl_aneg_bneg = mask_intersection2(k, negimagename, char(imgs, mask), str);


cluster_orthviews(cl_apos_bpos, {[1 0 0]}, 'overlay', overlay);
cluster_orthviews(cl_aneg_bneg, {[0 0 1]}, 'add');


if doslices && ~isempty(cl_apos_bpos) && ~isempty(cl_aneg_bneg)
    cluster_orthviews_showcenters([cl_apos_bpos cl_aneg_bneg], 'axial', overlay, 0, 1);
    h = findobj(gcf,'Type', 'text','Color', 'g');
    delete(h)
    h = findobj(gcf,'Type', 'text','Color', 'k');
    delete(h)
    enlarge_axes(gcf, 1.15)
    enlarge_axes(gcf, 1, 1.05)
    if dosave
        saveas(gcf,['intersect_ ' mystr '_a_and_b_axial'], 'png')
    end

    cluster_orthviews_showcenters([cl_apos_bpos cl_aneg_bneg], 'sagittal', overlay, 0, 1);
    h = findobj(gcf,'Type', 'text','Color', 'g');
    delete(h)
    h = findobj(gcf,'Type', 'text','Color', 'k');
    delete(h)
    enlarge_axes(gcf, 1.15)
    enlarge_axes(gcf, 1, 1.05)
    if dosave
        saveas(gcf,['intersect_' mystr '_a_and_b_sagittal'], 'png')
    end
    
    cluster_orthviews_showcenters([cl_apos_bpos cl_aneg_bneg], 'coronal', overlay, 0, 1);
    h = findobj(gcf,'Type', 'text','Color', 'g');
    delete(h)
    h = findobj(gcf,'Type', 'text','Color', 'k');
    delete(h)
    enlarge_axes(gcf, 1.15)
    enlarge_axes(gcf, 1, 1.05)
    if dosave
        saveas(gcf,['intersect_' mystr '_a_and_b_coronal'], 'png')
    end
end
%% inconsistent

str = ['i1>0 & i2< ' num2str(p) ' & i3<0 & i4<' num2str(p) ' & ~isnan(i1) & ~isnan(i2) & i5 > .5'];
cl_apos_bneg = mask_intersection2(k, 'intersect_apos_bneg.img', char(imgs, mask), str);

str = ['i1<0 & i2< ' num2str(p) ' & i3>0 & i4<' num2str(p) ' & ~isnan(i1) & ~isnan(i2) & i5 > .5'];
cl_aneg_bpos = mask_intersection2(k, 'intersect_aneg_bpos.img', char(imgs, mask), str);

cluster_orthviews(cl_apos_bneg, {[0 1 0]}, 'add');
cluster_orthviews(cl_aneg_bpos, {[.5 1 .5]}, 'add');

if dosave
    save intersect_a_and_b_clusters cl_*
end

% Extract

% fix for diffs in how single-level and multi-level data are stored
if ~isfield(SETUP, 'data')
    SETUP.data.M = {SETUP.M};
end

if ~isempty(cl_apos_bpos)
    cl_apos_bpos = extract_raw_data([], cl_apos_bpos, 'extract_from', SETUP.data.M, 'noraw');
end

if ~isempty(cl_aneg_bneg)
    cl_aneg_bneg = extract_raw_data([], cl_aneg_bneg, 'extract_from', SETUP.data.M, 'noraw');
end

if ~isempty(cl_apos_bneg)
    cl_apos_bneg = extract_raw_data([], cl_apos_bneg, 'extract_from', SETUP.data.M, 'noraw');
end

if ~isempty(cl_aneg_bpos)
    cl_aneg_bpos = extract_raw_data([], cl_aneg_bpos, 'extract_from', SETUP.data.M, 'noraw');
end

if ~isfield(cl_apos_bpos, 'all_data')
    % we have not extracted data successfully; skip by doing nothing
else
    for i = 1:length(cl_apos_bpos)
        k = size(cl_apos_bpos(i).all_data, 2);

        cl_apos_bpos(i).timeseries = cell(1, k);

        for j = 1:k
            cl_apos_bpos(i).timeseries{j} = cl_apos_bpos(i).all_data(:, j);
        end
    end
end

if ~isfield(cl_aneg_bneg, 'all_data')
    % we have not extracted data successfully; skip by doing nothing
else
    for i = 1:length(cl_aneg_bneg)
        k = size(cl_aneg_bneg(i).all_data, 2);

        cl_aneg_bneg(i).timeseries = cell(1, k);

        for j = 1:k
            cl_aneg_bneg(i).timeseries{j} = cl_aneg_bneg(i).all_data(:, j);
        end
    end
end

if ~isfield(cl_apos_bneg, 'all_data')
    % we have not extracted data successfully; skip by doing nothing
else
    for i = 1:length(cl_apos_bneg)
        k = size(cl_apos_bneg(i).all_data, 2);

        cl_apos_bneg(i).timeseries = cell(1, k);

        for j = 1:k
            cl_apos_bneg(i).timeseries{j} = cl_apos_bneg(i).all_data(:, j);
        end
    end
end

if ~isfield(cl_aneg_bpos, 'all_data')
    % we have not extracted data successfully; skip by doing nothing
else
    for i = 1:length(cl_aneg_bpos)
        k = size(cl_aneg_bpos(i).all_data, 2);

        cl_aneg_bpos(i).timeseries = cell(1, k);

        for j = 1:k
            cl_aneg_bpos(i).timeseries{j} = cl_aneg_bpos(i).all_data(:, j);
        end
    end
end

if dosave
    save intersect_a_and_b_clusters -append cl*
end

end % main function