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

mask_obj = fmri_data(mask, 'noverbose');

effect_obj = fmri_data(char('X-M_effect.img', 'M-Y_effect.img'), mask,'noverbose');
p_obj = fmri_data(char('X-M_pvals.img', 'M-Y_pvals.img'), mask,'noverbose');
%p_obj.dat(p_obj.dat == 0) = 1;  % resampling will otherwise make some 0 values look valid

% apply mask
% mask_obj = resample_space(mask_obj, effect_obj);
% effect_obj = apply_mask(effect_obj, mask_obj);
% p_obj = apply_mask(p_obj, mask_obj);

% compute overlap
bothsig = all(double(p_obj.dat) < p & double(p_obj.dat) > 0, 2); % vox where p-vals are low enough in both images

bothpos = bothsig & all(effect_obj.dat > 0, 2);
obj = effect_obj;
obj.dat = single(bothpos);
r = region(obj, 'noverbose');
wh_omit = cat(1, r.numVox) < k;
r(wh_omit) = [];
cl_apos_bpos = r;

bothneg = bothsig & all(effect_obj.dat < 0, 2);
obj = effect_obj;
obj.dat = single(bothpos);
r = region(obj, 'noverbose');
wh_omit = cat(1, r.numVox) < k;
r(wh_omit) = [];
cl_aneg_bneg = r;

aposbneg = bothsig & effect_obj.dat(:, 1) > 0 & effect_obj.dat(:, 2) < 0;
obj = effect_obj;
obj.dat = single(aposbneg);
r = region(obj, 'noverbose');
wh_omit = cat(1, r.numVox) < k;
r(wh_omit) = [];
cl_apos_bneg = r;

bposaneg = bothsig & effect_obj.dat(:, 2) > 0 & effect_obj.dat(:, 1) < 0;
obj = effect_obj;
obj.dat = single(bposaneg);
r = region(obj, 'noverbose');
wh_omit = cat(1, r.numVox) < k;
r(wh_omit) = [];
cl_aneg_bpos = r;

%
% effect_obj = resample_space(effect_obj, mask_obj);
% p_obj = resample_space(p_obj, mask_obj);

%% Plot

if doslices
    
    % Plot 1
    
    create_figure('overlap');
    axis off
    o2 = canlab_results_fmridisplay;
    
    o2 = montage(cl_apos_bpos, 'o2', o2, 'color', [1 0 0]);
    o2 = montage(cl_aneg_bneg, 'o2', o2, 'color', [0 0 1]);
    
    o2 = montage(cl_apos_bneg, 'o2', o2, 'color', [1 1 0]);
    o2 = montage(cl_aneg_bpos, 'o2', o2, 'color', [0 1 0]);
    
    if dosave
        saveas(gcf, ['intersect_' mystr '_a_and_b'], 'png')
    end

    % Plot 2
    
    o2 = montage(cl_apos_bpos, 'color', [1 0 0], 'montagetype', 'regioncenters');
      
    if dosave
        saveas(gcf, ['intersect_' mystr '_aposbpos_regioncenters'], 'png')
    end
    
    % Plot 3
    
    o2 = montage(cl_aneg_bneg, 'color', [0 0 1], 'montagetype', 'regioncenters');
    
    if dosave
        saveas(gcf, ['intersect_' mystr '_anegbneg_regioncenters'], 'png')
    end
    
    
    
end


%% Save

if dosave
    
    pstr = strrep(num2str(p), '0.', '');
    
    myfilename = sprintf('intersect_a_and_b_clusters_%s_k%d', pstr, k);
    
    fprintf('Saving clusters with extracted data in:\n%s\n', myfilename)
    
    save(myfilename, 'cl_*');
    
end

%% Extract

% Different for single or multi-level
if isfield(SETUP, 'M') % single
    
    is_single_level = true;
    
elseif  isfield(SETUP, 'data') && isfield(SETUP.data, 'M') % multi
    
    is_single_level = false;
else
    error('Cannot find either SETUP.M or SETUP.data.M. Incomplete mediation results?')
end

if ~is_single_level
    disp('Skipping data extraction from multi-level analysis');
    return
end

switch SETUP.cmdstring
    case 'Search for mediators'
        data_obj = fmri_data(SETUP.M, mask, 'noverbose');
        
    case 'Search for indirect influences'
        data_obj = fmri_data(SETUP.X, mask, 'noverbose');
        
    case 'Search for mediated outcomes'
        data_obj = fmri_data(SETUP.Y, mask, 'noverbose');
    otherwise
        error('Unknown cmd string "%s".', cmdstring);
end

% resample to mask to make sure in same space

data_obj = resample_space(data_obj, mask_obj, 'noverbose');

cl_apos_bpos = extract_data(cl_apos_bpos, data_obj);
cl_aneg_bneg = extract_data(cl_aneg_bneg, data_obj);

cl_apos_bneg = extract_data(cl_apos_bneg, data_obj);
cl_aneg_bpos = extract_data(cl_aneg_bpos, data_obj);

  
if dosave
    save intersect_a_and_b_clusters -append cl*
end


%%
% 
% % fix for diffs in how single-level and multi-level data are stored
% if ~isfield(SETUP, 'data')
%     SETUP.data.M = {SETUP.M};
% end
% 
% if ~isempty(cl_apos_bpos)
%     cl_apos_bpos = extract_raw_data([], cl_apos_bpos, 'extract_from', SETUP.data.M, 'noraw');
% end
% 
% if ~isempty(cl_aneg_bneg)
%     cl_aneg_bneg = extract_raw_data([], cl_aneg_bneg, 'extract_from', SETUP.data.M, 'noraw');
% end
% 
% if ~isempty(cl_apos_bneg)
%     cl_apos_bneg = extract_raw_data([], cl_apos_bneg, 'extract_from', SETUP.data.M, 'noraw');
% end
% 
% if ~isempty(cl_aneg_bpos)
%     cl_aneg_bpos = extract_raw_data([], cl_aneg_bpos, 'extract_from', SETUP.data.M, 'noraw');
% end
% 
% if ~isfield(cl_apos_bpos, 'all_data')
%     % we have not extracted data successfully; skip by doing nothing
% else
%     for i = 1:length(cl_apos_bpos)
%         k = size(cl_apos_bpos(i).all_data, 2);
% 
%         cl_apos_bpos(i).timeseries = cell(1, k);
% 
%         for j = 1:k
%             cl_apos_bpos(i).timeseries{j} = cl_apos_bpos(i).all_data(:, j);
%         end
%     end
% end
% 
% if ~isfield(cl_aneg_bneg, 'all_data')
%     % we have not extracted data successfully; skip by doing nothing
% else
%     for i = 1:length(cl_aneg_bneg)
%         k = size(cl_aneg_bneg(i).all_data, 2);
% 
%         cl_aneg_bneg(i).timeseries = cell(1, k);
% 
%         for j = 1:k
%             cl_aneg_bneg(i).timeseries{j} = cl_aneg_bneg(i).all_data(:, j);
%         end
%     end
% end
% 
% if ~isfield(cl_apos_bneg, 'all_data')
%     % we have not extracted data successfully; skip by doing nothing
% else
%     for i = 1:length(cl_apos_bneg)
%         k = size(cl_apos_bneg(i).all_data, 2);
% 
%         cl_apos_bneg(i).timeseries = cell(1, k);
% 
%         for j = 1:k
%             cl_apos_bneg(i).timeseries{j} = cl_apos_bneg(i).all_data(:, j);
%         end
%     end
% end
% 
% if ~isfield(cl_aneg_bpos, 'all_data')
%     % we have not extracted data successfully; skip by doing nothing
% else
%     for i = 1:length(cl_aneg_bpos)
%         k = size(cl_aneg_bpos(i).all_data, 2);
% 
%         cl_aneg_bpos(i).timeseries = cell(1, k);
% 
%         for j = 1:k
%             cl_aneg_bpos(i).timeseries{j} = cl_aneg_bpos(i).all_data(:, j);
%         end
%     end
% end

end % main function