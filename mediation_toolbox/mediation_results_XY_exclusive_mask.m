
% Create exclusively masked images
% Take two p-value images and create image with non-zero p-values where the
% first has low p-values (significant p<.002) and the second has high
% p-values (nonsignificant, p > .05)

maximg = .002;  % p-values must be this low or lower

minmask = .05;  % p-values in EXCLUDE images must be this high or higher

%% B effect, Not A effect

img = 'M-Y_pvals.img'; mask = 'X-M_pvals.img'; outname = 'B_p002_notA_p05.img';
mask_image(img, mask, outname, 'maximg', maximg,  'minmask', minmask);
spm_image('init', outname);

cl = mask2clusters('B_p002_notA_p05.img');
mediation_brain_results('a', 'thresh', [.05], 'size', [1], 'overlay', ovl, 'mask', 'mask.img');
cluster_orthviews(cl, {[0 1 0]},'overlay', ovl, 'add');

wh = (cat(1, cl.numVox) >= 10);
sum(wh)
cl = cl(wh);
cl_bnota = cl;

cl_bnota = mediation_extract_data(cl_bnota);

save cl_B_p002_notA_p05_k10 cl_bnota

cluster_orthviews(cl_bnota, {[0 1 0]},'overlay', ovl, 'solid');

disp('B not A: Report, not stim');
cluster_table(cl_bnota{1}, 1, 0, 'writefile','table_cl_bnota');

%%

img = 'X-M_pvals.img'; mask = 'M-Y_pvals.img'; maximg = .002; minmask = .05; outname = 'A_p002_notB_p05.img';
mask_image(img, mask, outname, 'maximg', maximg,  'minmask', minmask);
%spm_image('init', outname);

cl = mask2clusters('A_p002_notB_p05.img');
%mediation_brain_results('b', 'thresh', [.05], 'size', [1], 'overlay', ovl, 'mask', 'mask.img');
%cluster_orthviews(cl, {[0 1 0]},'overlay', ovl, 'add', 'solid');

wh = (cat(1, cl.numVox) >= 10);
sum(wh)
cl = cl(wh);
cl_anotb = cl;

cl_anotb = mediation_extract_data(cl_anotb);

save cl_A_p002_notB_p05_k10 cl_anotb

cluster_orthviews(cl_anotb{1}, {[1 0 0]},'overlay', ovl, 'solid', 'add');

disp('A not B: Stim, not report');
cluster_table(cl_anotb{1}, 1, 0, 'writefile','table_cl_anotb');
%% Mediators...

% % % % Create an image with non-zero numbers only where p-values in an image are greater than .05
% % % img = 'X-M_pvals.img'; mask = 'X-M_pvals.img'; maxmask = .05; outname = 'notA_p05.img';
% % % mask_image(img, mask, outname, 'reverse', 'maxmask', maxmask);
% % % spm_image('init', outname);
% % % 
% % % img = 'M-Y_pvals.img'; mask = 'M-Y_pvals.img'; maxmask = .05; outname = 'notB_p05.img';
% % % mask_image(img, mask, outname, 'reverse', 'maxmask', maxmask);
% % % spm_image('init', outname);

img = 'X-M-Y_pvals.img'; mask = 'M-Y_pvals.img'; maximg = .002; maxmask = .05; outname = 'tmp.img';
mask_image(img, mask, outname, 'maximg', maximg,  'maxmask', maxmask);
%spm_image('init', outname);

img = 'tmp.img'; mask = 'X-M_pvals.img';  maxmask = .05; outname = 'mediators_p002_p05a_and_b.img';
mask_image(img, mask, outname, 'maxmask', maxmask);

%%


cl = mask2clusters('mediators_p002_p05a_and_b.img');
wh = (cat(1, cl.numVox) >= 10);
sum(wh)
cl = cl(wh);
cl_mediators = cl;

cluster_orthviews(cl_mediators, {[1 1 0]},'overlay', ovl, 'add', 'solid');
montage_clusters(ovl, cl_mediators, {'y'});

%
cl_mediators = mediation_extract_data(cl_mediators);

save cl_mediators_p002_AandBp05_k10 cl_mediators

disp('Mediators and A and B');
cluster_table(cl_mediators{1}, 1, 0);


%% SHORTCUT IF YOU'VE DONE THE ABOVE

load('cl_A_p002_notB_p05_k10.mat')
load('cl_B_p002_notA_p05_k10.mat')
load('cl_mediators_p002_AandBp05_k10.mat')

cluster_orthviews(cl_anotb{1}, {[0 1 0]}, 'overlay', ovl, 'solid');
cluster_orthviews(cl_bnota{1}, {[1 0 0]}, 'overlay', ovl, 'solid', 'add');
cluster_orthviews(cl_mediators{1}, {[1 1 0]}, 'overlay', ovl, 'solid', 'add');

%% XY Sort Plot

load mediation_SETUP
xdata = SETUP.data.X;
ydata = SETUP.data.Y;
cl = {};
for i = 1:length(cl_anotb)
    cl{i} = [cl_anotb{i} cl_bnota{i} cl_mediators{i}];
end

out = mediation_sort_xy_proximity_plot(cl, xdata, ydata)
