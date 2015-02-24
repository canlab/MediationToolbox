function mediation_X_search_mask

image_names = spm_get(2,'*img','Select X->M and M->Y p-value images.');

[dat,volInfo] = iimg_threshold(image_names,'thr',[0 .05],'outnames','mask.img');

return

