function [cl, z, volInfo] = mediation_results_p2z(pthr,kthr,inname,signname,threshname,zname)
% pthr = .005; kthr = 5;
% inname = 'X-Y_direct_pvals.img'; outname = 'X-Y_direct_pvals_thr.img'; 
% signimg = 'X-Y_direct_effect.img';
% zname = 'X-Y_direct_zscore_thr.img';
%
% [cl, z, volInfo] = mediation_results_p2z(pthr,kthr,inname,signimg,outname,zname);

volInfo = iimg_read_img(inname,2);
[dat,volInfo2] = iimg_threshold(inname,'thr',[0 pthr],'k',kthr,'outnames',threshname);
fprintf(1,'\nVoxels: %3.0f\n',sum(dat>0));

% get z-scores
dat(dat==0) = 1; z = norminv(1 - dat); z(isinf(z)) = 0;

% get sign info
dat2 = iimg_mask(dat,signname,volInfo);
z = z .* sign(dat2);

% write output
voldata = iimg_reconstruct_3dvol(z,volInfo,'outname',zname);
cl = mask2clusters(voldata,volInfo.mat);

if ~isempty(cl)
sz = cat(1,cl.numVox);
cl(sz < kthr) = [];
end

fprintf(1,'Suprathreshold clusters: %3.0f',length(cl));

return