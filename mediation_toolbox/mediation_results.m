function [cl, z, volInfo] = mediation_results(pthr,kthr,keywd)
%[cl, z, volInfo] = mediation_results(pthr,kthr,keywd)
% keywd is 'direct' or 'mediation'


switch keywd
    case 'direct'
        inname = 'X-Y_direct_pvals.img';
        outname = 'X-Y_direct_pvals_thr.img';
        signimg = 'X-Y_direct_effect.img';
        zname = 'X-Y_direct_zscore_thr.img';

    case 'mediation'
        inname = 'X-M-Y_pvals.img';
        outname = 'X-M-Y_pvals_thr.img';
        signimg = 'X-M-Y_effect.img';
        zname = 'X-M-Y_zscore_thr.img';


    otherwise, error('Unknown keyword')
end

[cl, z, volInfo] = mediation_results_p2z(pthr,kthr,inname,signimg,outname,zname);

if ~isempty(cl)
    cluster_orthviews(cl,'bivalent');
    montage_clusters([],cl,[2 2]);
end

return

