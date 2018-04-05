function [dispobj, f] = plotPDM(pdm,maskobj)

% plots thresholded PDM maps on brain slices using
% canlab_results_fmridisplay
%
% **Input:**
%
%   **pdm*
%       PDM structure with thresholded results (thresholdPDM.m)
%         
%   **maskobj**
%       an fmri_data-object in the same space/dim as the PDMs
%         



WfPlot = [];
WfLabel = [];

if isfield(pdm,'WfullThresh')
    WfPlot  = pdm.WfullThresh;
    WfLabel = strsplit(strtrim(sprintf('PDM%2d\n',1:numel(WfPlot))),'\n');
end
if isfield(pdm,'WfullJointThresh')
    WfPlot = horzcat(WfPlot, pdm.WfullJointThresh);
    WfLabel = horzcat(WfLabel,'Joint PDM');
end

% plot PDMs
for k=1:numel(WfPlot)
   
    pdmplot = maskobj;
    pdmplot.dat = WfPlot{k};
    
    f(k) = figure(k); clf; set(f(k),'color','w','name',WfLabel{k});
    
    dispobj{k} = canlab_results_fmridisplay(pdmplot);
end