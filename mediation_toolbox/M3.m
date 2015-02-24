% Opening function for the M3 toolbox. For info on the GUI layout, see spm_config_mediation.m

function M3(varargin)
    addpath(fullfile(spm('dir'),'toolbox','M3'));

    spm_jobman('interactive','','jobs.tools.mediation');
    return
end