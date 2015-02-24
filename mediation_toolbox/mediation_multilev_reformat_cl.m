% output_cl = mediation_multilev_reformat_cl(input_cl)
%
% This function takes as input a multi-subject cl structure with one cell
% per subject (this is the format returned in clpos_data and clneg_data of
% mediation_brain_results) and re-arranges the data fields so that the
% output clusters structure is a "typical" cl structure with one element
% per region.  The data from each subject is entered in a cell array in the
% cl(x).timeseries field instead.
%
% This means that one can take cl(x).timeseries from the output clusters
% and run mediations on them, or use them as inputs to
% mediation_brain_multilev, directly.


function output_cl = mediation_multilev_reformat_cl(input_cl)

if ~iscell(input_cl)
    help mediation_multilev_reformat_cl
    error('Clusters must be a cell vector, one per subject, of clusters structures.')
end

if ~isfield(input_cl{1}, 'all_data') || ~isfield(input_cl{1}, 'timeseries')
    disp('WARNING! all_data or timeseries fields are missing from clusters.')
    help mediation_multilev_reformat_cl
    disp('Warning: input clusters do not have .all_data and .timeseries fields.')
    
    output_cl = input_cl;
end


nsubjects = length(input_cl);

output_cl = input_cl{1};

% another convenient format
for i = 1:length(output_cl)
    output_cl(i).all_data = cell(1, nsubjects);
    output_cl(i).timeseries = cell(1, nsubjects);
    for s = 1:nsubjects
        if isfield(input_cl, 'all_data')
            output_cl(i).all_data{s} = input_cl{s}(i).all_data;
        else
            output_cl(i).all_data{s} = [];
        end
        
        if isfield(input_cl, 'timeseries')
            output_cl(i).timeseries{s} = input_cl{s}(i).timeseries;
        else
            output_cl(i).timeseries{s} = [];
        end
        
        
    end
end

end