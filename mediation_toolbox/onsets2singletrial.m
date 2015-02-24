function [onsets2 runlabels eventlabels rundesign eventdesign ratings] = onsets2singletrial(onsets,nscans,timeval,TR,eventnames, varargin)
    
    %Converts 'structure' format onsets to 'cell' format suitable for
    %single-trial analysis. Note that the output will always be in TRs.
    %
    %USAGE: [onsets2 runlabels eventlabels rundesign eventdesign] = onsets2singletrial(onsets,nscans,TR,[timeval],[eventnames], ...
    %       ['noaddcumulative'],['ratingfields',ratingfields]);
    %
    % e.g.,
    % eventnames{1} = {'warm' 'low' 'medium' 'high'};
    % ratingnames{1} = {'warmrate' 'lowrate' 'medrate' 'highrate'};
    % [onsets2 runlabels eventlabels rundesign eventdesign ratings] = ...
    % onsets2singletrial(subjects,EXPT.SPM.nscan,'TR',2,eventnames,'noadd','ratingnames',ratingnames);
    %
    %
    %INPUTS:
    %   onsets - onsets in 'structure' format (see below)
    %   nscans - vector containing the number of scans in each run
    %   TR - you should know what this means
    %   timeval (optional) - 'sec' or 'TR' (default: 'TR')
    %   eventnames (optional) - this feature gives you flexibility with regard
    %   to which events you want to include. Ultimately, eventnames will be a cell
    %   array where the rows correspond to subjects and the columns correspond to
    %   runs, with each cell containing a cell-vector of strings (the names of the
    %   event fields in the onsets structure). However, you can specify it several ways:
    %       (1) If eventnames is a row vector containing the event cells for
    %       each run for a single subject, it will replicate this for all subs
    %       (2) If eventnames is a column vector containing the event cells for
    %       each subject for a single run, it will replicate this for all runs
    %       (3) If there is only 1 cell, it will assume all runs/subs have the
    %       same events
    %       (4) The most detailed specification is to feed in the complete cell
    %       array (as described above); this is the default
    %
    %OUTPUTS:
    %   onsets2 - onsets in single-trial cell format
    %   runlabels - labels for each onset specifying which run it's in
    %   eventlabels - labels for each onset specifying which event-type it is
    %                 (numbering is based on the ordering of names in eventnames)
    %   rundesign - design matrix containing run regressors
    %   eventdesign - design matrix containing event-type regressors
    %       (these can be combined simply by typing: [eventdesign rundesign])
    %
    %EXAMPLE: [onsets2 runlabels eventlabels rundesign eventdesign] = onsets2singletrial(onsets,[175 175 175 175 175 175],'sec',2,reward_eventnames);
    %
    %A NOTE ABOUT ONSETS:
    %The SCAN lab uses two different formats for its onsets. The
    %'structure' format contains a different field for each subject, and a
    %different sub-field for each run, with sub-fields for each condition that
    %contain the onsets. The 'cell' format is a cell vector, where each cell
    %contains 1 subject; within each subject-cell is another cell vector, where
    % each cell contains onsets for a particular run/condition. If, for example, you have 2
    %conditions and 2 runs, the cell array will be organized as follows:
    %The 1st cell will contain the onsets of all events from condition 1 in
    %run 1, the 2nd cell will contain the onsets of events from condition 2 in
    %run 1, the 3rd cell will contain onsets of events from condition 1 in run
    %2, etc. Specifically for single-trial analysis, the onsets must be in cell
    %format (and in TRs), with the added requirement that the onsets be adjusted
    %so that the first onset of each run is equal to the onset of the last scan 
    %(*not event*) of the previous run. For example, if you have 100 scans in
    %each run, the first onset of the 1st run should be 0, the first onset of
    %the 2nd run should be 100, etc.
    %
    % More optional inputs:
    % 'noaddcumulative' : if entered, DO NOT adjust onsets to be in times
    % since start of session (they should be in times since start of
    % session for single trial analysis, but sometimes they may already
    % be.)
    %
    % 'ratingnames', followed by field names :
    % if entered, sorts ratings by the same order as onsets, so that
    % ratings/other behavioral/physio/etc. data entered is in same format
    % as onsets.  data should be one observation per event.  enter cell
    % array of field names here, same format as eventnames
    % 
    % By Sam Gershman, modified 3/07 by Tor Wager
    %
    % Examples:
    % do for multi-subject onsets structure; do not adjust times (already
    % done)
    % [onsets2 runlabels eventlabels rundesign eventdesign] = onsets2singletrial(subjects,EXPT.SPM.nscan,'TR',2,eventnames,'noadd');
    
    addcumulative = 1;
    for i = 1:length(varargin)
        if ischar(varargin{i})
            
            switch varargin{i}
                case {'noadd', 'noaddcumulative'}, addcumulative = 0;
                    
                case {'ratingnames', 'ratingfields'}, ratingnames = varargin{i+1};
                    
                otherwise 
                    fprintf('Unknown string input: %s\n', varargin{i})
            end
        end
    end
    
    submat = fieldnames(onsets);    %subject names
    nsubs = length(submat);         %number of subjects
    nruns = length(nscans);         %number of runs
    
    %set up event-names (and ratings if entered)
    
get_eventnames;

onsets2 = cell(1, nsubs);
runlabels = cell(1, nsubs);
eventlabels = cell(1, nsubs);

    %get onsets for this subject
    for s = 1:nsubs
        
        onsets2{s} = []; runlabels{s} = []; eventlabels{s} = []; ratings{s} = [];
        
        for r = 1:nruns
            
            runstr = ['run' num2str(r)];
            evt = eventnames{s,r};
            
            if exist('ratingnames','var')
                rat = ratingnames{s,r};
            end
            
            for e = 1:length(evt);
                
              
                % onsets for this condition
                ons = onsets.(submat{s}).(runstr).(evt{e});
                
                if exist('ratingnames','var')
                    rating_this_evt = onsets.(submat{s}).(runstr).(rat{e});
                end
                %ons = eval(['onsets.',submat{s},'.run{',num2str(r),'}.',evt{e}]);   %all onsets for subject/run/event-type
                
                switch lower(timeval)
                    case 'sec'
                        ons = round(ons / TR);
                        
                    case 'tr'
                        % do nothing
                        
                    otherwise, error('timeval must be ''sec'' or ''tr'''); 
                end
                        
                %if strcmp(timeval,'sec'); ons = round(ons / TR); end;               %if in seconds, convert to TR
                onsets2{s} = cat(1,onsets2{s},ons);
                runlabels{s} = cat(1,runlabels{s},repmat(r,length(ons),1));
                eventlabels{s} = cat(1,eventlabels{s},repmat(e,length(ons),1));
                
                if exist('ratingnames','var')
                    ratings{s} = cat(1,ratings{s},rating_this_evt);
                end
                
            end
        end
    
        if addcumulative
            % adjust so that we have cum. times from start of scanning

            [onsets2{s} eventlabels{s} ind] = adjust_onsets(onsets2{s},runlabels{s},eventlabels{s},nscans); 
            %adjust onset timing
           
                
        else
            % don't adjust, already done
            [onsets2{s},ind] = sort(onsets2{s});
            eventlabels{s} = eventlabels{s}(ind);
        end
                
         if exist('ratingnames','var')
             ratings{s} = ratings{s}(ind);
         end
         
        rundesign{s} = label2design(runlabels{s}); 
        eventdesign{s} = label2design(eventlabels{s});

    end
    
    
% ------------------------------------------------
% inline
% ------------------------------------------------

    function get_eventnames

        doratings = 0;
        if exist('ratingnames', 'var')
            doratings = 1;
            disp('  Found ratings/other data...will sort in format of onsets2');
        end

        if exist('eventnames', 'var')      %event-names are provided
            
            if numel(eventnames) == 1
                eventnames = repmat(eventnames,nsubs,nruns);

                if doratings
                    ratingnames = repmat(ratingnames,nsubs,nruns);
                end
            end

            if isvector(eventnames)

                if size(eventnames,1) == 1;      %has the run info for 1 subject; repeat it for all subjects
                    eventnames = repmat(eventnames,nsubs,1);

                    if doratings
                        ratingnames = repmat(ratingnames,nsubs,1);
                    end

                elseif size(eventnames,2) == 1;  %has the info for 1 run for all subjects; repeat it for all runs
                    eventnames = repmat(eventnames,1,nruns);

                    if doratings
                        ratingnames = repmat(ratingnames,1,nruns);
                    end
                end
            end

        else              %get the event-names automatically

            for s = 1:nsubs
                for r = 1:nruns
                    eventnames{s,r} = fieldnames(eval(['onsets.',submat{s},'.run{',num2str(r),'}']));

                    if doratings
                        error('Impossible...cannot enter rating field names without specifying onset field names explicitly.');
                    end
                end
            end
        end

    end

    
end   % end main function



function [onsets eventlabels ind] = adjust_onsets(onsets,runlabels,eventlabels,nscans)

    %Adjusts timing of onsets so that the first onset of each run corresponds to
    %the timing of the last scan of the previous run. See onsets2singletrial.m
    
    nruns = length(nscans);

    num2add = [0 cumsum(nscans)]; 
    num2add = num2add(1:end-1);

    onsets = onsets + num2add(runlabels')';

    [onsets,ind] = sort(onsets);
    eventlabels = eventlabels(ind);
    
end

    
function design = label2design(labels)
    lab = unique(labels);
    nlab = length(lab);
    design = zeros(length(labels),nlab);
    for lb = 1:length(labels);
        design(lb,labels(lb)) = 1;
    end
end