function [paths, means, stes, pvals] = dcm_sim(SPM, psych_input, X, M, Y)
    num_subjects = size(X, 2);
    for i=1:num_subjects
        XVOI = fakeVOI(X(:,i), SPM, 'X', 'mm_location', [0 -30 0], 'filter_data', 0);
        MVOI = fakeVOI(M(:,i), SPM, 'M', 'mm_location', [0 0 20], 'filter_data', 0);
        YVOI = fakeVOI(Y(:,i), SPM, 'Y', 'mm_location', [0 30 0], 'filter_data', 0);

        DCM_filename = save_test_DCM(psych_input, SPM, XVOI, MVOI, YVOI);
        spm_dcm_estimate(DCM_filename);
        load(DCM_filename);

        a = DCM.A(2,1);
        b = DCM.A(3,2);
        c = DCM.A(3,1);

        paths(i,:) = [a b c a*b];
    end
    means = mean(paths);
    stes = std(paths) / sqrt(num_subjects);
    [h, pvals] = ttest(paths);
end



function DCM_filename = save_test_DCM(psych_input, SPM, XVOI, MVOI, YVOI) %taken from spm_dcm_ui
    swd = pwd();
    name  = 'sim_last_run';

    % outputs
    %===================================================================

    % get cell array of region structures
    %-------------------------------------------------------------------
    m     = 3;

    xY(1) = XVOI;
    xY(2) = MVOI;
    xY(3) = YVOI;

    % inputs
    %===================================================================

    % get 'causes' or inputs U
    %-------------------------------------------------------------------
    Sess   = SPM.Sess(xY(1).Sess);
    U.dt   = Sess.U(1).dt;

    U.name = {'psych_input'};
    U.u = interp(psych_input, SPM.xY.RT/U.dt);
    U.u = U.u(1:(length(psych_input)*SPM.xY.RT/U.dt));


    % graph connections
    %===================================================================
    n     = size(U.u,2);
    a = tril(ones(3,3));
    c = [1 0 0]';
    b     = zeros(m,m,n);
    dx    = 20;

    % slice timing
    %-------------------------------------------------------------------
    delays = SPM.xY.RT*ones(1, m);

    % confounds (NB: the data have been filtered and whitened)
    %-------------------------------------------------------------------
    v     = size(xY(1).u,1);
    X0    = xY(1).X0;

    % response variables
    %-------------------------------------------------------------------
    n     = length(xY);
    Y.dt  = SPM.xY.RT;
    Y.X0  = X0;
    for i = 1:n

        % regional responses
        %-----------------------------------------------------------
        Y.y(:,i)  = xY(i).u;
        Y.name{i} = xY(i).name;
    end

    % error precision components (one for each region) - i.i.d. (because of W)
    %-------------------------------------------------------------------
    Y.Q = spm_Ce(ones(1,n)*v);

    DCM.a=a;
    DCM.b=b;
    DCM.c=c;
    DCM.U=U;
    DCM.Y=Y;
    DCM.xY=xY;
    DCM.v=v;
    DCM.n=n;
    DCM.delays = delays;

    %-Save and reset title
    %-------------------------------------------------------------------
    DCM_filename = fullfile(swd,['DCM_' name]);

    if spm_matlab_version_chk('7') >= 0
        save(DCM_filename,'-V6','DCM');
    else
        save(DCM_filename,'DCM');
    end;
end