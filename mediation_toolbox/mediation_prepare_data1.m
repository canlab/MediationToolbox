cd('/Users/tor/Documents/Tor_Documents/PublishedProjects/2005_VNL_visual_nonlinearity/VNL_OUT5')
Lv = getfullpath(filenames('nlfit*Lvis.mat',1));
Rv = getfullpath(filenames('nlfit*Rvis.mat',1));
Lm = getfullpath(filenames('nlfit*Lmot.mat',1));
Rm = getfullpath(filenames('nlfit*Rmot.mat',1));

data.nsubs = size(Lv,1);
data.Lv = Lv; data.Rv = Rv; data.Lm = Lm; data.Rm = Rm;

for i=1:data.nsubs
    clear ROI, load(Lv(i,:)); y = ROI.adjustedy; y(1:860:end) = 0;
    data.data{i}(:,1) = y;
    
    clear ROI, load(Rv(i,:)); y = ROI.adjustedy; y(1:860:end) = 0;
    data.data{i}(:,2) = y;

    clear ROI, load(Lm(i,:)); y = ROI.adjustedy; y(1:860:end) = 0;
    data.data{i}(:,3) = y;
    
    clear ROI, load(Rm(i,:)); y = ROI.adjustedy; y(1:860:end) = 0;
    data.data{i}(:,4) = y;
    
end

data.condf = vnl_condf;

dd = double(vnl_condf > 0);

data.onsets = find(dd);

save vnl_mediation_data data


% get design

[X,px] = trial_level_beta3('onsets',onsets,'rows',6880,'output','design');
data.X = X;
data.px = px;

[X4,px4] = trial_level_beta3('onsets',onsets,'rows',6020,'output','design');
data.X4 = X4;
data.px4 = px4;

% get betas -- h (amplitude) for each trial

for i=1:data.nsubs
    for j = 1:4

        if i == 4
            % one subject has a different model
            [betas,f,X2,px2,h,t,w] = trial_level_beta3('X',X4,'pinvx',px4,'output','betas','data',data.data{i}(:,j));
        else
            [betas,f,X2,px2,h,t,w] = trial_level_beta3('X',X,'pinvx',px,'output','betas','data',data.data{i}(:,j));

        end



        data.betas{i}(:,j) = betas;
        data.h{i}(:,j) = h;
        data.t{i}(:,j) = t;
        data.w{i}(:,j) = w;

    end
end

save vnl_mediation_data data



% get autocorrelation

data.h_xc = [];
for i=1:data.nsubs
    xc = [];
    for j = 1:4
        xc(j,:) = getV('get',data.h{i}(:,j));
    end
    
    if i > 1
        tmp = pad(mean(xc)',data.h_xc(1,:)')';
        data.h_xc(i,:) = tmp;
    else
        data.h_xc(i,:) = mean(xc);
    end
end

% get correlation between height, time, width
data.htw_corr = []; clear c
for i=1:data.nsubs
    for j = 1:4
        c(:,:,j) = corrcoef([data.h{i}(:,j) data.t{i}(:,j) data.w{i}(:,j)]);
    end
    xc = squeeze(mean(c,3));
    data.htw_corr(:,:,i) = xc;
end

[mean,t,sig,out] = ttest3d(data.htw_corr);


% get correlations in data
data.corr = []; clear c
for i=1:data.nsubs
    xc = corrcoef([data.h{i}]);
    data.corr(:,:,i) = xc;
end

[mean,t,sig,out] = ttest3d(data.htw_corr);


% collect data for mediation

for i=1:data.nsubs,mediation_data{i}=data.h{i}(:,[1 3 2]);,end
data.mediation_data = mediation_data;
data.mediation_data_names = {'Lvis' 'Lmot' 'Rvis'}

for i=1:data.nsubs,X{i}=data.h{i}(:,1);,end
for i=1:data.nsubs,Y{i}=data.h{i}(:,3);,end
for i=1:data.nsubs,M{i}=data.h{i}(:,2);,end
data.X = X; data.Y = Y; data.M = M;
data.names = {'L V1' 'L M1' 'R V1'};
% run
[paths,stats2] = mediation_first_level(data.X,data.Y,data.M,'stats','plots','robust','verbose','names',data.names);


for i=1:data.nsubs,X{i}=data.h{i}(:,1);,end
for i=1:data.nsubs,Y{i}=data.h{i}(:,4);,end
for i=1:data.nsubs,M{i}=data.h{i}(:,3);,end
data.X2 = X; data.Y2 = Y; data.M2 = M;
data.names2 = {'L V1' 'R M1' 'L M1'};
[paths,stats2] = mediation_first_level(data.X2,data.Y2,data.M2,'stats','plots','robust','verbose','names',{'L V1' 'R M1' 'L M1'});
