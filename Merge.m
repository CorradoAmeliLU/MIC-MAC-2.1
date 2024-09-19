
%% PARAMETERS

folder_path = '.\extracted\'; % Output folder.

selected_features = [1:4 6 7 9:12 28 41 42 43 45 48];

feat_names = ["Volume","Volume/nEdges","Average Curvature","Max Curvature","Compactness","Polarity","Link Density","Average Node degree","Mean Edge Length","Max Edge Length","S-Metric","Ending Nodes Density","1st Largest Bound","2nd Largest Bound","1 over 2","Node density"];

%% OPEN DIRECTORY 

addpath(genpath('.\Library'));

dispstat('           -- MIC-MAC 2.0 -- ','keepthis');
dispstat('              -Merge- ','keepthis');
dispstat('          ','keepthis');

fold = dir(folder_path);
[n_objs, ~] = size(fold);
samples_attribute_name = [""];
samples_volume_name = [""];

n_samps_attr = 0;
n_samps_vol = 0;

for ii = 1 : n_objs
    if contains(fold(ii).name,'attributes.mat')
        n_samps_attr = n_samps_attr+1;
        samples_attribute_name(n_samps_attr) = string(fold(ii).name);
    elseif contains(fold(ii).name,'volumes.mat')
        n_samps_vol = n_samps_vol+1;
        samples_volume_name(n_samps_vol) = string(fold(ii).name);
    else
        continue;
    end
end

if n_samps_attr ~= n_samps_vol
    error('For some sample/s at least one file (attributes.mat, volumes.mat) is missing.');
else
    dispstat(['Found ' num2str(n_samps_attr) ' samples.'],'timestamp','keepthis');
    n_samps = n_samps_attr;
    clear n_samps_attr n_samps_vol n_samps_meta n_objs fold
end

%% MERGING ATTRIBUTE DATASET

data = [];
structsPerSample = zeros(1,n_samps);
samples_Conditions = zeros(1,n_samps);

for ii = 1: n_samps
        
    load(strcat(folder_path,samples_attribute_name(ii)));
    
    [structsPerSample(ii),~] = size(attributes); 
      
    data = vertcat(data,attributes);
end

dispstat('Merging Cell Attributes Complete','timestamp','keepthis');

data = data(:,selected_features);

%% EXPORT

writematrix(data)
writematrix(structsPerSample)
writematrix(feat_names)
writematrix(samples_attribute_name)

%% MERGING VOLUME DATASET

dataV = [];

for ii = 1: n_samps
        
    load(strcat(folder_path,samples_volume_name(ii)));
    
     temp = split(samples_volume_name(ii),".");
     sample = temp(1);
     structs.Sample = repmat(sample,structsPerSample(ii),1);

    dataV = vertcat(dataV,structs);
end

dispstat('Merging Cell Volumes Complete','timestamp','keepthis');

samples = unique(dataV.Sample);
