function [sig_depth, sig_raw2, rawNeuronLength] = getPval_long(p_threshold,folders)
%
% folders={'auc_ishizu','fof_ishizu','mpfc_ishizu'};
%

analysis_dir = eval(folders);
p_task_long=cell(length(analysis_dir),1);
depth_neuron=cell(length(analysis_dir),1);
rawNeuronLength=zeros(length(analysis_dir),1);
for i = 1:length(analysis_dir)
    output = pVal_depth(analysis_dir{i}); 
    p_task_long{i}=output.p_long; 
    depth_neuron{i}  = output.depth;
    rawNeuronLength(i)=length(output.depth);
end
p_task_long2=cell2mat(p_task_long);
depth_ALL = cell2mat(depth_neuron);
p_task_tmp = p_task_long2(depth_ALL==1,:);
p_task = (p_task_tmp(:,15:25) < p_threshold);
sig_depth = find(sum(p_task,2)>0);
use_neuron= check_nandata(folders);
sig_depth = sig_depth(use_neuron);

% use_neuron_longtone = check_nandata(folders);
% sig_neuron_longtone = tasksig_all(use_neuron_longtone);
sig_raw = zeros(length(depth_ALL),1);
tmp=find(depth_ALL==1);
sig_raw(tmp(sig_depth))=1;
sig_raw2 = find(sig_raw==1);
end

function  output =  pVal_depth(pathname)

cd(pathname);

temp = dir('sig2_task_neurons_2022*');
load(temp.name,'p_task_long');

temp = dir('depth_spike_20220517*');
if length(temp) == 1
    load(temp.name);
    depth_neuron = find(spike_depth <= def_depth);
elseif length(temp) == 0
    depth_neuron = 1:size(p_task_long,1); %Use all the neurons
else
end

% output.p_long=p_task_long(depth_neuron,:);
output.p_long=p_task_long;
depth = zeros(size(p_task_long,1),1);
depth(depth_neuron)=1;
output.depth=depth;
end

function use_neuron =check_nandata(folders)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determine the sig_neuron at certain time window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sig_time_window = 15; %Before sound neurons 
%sig_time_window = 25; %Only during sound neurons
sig_time_window = 0; %Before + During sound neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

analysis_dir = eval(folders);
use_frame = 6:40;

for i = 1:40
    for j = 1:6
        evi_prefer_all(i,j).matrix = [];
        evi_nonprefer_all(i,j).matrix = [];
    end
end
for i = 1:length(analysis_dir)
   
    [~,~,~,~,~,~,~,~,~,~,~,~,evi_prefer,evi_nonprefer,...
        ~,~,~,~,~,~,~,~,~,~,~,~] = ...
        Task_kaiseki_tokyo1_20230711_sound_choice_process4A_depth(analysis_dir{i},sig_time_window);
    
    for j = 1:40      
        for k = 1:6
            evi_prefer_all(j,k).matrix = [evi_prefer_all(j,k).matrix; evi_prefer(j,k).matrix];  %Correct trials with evidence prefer block
            evi_nonprefer_all(j,k).matrix = [evi_nonprefer_all(j,k).matrix; evi_nonprefer(j,k).matrix]; %Correct trials with evidence nonprefer block
        end
    end
end
% delete(gcp('nocreate'))

for j = 1:length(use_frame)
    temp_frame = use_frame(j);    
    %Get the nan frame
    for k = 1:6
        %moto_data = evi_trace_all(temp_frame,k).matrix;
        %nan_check(:,k) = isnan(moto_data); %detect_non_nan
        moto_data = evi_prefer_all(temp_frame,k).matrix;
        nan_check_prefer(:,k) = isnan(moto_data); %detect_non_nan
        moto_data = evi_nonprefer_all(temp_frame,k).matrix;
        nan_check_nonprefer(:,k) = isnan(moto_data); %detect_non_nan
    end
    nan_check_prefer = max(nan_check_prefer,[],2);
    nan_check_nonprefer = max(nan_check_nonprefer,[],2);
    nan_check = max([nan_check_prefer,nan_check_nonprefer],[],2);
    use_neuron = find(nan_check == 0);
end
end
