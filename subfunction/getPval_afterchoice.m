function [sig_neuron_longtone,sig_neuron_sustain,sig_neuron_choice,sig_neuron_prior,p_tone,p_choice]=...
    getPval_afterchoice(p_threshold,folders)
%
% folders={'auc_ishizu','fof_ishizu','mpfc_ishizu'};
%

analysis_dir = eval(folders);
p_task_long  =cell(length(analysis_dir),1);
p_task2_long =cell(length(analysis_dir),1);
for i = 1:length(analysis_dir)
    [p_task_long{i},p_task2_long{i}] = pVal_depth(analysis_dir{i});   
end
p_task_long2  =cell2mat(p_task_long);
p_task2_long2 =cell2mat(p_task2_long);

p_longtone = (p_task_long2(:,15:25) < p_threshold);
p_tone = (p_task_long2 < p_threshold);
p_prior = (p_task_long2(:,15) < p_threshold);

p_sustain  = (p_task2_long2(:,15:24) < p_threshold);
p_choice= (p_task2_long2< p_threshold);

sig_neuron_longtone = find(sum(p_longtone,2)>0);
use_neuron_longtone = check_nandata(folders);
sig_neuron_longtone = sig_neuron_longtone(use_neuron_longtone);

sig_neuron_prior = intersect(find(p_prior==1),sig_neuron_longtone);

sig_neuron_sustain = find(sum(p_sustain,2)==size(p_sustain,2));

sig_neuron_choice = find(sum(p_choice, 2)>0);
end

function [p_task_long, p_task2_long] = pVal_depth(pathname)

cd(pathname);

temp = dir('sig2_task_neurons_2022*');
load(temp.name,'p_task_long');

temp = dir('sig4_task_neurons_2022*');
load(temp.name,'p_task2_long');

temp = dir('depth_spike_20220517*');
if length(temp) == 1
    load(temp.name);
    depth_neuron = find(spike_depth <= def_depth);
elseif length(temp) == 0
    depth_neuron = 1:size(p_task_long,1); %Use all the neurons
else
end

p_task_long  = p_task_long(depth_neuron,:);
p_task2_long = p_task2_long(depth_neuron,:);
end

function use_neuron =check_nandata(folders)
% function use_neuron =check_nandata(folders,ia,ib)

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
