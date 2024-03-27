
%{
----------------------------------------------------------------------------
Determine the time window for task relevant neurons
%Start of task
%Sound on
%Sound off
%Before choice
%After choice (0sec)
%After choice (1sec)
%After choice (2sec)
%p = 0.001
%Each epoch, predict the prior, sensory and choice (integration)
----------------------------------------------------------------------------
%}
function process_20220516_make_depth_file(folders)

analysis_dir = eval(folders);
analysis_dir

for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]
    
    [elec_depth, spike_depth, length_neuron] = get_depth_data(analysis_dir{i});
    
    temp = dir('spike_ch*');
    if length(temp) ~= 1
        hoge
    else
        cd(temp.name)
    end
    temp = dir('task_spike_stripe*');
    length_neuron2 = length(temp);
    if length_neuron ~= length_neuron2
        [length_neuron, length_neuron2]
        hoge
    end
    cd ../
    
    %Depth of electrode and limit for the target areas
    %elec_depth
    %Depth of spikes from the tip of electrode
    %meaning that... 20 is the most deep place of electrode (20 becomes 0)
    spike_depth = spike_depth - 20;
    %Subtract with inserted depth
    spike_depth = elec_depth(1) - spike_depth;
    
    def_depth = elec_depth(2);
    
    save depth_spike_20220517 spike_depth def_depth length_neuron
end
delete(gcp('nocreate'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [elec_depth, spike_depth, length_neuron] = get_depth_data(pathname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(pathname)
cd ../ %Get 
elec_depth = dir('depth_info*');
if length(elec_depth) ~= 1
    hoge
else
    elec_depth = load(elec_depth.name);
end

neuron_depth = dir('spike_data_all*');
if length(neuron_depth) ~= 1
    hoge
else
    neuron_depth = load(neuron_depth.name);
end
length_neuron = length(neuron_depth.spike);
for i = 1:length_neuron
    spike_depth(i,1) = neuron_depth.spike(i).depth;
end

cd(pathname) %Back to the original dir
return


