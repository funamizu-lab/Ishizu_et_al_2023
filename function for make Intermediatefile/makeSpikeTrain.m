function makeSpikeTrain(folders)

analysis_dir = eval(folders);
disp(folders);
for i = 1:length(analysis_dir)
    disp([num2str(i),'/',num2str(length(analysis_dir))]);
    workfolder = analysis_dir{i};
    cd(workfolder);
    mkdir('spike_ch');
    savefolder=[workfolder,'\spike_ch'];
    
    %Make spike_density_function & spike train dataset    
    makeSpikeDensityFunc(savefolder);
end

end

function makeSpikeDensityFunc(savefolder)

load('task_frame_tokyo_ephys_20220210.mat','ave_velocity');
max_time1 = length(ave_velocity); %ms

load('spikedata.mat','spikedata');
max_time2 = zeros(length(spikedata),1);
for i = 1:length(spikedata)
    tmp = spikedata(i).firetiming;
    if(~isempty(tmp))
        max_time2(i) = tmp(end);
    else
    end
end

max_time = max([max_time1; max_time2]);

%Make vector for spike_density_function
SDF_filter = struct();
gauss_std = 100; %ms
width = 50; %ms
for i = 1:max_time
    use_time = (i-width) : (i + width);
    use_time = use_time(use_time > 0 & use_time <= max_time);
    
    gauss_pdf = normpdf(use_time,i,gauss_std);
    gauss_pdf = gauss_pdf ./ sum(gauss_pdf);
    SDF_filter(i).time = use_time;
    SDF_filter(i).gauss_pdf = gauss_pdf;
end

for id = 1:length(spikedata)
    spike_mark = zeros(1,max_time); %ms
    spike_filter= zeros(1,max_time); %ms
    
    spike_mark(spikedata(id).firetiming) =1;
    
    %Based on the spike_mark, make spike density function
    for j = 1:max_time
        spike_filter(j) = sum(spike_mark(SDF_filter(j).time) .* SDF_filter(j).gauss_pdf);
    end
    
    save_file = sprintf('task_spike_stripe20210520_%d.mat',id);
    save_fullfile = fullfile(savefolder,save_file);
    save(save_fullfile, 'spike_mark', 'spike_filter', 'gauss_std');
end
end


% test_mark=spike_mark(1:length(spike_mark_test)) -  spike_mark_test;
% test_filter=spike_filter(1:length(spike_filter_test)) -  spike_filter_test;
% close all; plot(spike_filter); figure; plot(test_filter)