
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
function Figure6c_process_20230920_0915_SLR_long_short_integrate(folders)

%frame_long(i,:)  = [1 2 6 7 8 9]; %before,init,end,after1,after2,choice
[pre_s_choice, pre_l_choice] = ...
    process_20220520_0912_SLR_20230915_choice_long_short(folders);
close all

%selected_window = [1:8];%4 is choice onset
[post_s_choice, post_l_choice] = ...
    process_20220520_0912_SLR_20230920_choice_long_short(folders);
close all

%Select the use numbers
%last 200ms of sound
integ_s_choice(:,1) = pre_s_choice(:,3); 
integ_l_choice(:,1) = pre_l_choice(:,3); 

%-0.4 - -0.2 from choice
integ_s_choice(:,2) = post_s_choice(:,2); 
integ_l_choice(:,2) = post_l_choice(:,2); 
%-0.2 - 0 from choice
integ_s_choice(:,3) = post_s_choice(:,3); 
integ_l_choice(:,3) = post_l_choice(:,3); 
%0 - 0.2 from choice
integ_s_choice(:,4) = post_s_choice(:,4); 
integ_l_choice(:,4) = post_l_choice(:,4); 

%% Fig 6c
cd('G:\upload_code\Figure6\Fig6c');
name = extractBefore(folders,'_');

figure; hold on
sdata = struct();% source data 
[median_trace,~,se_trace]= plot_median_se_moto(integ_s_choice, [0 0 1], 2);
sdata.short_median = median_trace';
sdata.short_se = se_trace';

[median_trace,~,se_trace]= plot_median_se_moto(integ_l_choice, [1 0 0], 2);
sdata.long_median = median_trace';
sdata.long_se = se_trace';
set(gca,'xlim',[0.5 4.5])
set(gca,'ylim',[0.65 1])

for i = 1:4
    signrank(integ_s_choice(:,i),integ_l_choice(:,i))
end

T = struct2table(sdata);
writetable(T, ['source fig6c ',name,'.csv']);
