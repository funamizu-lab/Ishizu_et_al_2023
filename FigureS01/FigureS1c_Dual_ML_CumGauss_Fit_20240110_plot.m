%{
----------------------------------------------------------------------------
%Stimulus bias changes by block
%Just model analysis
%Based on Rao 2010 Front Comp Neurosci
%Value updating + Prior updating
----------------------------------------------------------------------------
%}

function FigureS1c_Dual_ML_CumGauss_Fit_20240110_plot

filename1=uigetfile('Dual_ML_CumGauss_Fit_2024011*.mat');
load(filename1)

max_mouse = max(mouse_number);
for i = 1:max_mouse
    temp = find(mouse_number == i);    
    temp_log_likeli = log_likeli_all(temp,:);
    mouse_log_likeli(i,:) = mean(temp_log_likeli);
end

[~,size_x] = size(log_likeli_all);

plot_history = 1:past_trial;


%Make new figure for sabun_log_likeli
for i = 1:size_x-1
    sabun_new_likeli(:,i) = log_likeli_all(:,i+1) - log_likeli_all(:,i);
    sabun_new_mouse_likeli(:,i) = mouse_log_likeli(:,i+1) - mouse_log_likeli(:,i);
end
figure
hold on
plot(sabun_new_mouse_likeli', 'r');
plot_mean_se_moto(sabun_new_likeli(:,plot_history), [0 0 0], 2);
set(gca,'xlim',[0.5 max(plot_history)+0.5])

%log_likeli test
ave_mouse_log = mean(mouse_log_likeli);
for i = 1:size_x-1
    [~,p(i)] = lratiotest(ave_mouse_log(i+1),ave_mouse_log(i),1);
end
p

