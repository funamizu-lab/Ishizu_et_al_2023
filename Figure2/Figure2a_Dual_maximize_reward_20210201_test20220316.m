%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Figure2a_Dual_maximize_reward_20210201_test20220316
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Noise is determined by gaussian
%With different std and different decision boundary

bin = 0.001;
bin_x = 0:bin:1;

%para_bias: 0.6318    0.6296    0.3851    0.3871
%para_sense: 0.2374    0.3181    0.2182    0.2950

Sense_std = 0.01:0.01:0.5;
reward = [3.8,1]; %left-right
stim_prob = [0.5, 0.5];

std_long  = 0.2574; %average STD for short
std_short = 0.3209; %average STD for long

%%% fig 2a right %%%
cd('G:\upload_code\Figure2\Fig2a');
sdata = struct();
[sdata.long_x, sdata.long_y]=draw_tranc_gauss_psycho05(std_long);
[sdata.short_x,sdata.short_y]=draw_tranc_gauss_psycho05(std_short);
T = struct2table(sdata);
writetable(T, 'source fig2a.csv');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_value,temp_y]=draw_tranc_gauss_psycho05(sense)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bin = 0.01;
x_value = 0:bin:1;

%Tone with 0.5
temp_y = normpdf(x_value,0.5,sense);
temp = normcdf(1,0.5,sense) - normcdf(0,0.5,sense);
temp_y = temp_y ./ temp;

figure;hold on
plot(x_value,temp_y,'k')
plot([0 0],[0 temp_y(1)],'k')
plot([1 1],[0 temp_y(end)],'k')
set(gca,'xlim',[-0.1 1.1])

return

