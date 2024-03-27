
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
%Compare sound, choice prior during long sound (-100 ms to 1500 ms for long)
%Compare sound, choice prior during short sound (-100 ms to 700 ms for short)
%p = 0.001
%Each epoch, predict the prior, sensory and choice (integration)
----------------------------------------------------------------------------
%}
function Task_kaiseki_tokyo_encode_20220912_all_all(folders)
%p = genpath('/home/funamizu/Tokyo_ephys/Tokyo_ephys_program')
%addpath(p);

analysis_dir = eval(folders);
analysis_dir

for i = 1:length(analysis_dir)
    Task_kaiseki_tokyo_encode_20220912_glmnet(analysis_dir{i});
    Task_kaiseki_tokyo_encode_20220912_glmnet_WO_choice(analysis_dir{i});
    Task_kaiseki_tokyo_encode_20220912_glmnet_WO_prior(analysis_dir{i});
    Task_kaiseki_tokyo_encode_20220912_glmnet_WO_sound(analysis_dir{i});
end