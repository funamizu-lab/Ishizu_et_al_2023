
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
function [analysis_dir,short_length] = fof_folders

analysis_dir = {
    'G:\Ishizu_data\IntermediateFiles\fof\a04\2021-02-27_a04_FOF_Right_chunkerror_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\fof\a04\2021-02-28_a04_FOF_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\fof\a04\2021-03-01_a04_FOF_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\fof\a04\2021-03-10_a04_FOF_Left_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\fof\a04\2021-03-11_a04_FOF_Left_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\fof\a04\2021-03-12_a04_FOF_Left_OK\recording1_task'
    
    'G:\Ishizu_data\IntermediateFiles\fof\a08\2021-02-12_a08_FOF_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\fof\a08\2021-02-13_a08_FOF_Right_OK\recording2_task'
    'G:\Ishizu_data\IntermediateFiles\fof\a08\2021-02-19_a08_FOF_Left_OK\recording2_task'
    'G:\Ishizu_data\IntermediateFiles\fof\a08\2021-02-20_a08_FOF_Left_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\fof\a08\2021-02-23_a08_FOF_Left_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\fof\a08\2021-02-26_a08_FOF_Left_OK\recording3_task'
    
    'G:\Ishizu_data\IntermediateFiles\fof\i20\2021-04-29_i20_FOF_Right_chunkerror_OK\recording4_task'
    'G:\Ishizu_data\IntermediateFiles\fof\i20\2021-04-30_i20_FOF_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\fof\i20\2021-05-13_i20_FOF_Left_OK\recording2_task'
    'G:\Ishizu_data\IntermediateFiles\fof\i20\2021-05-14_i20_FOF_Left_chunkerror_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\fof\i20\2021-05-15_i20_FOF_Left_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\fof\i20\2021-05-16_i20_FOF_Left_OK\recording2_task'
    
    'G:\Ishizu_data\IntermediateFiles\fof\i24\2021-08-03_i24_FOF_Right_OK\recording2_task'
    'G:\Ishizu_data\IntermediateFiles\fof\i24\2021-08-04_i24_FOF_Right_OK\recording2_task'
    'G:\Ishizu_data\IntermediateFiles\fof\i24\2021-08-07_i24_FOF_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\fof\i24\2021-09-29_i24_FOF_left_OK\recording2_task'
    'G:\Ishizu_data\IntermediateFiles\fof\i24\2021-10-01_i24_FOF_left_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\fof\i24\2021-10-02_i24_FOF_left_OK\recording1_task'
    
    'G:\Ishizu_data\IntermediateFiles\fof\i34\2021-11-10_i34_FOF_left_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\fof\i34\2021-12-01_i34_FOF_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\fof\i34\2021-12-02_i34_FOF_Right_OK\recording2_task'
    
    'G:\Ishizu_data\IntermediateFiles\fof\i35\2021-11-03_i35_FOF_Left_OK\recording2_task'
    'G:\Ishizu_data\IntermediateFiles\fof\i35\2021-11-04_i35_FOF_Left_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\fof\i35\2021-11-26_i35_FOF_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\fof\i35\2021-11-27_i35_FOF_Right_OK\recording1_task'
    
    'G:\Ishizu_data\IntermediateFiles\fof\i43\2022-06-22_i43_FOF_Left\recording2_task'
    
    'G:\Ishizu_data\IntermediateFiles\fof\i46\2022-06-01_i46_FOF_Left\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\fof\i46\2022-06-02_i46_FOF_Left\recording2_task'
    };
%short_length = [zeros(6,1);zeros(6,1);zeros(6,1);zeros(6,1);ones(3,1);zeros(4,1)];
short_length = [zeros(6,1);zeros(6,1);zeros(6,1);zeros(6,1);ones(3,1);zeros(4,1);ones(1,1);ones(1,1);zeros(1,1)];
if length(analysis_dir) ~= length(short_length)
    hoge
end