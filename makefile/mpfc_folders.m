
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
function [analysis_dir,short_length] = mpfc_folders

analysis_dir = {
    'G:\Ishizu_data\IntermediateFiles\mpfc\a04\2021-03-04_a04_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\a04\2021-03-05_a04_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\a04\2021-03-06_a04_mPFC_Right_OK\recording2_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\a04\2021-03-07_a04_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\a04\2021-03-15_a04_mPFC_Left_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\a04\2021-03-16_a04_mPFC_Left_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\a04\2021-03-17_a04_mPFC_Left_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\a04\2021-03-18_a04_mPFC_Left_OK\recording1_task'

    'G:\Ishizu_data\IntermediateFiles\mpfc\a08\2021-02-16_a08_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\a08\2021-02-17_a08_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\a08\2021-02-21_a08_mPFC_Left_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\a08\2021-02-22_a08_mPFC_Left_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\a08\2021-02-24_a08_mPFC_Left_OK\recording1_task'
    
    'G:\Ishizu_data\IntermediateFiles\mpfc\i20\2021-05-08_i20_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\i20\2021-05-09_i20_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\i20\2021-05-10_i20_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\i20\2021-05-18_i20_mPFC_Left_OK\recording2_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\i20\2021-05-19_i20_mPFC_Left_OK\recording3_task'
    
    'G:\Ishizu_data\IntermediateFiles\mpfc\i24\2021-08-05_i24_mPFC_Right_OK\recording2_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\i24\2021-08-06_i24_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\i24\2021-09-30_i24_mPFC_left_OK\recording1_task'
    
    'G:\Ishizu_data\IntermediateFiles\mpfc\i34\2021-11-14_i34_mPFC_Left_OK\recording2_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\i34\2021-12-03_i34_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\i34\2021-12-05_i34_mPFC_Right_OK\recording1_task'
    
    'G:\Ishizu_data\IntermediateFiles\mpfc\i35\2021-11-05_i35_mPFC_Left_OK\recording2_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\i35\2021-11-06_i35_mPFC_Left_OK\recording2_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\i35\2021-11-24_i35_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\i35\2021-11-25_i35_mPFC_Right_OK\recording1_task'
    
    'G:\Ishizu_data\IntermediateFiles\mpfc\i43\2022-06-03_i43_mPFC_Right\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\i43\2022-06-04_i43_mPFC_Right\recording4_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\i43\2022-06-24_i43_mPFC_Left\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\i43\2022-06-27_i43_mPFC_Left\recording1_task'
    
    'G:\Ishizu_data\IntermediateFiles\mpfc\i46\2022-06-03_i46_mPFC_Left\recording2_task'
    'G:\Ishizu_data\IntermediateFiles\mpfc\i46\2022-06-04_i46_mPFC_Left\recording1_task'
    };
%short_length = [zeros(8,1);zeros(5,1);zeros(5,1);zeros(3,1);ones(3,1);zeros(4,1)];
short_length = [zeros(8,1);zeros(5,1);zeros(5,1);zeros(3,1);ones(3,1);zeros(4,1);ones(4,1);zeros(2,1)];
if length(analysis_dir) ~= length(short_length)
    hoge
end
