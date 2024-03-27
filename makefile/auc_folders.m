
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
function [analysis_dir,short_length] = auc_folders

analysis_dir = {
    'G:\Ishizu_data\IntermediateFiles\auditory\a04\2021-02-13_a04_1PPC2AC_Left_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\a04\2021-02-18_a04_1AC2PPC_Right\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\a04\2021-02-19_a04_1AC2PPC_Right\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\a04\2021-02-20_a04_AC_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\a04\2021-02-21_a04_AC_Right_OK\recording1_task'

    'G:\Ishizu_data\IntermediateFiles\auditory\a08\2021-02-09_a08_1PPC2AC_Left_OK\recording2_task'
    
    'G:\Ishizu_data\IntermediateFiles\auditory\i24\2021-10-05_i24_AC_Left_chunk_error_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\i24\2021-10-06_i24_AC_left_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\i24\2021-10-07_i24_AC_left_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\i24\2021-10-08_i24_AC_left_OK\recording2_task'

    'G:\Ishizu_data\IntermediateFiles\auditory\i34\2021-11-17_i34_AC_left_OK\recording2_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\i34\2021-11-18_i34_AC_left_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\i34\2021-11-19_i34_AC_left_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\i34\2021-12-08_i34_AC_Right_OK\recording2_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\i34\2021-12-09_i34_AC_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\i34\2021-12-10_i34_AC_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\i34\2021-12-11_i34_AC_Right_OK\recording1_task'

    'G:\Ishizu_data\IntermediateFiles\auditory\i35\2021-11-10_i35_AC_Left_OK\recording2_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\i35\2021-11-11_i35_AC_left_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\i35\2021-11-12_i35_AC_left_OK\recording3_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\i35\2021-11-13_i35_AC_Left_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\i35\2021-12-01_i35_AC_Right_OK\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\i35\2021-12-02_i35_AC_Right_OK\recording3_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\i35\2021-12-03_i35_AC_Right_OK\recording1_task'
    
    'G:\Ishizu_data\IntermediateFiles\auditory\i43\2022-06-11_i43_1AC2PPC_Right\recording2_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\i43\2022-06-15_i43_AC_Left\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\i43\2022-06-16_i43_1AC2PPC_Left\recording1_task'
    'G:\Ishizu_data\IntermediateFiles\auditory\i43\2022-06-17_i43_1AC2PPC_Left\recording1_task'
    
    'G:\Ishizu_data\IntermediateFiles\auditory\i46\2022-06-08_i46_AC_Left\recording1_task'
    };

%short_length = [zeros(5,1);zeros(1,1);zeros(4,1);ones(7,1);zeros(7,1)];
short_length = [zeros(5,1);zeros(1,1);zeros(4,1);ones(7,1);zeros(7,1);ones(4,1);zeros(1,1);];
if length(analysis_dir) ~= length(short_length)
    hoge
end
