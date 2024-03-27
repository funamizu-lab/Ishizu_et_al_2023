%{
----------------------------------------------------------------------------
%Stimulus bias changes by block
%Just model analysis
%Based on Rao 2010 Front Comp Neurosci
%Value updating + Prior updating
----------------------------------------------------------------------------
%}

function Dual_RL_model_block1_20220818_para_search_all
clear all

% [filename1, pathname1,findex]=uigetfile('*.mat','Block_mat','Multiselect','on');
% filename1
filename1 = dir('Bpod_mat*');

for filecount = 1 : length(filename1)
%     temp_filename = filename1(filecount) 
%     temp_filename = cell2mat(temp_filename);
%     temp_path = pathname1;
%     fpath = fullfile(temp_path, temp_filename);
    temp_filename = filename1(filecount).name
    fpath = temp_filename;
    
%     [ave_likeli(filecount,:),BIC(filecount,:),log_likeli(filecount,:),...
%         para_max(filecount).matrix,N_trial(filecount,1)] = Dual_RL_model_block1_20220118_para_search(fpath);
    [ave_likeli,BIC,log_likeli,para_max,N_trial] = Dual_RL_model_block1_20220818_para_search(fpath);
    
    save_file = ['RL_20220818_confidence2_',temp_filename];
    save(save_file, 'ave_likeli' ,'BIC','log_likeli','para_max','N_trial')
end
ave_likeli

delete(gcp('nocreate'))