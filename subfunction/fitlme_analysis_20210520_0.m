%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lme = fitlme_analysis_20210520_0(data,all_subject)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%make table
tbl = table(data,all_subject,'VariableNames',{'value','random1'});
lme0 = fitlme(tbl,'value ~ 1 + random1'); %linear regression
lme1 = fitlme(tbl,'value ~ 1 + (1|random1)'); %linear regression

% tbl
% data_name
% lme0 = fitlme(tbl,'value ~ fixed + random1'); %linear regression
% lme1 = fitlme(tbl,'value ~ fixed+(1|random1)');
% lme2 = fitlme(tbl,'value ~ fixed+(fixed|random1)');
% lme3 = fitlme(tbl,'value ~ fixed+(1|random1)+(fixed-1|random1)');
% %lme4 = fitlme(tbl,'value ~ 1 + fixed+(1|random1)'); % same as lme4

lme(1).lme = lme0;
lme(2).lme = lme1;

end