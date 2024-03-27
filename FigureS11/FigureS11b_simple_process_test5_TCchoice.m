
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
function FigureS11b_simple_process_test5_TCchoice(folders)

close all
analysis_dir = eval(folders);

for i = 1:40
    for j = 1:6
        evi_trace_all(i,j).matrix = [];
        evi_trace_e_all(i,j).matrix = [];
    end
end
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]
    
    [evi_trace, evi_trace_e] = ...
        Task_kaiseki_tokyo1_20220516_process4_ishizu_TCchoice(analysis_dir{i});
    
    for j = 1:40
        for k = 1:6
            evi_trace_all(j,k).matrix = [evi_trace_all(j,k).matrix; evi_trace(j,k).matrix];
            evi_trace_e_all(j,k).matrix = [evi_trace_e_all(j,k).matrix; evi_trace_e(j,k).matrix];
        end
    end
end


%%% FigS10b %%%
ID=25;
Ylim=[0. 0.12]; %auc
drawData(evi_trace_all,evi_trace_e_all,ID,Ylim);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawData(pref,npref,id,Ylim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

evidence = [0,0.25,0.45,0.55,0.75,1];
% pval = p_choice(neuronID,:);

Pref  = zeros(length(id),length(evidence));
Pref_se = zeros(length(id),length(evidence));
Nonpref = zeros(length(id),length(evidence));
Nonpref_se = zeros(length(id),length(evidence));
pval=zeros(length(id),length(evidence));
neuron=zeros(length(id),1);
neuron2=zeros(length(id),1);
for i=1:length(id)
    for j=1:6
        p=pref(id(i),j).matrix;
        n=npref(id(i),j).matrix;
        neuron(i)=length(p);
        nanid=union(find(isnan(p)==1),find(isnan(n)==1));
        p(nanid)=[];
        n(nanid)=[];
        neuron2(i)=length(p);
        
        pval(i,j)=signrank(p,n);
        
        Pref(i,j)  = median(p);
        Nonpref(i,j) = median(n);
        
        Pref_se(i,j)   = std(p)/sqrt(length(p));
        Nonpref_se(i,j)= std(n)/sqrt(length(n));
    end
end

%%% fig S11b %%%    
figure('Position',[100 100 1500 500]);
for i=1:length(id)
    subplot(1,length(id),i);hold on
    errorplot(evidence,Pref(i,:),Pref_se(i,:),Pref_se(i,:),'r',0.5,1);
    errorplot(evidence,Nonpref(i,:),Nonpref_se(i,:),Nonpref_se(i,:),'b',0.5,1);
    xticks(0:0.2:1)
    xlim([-0.1, 1.1])
    ylim(Ylim)
    title({[num2str(id(i)),'  neuron:', num2str(neuron(i)),'  for calc:', num2str(neuron2(i))],...
        num2str(pval(i,1:3)),num2str(pval(i,4:6))})
end

%%% source data %%%
cd('G:\upload_code\FigureS11');
sdata = struct();% source data 
sdata.x=evidence';
sdata.correct=Pref(i,:)';
sdata.correct_se=Pref_se(i,:)';
sdata.error=Nonpref(i,:)';
sdata.error_se=Nonpref_se(i,:)';
sdata.pval = pval';
T = struct2table(sdata);
writetable(T, 'source fig S11b.csv');

end