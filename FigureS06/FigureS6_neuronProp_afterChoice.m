function FigureS6_neuronProp_afterChoice

parent='G:/Ishizu_data';
outPath='/Revise_ishizu/output/neuronProp_afterChoice2';

% hulistic parameter %
p_threshold = 1.0e-10; 
figsaveTYPE='-dsvg';
%--------------------%



%% plot figure %%
% save folder setting %
% region = {'auditory','fof','mpfc'};
folders= {'auc_ishizu','fof_ishizu','mpfc_ishizu'};
savefolder =[parent,outPath];
if(~exist(savefolder,'dir')), mkdir(savefolder); end
cd(savefolder);

close all
CompareSequence('dataMPFC_choice.mat','dataMPFC_sound.mat',p_threshold,folders{3},savefolder,'MPFC_pie','MPFC_psycho','MPFC_sig_time', figsaveTYPE);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PrefvsNpref,sdata]=drawData(flag,name,pref,npref,neuronID,ps,id,Ylim,figsaveTYPE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

evidence = [0,0.25,0.45,0.55,0.75,1];
if(flag==1)
    P = ps(neuronID,:);
elseif(flag==2)
    P = ones(length(neuronID),size(ps,2));
end
% ID1 = find(P(:,id(1))==1);

Pref  = zeros(length(id),length(evidence));
Pref_se = zeros(length(id),length(evidence));
Nonpref = zeros(length(id),length(evidence));
Nonpref_se = zeros(length(id),length(evidence));
PrefvsNpref = cell(length(id),1);
hval = zeros(length(id),length(evidence));
p_anova= zeros(length(id),2);
neuron=zeros(length(id),1);
for i=1:length(id)
    tmp_pref=zeros(length(find(P(:,id(i))==1)),length(evidence));
    tmp_npref=zeros(length(find(P(:,id(i))==1)),length(evidence));
    ano_activ = [];
    ano_sound = [];
    ano_block = [];
    for j=1:6
        tmp_p=pref(id(i),j).matrix(neuronID);
        tmp_n=npref(id(i),j).matrix(neuronID);
        p=tmp_p(P(:,id(i))==1);
        n=tmp_n(P(:,id(i))==1);
        tmp_pref(:,j)=p;
        tmp_npref(:,j)=n;
        nanid=union(find(isnan(p)==1),find(isnan(n)==1));
        p(nanid)=[];
        n(nanid)=[];
        neuron(i)=length(p);
        
        try
        [pval(i,j),hval(i,j)] = signrank(p,n,'alpha',0.01);
        catch
        end
        
        Pref(i,j)  = median(p);
        Nonpref(i,j) = median(n);
%         Pref(i,j)  = mean(p);
%         Nonpref(i,j) = mean(n);
        
        Pref_se(i,j)   = std(p)/sqrt(length(p));
        Nonpref_se(i,j)= std(n)/sqrt(length(n));
                
        sound0 = ones(length(p),1).*j;
        sound1 = ones(length(n),1).*j;
        block0 = zeros(length(p),1);
        block1 = ones(length(n),1);
        
        ano_activ = [ano_activ;p;n];
        ano_sound = [ano_sound;sound0;sound1];
        ano_block = [ano_block;block0;block1];
    end
    p = anovan(ano_activ,{ano_sound,ano_block},'display','off');
    p_anova(i,:) = p';
    PvsN=mean(tmp_pref,2)-mean(tmp_npref,2);
    PvsN(isnan(PvsN))=[];    
    PrefvsNpref{i} = PvsN;
end

h=figure('Position',[100 100 1500 500]);
sdata = struct();% source data 
for i=1:length(id)
    subplot(1,length(id),i);hold on
    e=errorbar(evidence,Pref(i,:),Pref_se(i,:),'r');
    e.CapSize = 0;
    e=errorbar(evidence,Nonpref(i,:),Nonpref_se(i,:),'b');
    e.CapSize = 0;
    xticks(0:0.2:1)
    xlim([-0.1, 1.1])
    ylim(Ylim)
    title({[num2str(id(i)),'  neuron:', num2str(neuron(i))],...
        ['p<0.01: ', num2str(hval(i,:))],...
        ['sound: ', num2str(p_anova(i,1))],...
        ['block: ', num2str(p_anova(i,2))]})
end
sdata.x= evidence';
sdata.preferrd=Pref(i,:)';
sdata.preferrd_se=Pref_se(i,:)';
sdata.nonpreferrd=Nonpref(i,:)';
sdata.nonpreferrd_se=Nonpref_se(i,:)';
sdata.pval=pval(i,:)';

set(h,'PaperPositionMode','auto');
print(h,'-r0',name,figsaveTYPE);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CompareSequence(choicemat,soundmat,p_threshold,folders,savefolder,figname1,figname2,figname3,figsaveTYPE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% get target neuron id %%%
[sig_longtone,sig_sustain,sig_choice,sig_prior,p_tone,p_choice] =...
    getPval_afterchoice(p_threshold,folders);
%
sig_choice_notone = setdiff(sig_choice,sig_longtone);
sig_lonchoi = intersect(sig_choice,sig_longtone);

cd(savefolder);


%%% FigS6a top left panel :  pie chart %%%
sig_prisus   = intersect(sig_prior,sig_sustain);
sig_pri_nosus= setdiff(sig_prior,sig_sustain);
sig_sus_nopri= setdiff(sig_sustain,sig_prior);
X=[length(sig_pri_nosus),length(sig_prisus),length(sig_sus_nopri)];

h=figure('Position',[100 100 1000 500]);
subplot(1,2,1)
pie(X,'%.3f%%')
title({'pure prior / pri nad sustain / pure sustain',num2str(X)})
subplot(1,2,2)
pie(X,{'pure prior','pri nad sustain','pure sustain'})
set(h,'PaperPositionMode','auto');
print(h,'-r0',figname1,figsaveTYPE);

%%% FigS6b top panels : get bias data %%%
load(soundmat,'evi_prefer_all','evi_nonprefer_all');
ID=15;
Ylim=[0.02 0.065];    [PrefvsNpref_tone,sdata1] = drawData(1,['toneS',figname2],evi_prefer_all,evi_nonprefer_all,sig_longtone,p_tone,ID,Ylim,figsaveTYPE);
Ylim=[-0.055 -0.035]; [PrefvsNpref_choi,sdata2] = drawData(2,['susS',figname2],evi_prefer_all,evi_nonprefer_all,sig_sustain,p_tone,ID,Ylim,figsaveTYPE);

% source data 
T = struct2table(sdata1);
writetable(T, 'source fig S6b top left.csv');
T = struct2table(sdata2);
writetable(T, 'source fig S6b top right.csv');


%%% FigS6c top panel %%%
h=figure('Position',[100 100 1500 500]);
tmp=[cell2mat(PrefvsNpref_tone);cell2mat(PrefvsNpref_choi)];
Ymin = min(tmp);Ymax = max(tmp);
% for i=1:length(ID)
for i=1
    X=[PrefvsNpref_tone{i};PrefvsNpref_choi{i}];
    tag=[repmat({'beforeSound'},length(PrefvsNpref_tone{i}),1);...
        repmat({'ITIsustained'},length(PrefvsNpref_choi{i}),1)];
    p=ranksum(PrefvsNpref_tone{i},PrefvsNpref_choi{i});
    subplot(1,length(ID),i);
    boxplot(X,tag);
    ylim([Ymin Ymax])
    title({num2str(ID(i)),['beforeSound: ', num2str(length(PrefvsNpref_tone{i})),...
        ' / ITIsustained: ', num2str(length(PrefvsNpref_choi{i}))],...
        ['pval: ',num2str(p)]})
end
set(h,'PaperPositionMode','auto');
print(h,'-r0',[figname2,'diffSound'],figsaveTYPE);

sdata = struct();% source data 
sdata.data = X;
sdata.tag = tag;
T = struct2table(sdata);
writetable(T, 'source fig S6c top.csv');


%%% FigS6b bottom panels : get bias data %%%
load(choicemat,'evi_prefer_all','evi_nonprefer_all');

ID=20;
Ylim=[0 0.07]; [PrefvsNpref_choi,sdata1] = drawData(2,['susC',figname2],evi_prefer_all,evi_nonprefer_all,sig_sustain,p_choice,ID,Ylim,figsaveTYPE);
Ylim=[-0.06 -0.02];  [PrefvsNpref_tone,sdata2] = drawData(2,['prior',figname2],evi_prefer_all,evi_nonprefer_all,sig_prior,p_choice,ID,Ylim,figsaveTYPE);

% source data 
T = struct2table(sdata1);
writetable(T, 'source fig S6b bottom right.csv');
T = struct2table(sdata2);
writetable(T, 'source fig S6b bottom left.csv');


%%% FigS6c bottom panel %%%
h=figure('Position',[100 100 1500 500]);
tmp=[cell2mat(PrefvsNpref_tone);cell2mat(PrefvsNpref_choi)];
Ymin = min(tmp);Ymax = max(tmp);
for i=1:length(ID)
    X=[PrefvsNpref_tone{i};PrefvsNpref_choi{i}];
    tag=[repmat({'tone'},length(PrefvsNpref_tone{i}),1);...
        repmat({'choi'},length(PrefvsNpref_choi{i}),1)];
    p=ranksum(PrefvsNpref_tone{i},PrefvsNpref_choi{i});
    subplot(1,length(ID),i);
    boxplot(X,tag);
    ylim([Ymin Ymax])
    title({num2str(ID(i)),['tone: ', num2str(length(PrefvsNpref_tone{i})),...
        ' / choi: ', num2str(length(PrefvsNpref_choi{i}))],...
        ['pval: ',num2str(p)]})
end
set(h,'PaperPositionMode','auto');
print(h,'-r0',[figname2,'diffChoice'],figsaveTYPE);

sdata = struct();% source data 
sdata.data = X;
sdata.tag = tag;
T = struct2table(sdata);
writetable(T, 'source fig S6c bottom.csv');


%%% FigS6a %%%
load('spout.mat');
neuroact_sus = zeros(length(sig_sustain),size(evi_prefer_all,1));
neuroact_pri = zeros(length(sig_pri_nosus),size(evi_prefer_all,1));
for i=1:size(evi_prefer_all,1)
    tmp_sus =zeros(length(sig_sustain),12);
    tmp_pri =zeros(length(sig_pri_nosus),12);
    for j=1:6
        tmp_sus(:,j)=evi_prefer_all(i).matrix(sig_sustain);
        tmp_pri(:,j)=evi_prefer_all(i).matrix(sig_pri_nosus);
        tmp_sus(:,j+6)=evi_nonprefer_all(i).matrix(sig_sustain);
        tmp_pri(:,j+6)=evi_nonprefer_all(i).matrix(sig_pri_nosus);
    end
    neuroact_sus(:,i) = mean(tmp_sus,2,'omitnan');
    neuroact_pri(:,i) = mean(tmp_pri,2,'omitnan');
end

plot_sus = zeros(2,size(evi_prefer_all,1));
plot_pri = zeros(2,size(evi_prefer_all,1));
pval=zeros(1,size(evi_prefer_all,1));
for i=1:size(evi_prefer_all,1)
    plot_sus(1,i) = mean(neuroact_sus(:,i));
    plot_sus(2,i) = std(neuroact_sus(:,i))/sqrt(length(sig_sustain));
    plot_pri(1,i) = mean(neuroact_pri(:,i));
    plot_pri(2,i) = std(neuroact_pri(:,i))/sqrt(length(sig_pri_nosus));
    pval(i) = ranksum(neuroact_pri(:,i),neuroact_sus(:,i));
end

[~,maxOrder]=max(neuroact_sus,[],2);
[~,sortOrder]=sort(maxOrder);
X = neuroact_sus(sortOrder,:);
neuroact_sus2=zeros(size(X,1),size(X,2));
for i=1:size(X,1)
    neuroact_sus2(i,:)=rescale(X(i,:));
end
[~,maxOrder]=max(neuroact_pri,[],2);
[~,sortOrder]=sort(maxOrder);
X= neuroact_pri(sortOrder,:);
neuroact_pri2=zeros(size(X,1),size(X,2));
for i=1:size(X,1)
    neuroact_pri2(i,:)=rescale(X(i,:));
end

h=figure('Position',[100 100 1300 1100]);
subplot(2,2,1);
hold on
p1=errorplot(1:50,plot_sus(1,:),plot_sus(2,:),plot_sus(2,:),'b',0.5,1);
p2=errorplot(1:50,plot_pri(1,:),plot_pri(2,:),plot_pri(2,:),'r',0.5,1);
legend([p1,p2],{'sustain','pri no sus'})
xlim([0 50])
xticks([(500+ave_spout(1))/100,5,20,(500+ave_spout(2))/100,45]);
xticklabels({'s','0','1.5','s','4'});

%%% source data %%%
sdata = struct();% source data 
sdata.xtime_sec = ((1:50)'-5)/10;
sdata.ITIsustained = plot_sus(1,:)';
sdata.ITIsustained_se = plot_sus(2,:)';
sdata.BeforeSound = plot_pri(1,:)';
sdata.BeforeSound = plot_pri(2,:)';
sdata.pval = pval';
T = struct2table(sdata);
writetable(T, 'source fig S6a bottom left.csv');


subplot(2,2,3);
hold on
plot(1:50,-log10(pval))
plot([0,51],[10, 10])
xticks([(500+ave_spout(1))/100,5,20,(500+ave_spout(2))/100,45]);
xticklabels({'s','0','1.5','s','4'});
xlim([0 50])

ax=subplot(2,2,2);
imagesc(neuroact_sus2);
colorbar
title('sustain');
xticks([(500+ave_spout(1))/100,5,(500+ave_spout(2))/100,45]);
xticklabels({'s','0','s','4'});
yticks([1,size(neuroact_sus2,1)]);
ax.YDir='reverse';

%%% source data %%%
sdata = struct();% source data 
sdata.ITIsustained_time = neuroact_sus2;
T = struct2table(sdata);
writetable(T, 'source fig S6a heatmap bottom right.csv');

ax=subplot(2,2,4);
imagesc(neuroact_pri2);
colorbar
title('prior');
xticks([(500+ave_spout(1))/100,5,(500+ave_spout(2))/100,45]);
xticklabels({'s','0','s','4'});
yticks([1,size(neuroact_pri2,1)]);
ax.YDir='reverse';

%%% source data %%%
sdata = struct();% source data 
sdata.BeforeSound_time = neuroact_pri2;
T = struct2table(sdata);
writetable(T, 'source fig S6a heatmap top right.csv');

set(h,'PaperPositionMode','auto');
print(h,'-r0',[figname3,'_neuralact'],figsaveTYPE);

end
