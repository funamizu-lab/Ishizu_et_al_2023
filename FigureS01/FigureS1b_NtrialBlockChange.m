function FigureS1b_NtrialBlockChange

parent='G:/Ishizu_data';
workpath='/Revise_ishizu';
outpath='/Revise_ishizu/output/NtrialBlockChange';

% hulistic parameter %
xTrial=[-5,8];% range of plotting / xTrial=0: Block change timing
useES =2;% 1:All evidence strength / 2:moderate+difficult / 3:difficult
figsaveTYPE='-dsvg';
%--------------------%

ESname ={'All ES', 'Mode Diff', 'Diff'};
ESrage ={[0.5 1.05], [0.35 1.05], [0.2 1.05]};

cd([parent,workpath]);

%%% collect the choice data %%%
behavePath='/Tokyo_ephys_ishizu/only_all_behaviors';
cd([parent,behavePath]);
list = dir('*_behave');

xTrial=[xTrial(1):-1,1:xTrial(end)];
Choice_BlockChange = cell(length(list),1);
mousename = cell(length(list),2);
for i=1:length(list)
    cd(list(i).name);
    mousename{i,1} = list(i).name(1:3);
    mousename{i,2} = ['mouse',num2str(i)];
    Choice_BlockChange{i} = getChoice_BlockChange(xTrial,useES);
    cd([parent,behavePath]);
end


%% save folder %%%
if(~exist([parent,outpath],'dir')), mkdir([parent,outpath]); end
cd([parent,outpath]);
save('savemat.mat','Choice_BlockChange','mousename','xTrial','useES');

%% plot figure %%
close all;
es = length(useES);
EachAnimalLabel  = cell(length(Choice_BlockChange),1);
EachAnimalSession= cell(length(Choice_BlockChange),1);
EachAnimalBlock  = cell(length(Choice_BlockChange),1);
EachAnimalChoice = cell(2,length(Choice_BlockChange),es);
EachAnimalChoice_long = cell(2,length(Choice_BlockChange),es);
EachAnimalChoice_short= cell(2,length(Choice_BlockChange),es);
for i=1:length(Choice_BlockChange)
    data= Choice_BlockChange{i};
    
    blockchange=data.blockChangeNum_EachSession;
    EachAnimalLabel{i}=i*ones(sum(blockchange),1);
    tmpEachAnimalSession= cell(length(blockchange),1);
    tmpEachAnimalBlock  = cell(length(blockchange),1);
    for k=1:length(blockchange)
        tmpEachAnimalSession{k}=k*ones(blockchange(k),1);
        tmpEachAnimalBlock{k}  =transpose(1:blockchange(k));
    end
    EachAnimalSession{i}= cell2mat(tmpEachAnimalSession);
    EachAnimalBlock{i}  = cell2mat(tmpEachAnimalBlock);
    
    for j=1:es
        EachAnimalChoice{1,i,j}= matTransform2(data.prefer,xTrial,j);
        EachAnimalChoice_long{1,i,j}= matTransform2(data.prefer_long,xTrial,j);
        EachAnimalChoice_short{1,i,j}=matTransform2(data.prefer_short,xTrial,j);
        EachAnimalChoice{2,i,j}= matTransform2(data.nprefer,xTrial,j);
        EachAnimalChoice_long{2,i,j}= matTransform2(data.nprefer_long,xTrial,j);
        EachAnimalChoice_short{2,i,j}=matTransform2(data.nprefer_short,xTrial,j);
    end
end
label_animal =cell2mat(EachAnimalLabel);
label_session=cell2mat(EachAnimalSession);
label_block  =cell2mat(EachAnimalBlock);

[~,blockchange] = min(abs(xTrial));
xtriallabel=cell(1,length(xTrial));
for i=1:length(xTrial)
    xtriallabel{i} = num2str(xTrial(i));
end
for j=1:es
    h=figure('Position',[20 20 1800 1600]);
    pref = cell2mat(EachAnimalChoice(1,:,j));
    npref= cell2mat(EachAnimalChoice(2,:,j));
    [estVal,Pval] = sigTest(pref,npref,xTrial,label_animal,label_session,label_block);
    [p_pref, err1_pref, err2_pref] = getplotdata(pref);
    [p_npref,err1_npref,err2_npref]= getplotdata(npref);
   
    subplot(3,2,1); hold on;
    p1=errorplot(1:length(xTrial),p_pref', err1_pref', err2_pref','r',0.5,1.5);
    p2=errorplot(1:length(xTrial),p_npref',err1_npref',err2_npref','b',0.5,1.5);
    plot([blockchange+0.5,blockchange+0.5],[-1,1],'--k');
    legend([p1,p2],{'prefer','nonprefer'},'Location','southwest');
    axis([0 length(xTrial)+1 ESrage{useES(j)}(1) ESrage{useES(j)}(2)]);
    xticks(1:length(xTrial)); xticklabels(xtriallabel);
    title('long and short tone');
    xlabel('trail');
    ylabel('correct rate'); 
    
    subplot(3,2,2); hold on;
    plot(1:length(xTrial),estVal,'g');
    plot([blockchange+0.5,blockchange+0.5],[-1,1],'--k');
    plot([0.5,length(xTrial)+0.5],[0,0],'--k');
    xticks(1:length(xTrial)); xticklabels(xtriallabel);
    title(['p val: ',num2str(Pval')]);
    xlabel('trail');
    ylabel('pref - npref');
    
    pref_long = cell2mat(EachAnimalChoice_long(1,:,j));
    npref_long= cell2mat(EachAnimalChoice_long(2,:,j));
    [estVal,Pval] = sigTest(pref_long,npref_long,xTrial,label_animal,label_session,label_block);
    [p_pref, err1_pref, err2_pref] = getplotdata(pref_long);
    [p_npref,err1_npref,err2_npref]= getplotdata(npref_long);
    
    subplot(3,2,3); hold on;
    p1=errorplot(1:length(xTrial),p_pref', err1_pref', err2_pref','r',0.5,1.5);
    p2=errorplot(1:length(xTrial),p_npref',err1_npref',err2_npref','b',0.5,1.5);
    plot([blockchange+0.5,blockchange+0.5],[-1,1],'--k');
    legend([p1,p2],{'prefer','nonprefer'},'Location','southwest');
    axis([0 length(xTrial)+1 ESrage{useES(j)}(1) ESrage{useES(j)}(2)]);
    xticks(1:length(xTrial)); xticklabels(xtriallabel);
    title('long tone');
    xlabel('trail');
    ylabel('correct rate');
    
    subplot(3,2,4); hold on;
    plot(1:length(xTrial),estVal,'g');
    plot([blockchange+0.5,blockchange+0.5],[-1,1],'--k');
    plot([0.5,length(xTrial)+0.5],[0,0],'--k');
    xticks(1:length(xTrial)); xticklabels(xtriallabel);
    title(['p val: ',num2str(Pval')]);
    xlabel('trail');
    ylabel('pref - npref');
    
    pref_short = cell2mat(EachAnimalChoice_short(1,:,j));
    npref_short= cell2mat(EachAnimalChoice_short(2,:,j));
    [estVal,Pval] = sigTest(pref_short,npref_short,xTrial,label_animal,label_session,label_block);
    [p_pref, err1_pref, err2_pref] = getplotdata(pref_short);
    [p_npref,err1_npref,err2_npref]= getplotdata(npref_short);
    
    subplot(3,2,5); hold on;
    p1=errorplot(1:length(xTrial),p_pref', err1_pref', err2_pref','r',0.5,1.5);
    p2=errorplot(1:length(xTrial),p_npref',err1_npref',err2_npref','b',0.5,1.5);
    plot([blockchange+0.5,blockchange+0.5],[-1,1],'--k');
    legend([p1,p2],{'prefer','nonprefer'},'Location','southwest');
    axis([0 length(xTrial)+1 ESrage{useES(j)}(1) ESrage{useES(j)}(2)]);
    xticks(1:length(xTrial)); xticklabels(xtriallabel);
    title(['short tone']);
    xlabel('trail');
    ylabel('correct rate');
    
    subplot(3,2,6); hold on;
    plot(1:length(xTrial),estVal,'g');
    plot([blockchange+0.5,blockchange+0.5],[-1,1],'--k');
    plot([0.5,length(xTrial)+0.5],[0,0],'--k');
    xticks(1:length(xTrial)); xticklabels(xtriallabel);
    title(['p val: ',num2str(Pval')]);
    xlabel('trail');
    ylabel('pref - npref');
    
    set(h,'PaperPositionMode','auto');
    print(h,'-r0',ESname{useES(j)},figsaveTYPE);
end
end

%-------------------------------------------------------------------------%
function [phat, err1, err2] = getplotdata(data)
%-------------------------------------------------------------------------%
x = sum(data,2,'omitnan');
n = size(data,2)*ones(size(data,1),1) - sum(isnan(data),2);
[phat,pci] = binofit(x,n);
err1 = pci(:,2)-phat;
err2 = phat-pci(:,1);
end

%-------------------------------------------------------------------------%
function [estVal,Pval] = sigTest(pref,npref,xTrial,label_animal,label_session,label_block)
%-------------------------------------------------------------------------%
tmp1 =pref(xTrial<0,:)-npref(xTrial<0,:);
tmp2 =npref(xTrial>0,:)-pref(xTrial>0,:);
data=[tmp1;tmp2];

estVal= zeros(length(xTrial),1);
Pval  = zeros(length(xTrial),1);
for i=1:length(xTrial)
    test = data(i,:)';
    tbl = table(test ,label_animal,label_session,label_block,'VariableNames',{'value','mouse','session','block'});
    %linear regression
    lme = fitlme(tbl,'value ~ 1 + (1|mouse)+(1|session)+(1|block)');
    estVal(i)= lme.Coefficients.Estimate;
    Pval(i)  = lme.Coefficients.pValue;
end

end

%-------------------------------------------------------------------------%
function output = getChoice_BlockChange(xTrial,useES)
%-------------------------------------------------------------------------%
list = dir('Bpod_*');
Choicedata=cell(length(list),1);
sessionNum = length(list);
blockChangeNum_EachSession = zeros(sessionNum,1);
for i=1:length(list)
    load(list(i).name,'EvidenceStrength','Outcome','Correct_side','Chosen_side','StimDuration','TrialBlock','BlockReward','Tone_cloud');
    tmpChoicedata=blockChangeSquence(EvidenceStrength,Tone_cloud,Outcome,...
        Correct_side,Chosen_side,StimDuration,TrialBlock,BlockReward,xTrial,useES);
    Choicedata{i}=tmpChoicedata;
    blockChangeNum_EachSession(i)=size(tmpChoicedata.preferside{1},2);
end
es=length(useES);

pref =cell(es,length(list));
pref_long =cell(es,length(list));
pref_short=cell(es,length(list));
nonpref =cell(es,length(list));
nonpref_long =cell(es,length(list));
nonpref_short=cell(es,length(list));
for j=1:es
    for i=1:length(list)
        pref{j,i} =Choicedata{i}.preferside{j};
        pref_long{j,i} =Choicedata{i}.preferside_long{j};
        pref_short{j,i}=Choicedata{i}.preferside_short{j};
        nonpref{j,i} =Choicedata{i}.nonpreferside{j};
        nonpref_long{j,i} =Choicedata{i}.nonpreferside_long{j};
        nonpref_short{j,i}=Choicedata{i}.nonpreferside_short{j};
    end
end

output.sessionNum = sessionNum;
output.blockChangeNum_EachSession = blockChangeNum_EachSession;
output.prefer = cell2mat(pref);
output.prefer_long =cell2mat(pref_long);
output.prefer_short=cell2mat(pref_short);
output.nprefer=cell2mat(nonpref);
output.nprefer_long =cell2mat(nonpref_long);
output.nprefer_short=cell2mat(nonpref_short);
end

%-------------------------------------------------------------------------%
function output = blockChangeSquence(EvidenceStrength,Tone_cloud,...
    Outcome,Correct_side,Chosen_side,StimDuration,TrialBlock,BlockReward,xTrial,useEStype)
%-------------------------------------------------------------------------%

%%% get tone evidence data in each trial %%%
temp_evi = unique(EvidenceStrength);
temp_evi_low  = 0.5 - temp_evi/2;
temp_evi_high = 0.5 + temp_evi/2;

% Put tone evidence in all trials;
trial_evidence = zeros(length(EvidenceStrength),1);
low  = find(Correct_side == 0);
high = find(Correct_side == 1);
for i = 1:length(temp_evi)
    temp = find(EvidenceStrength == temp_evi(i));
    temp_left  = intersect(temp,low);
    temp_right = intersect(temp,high);
    trial_evidence(temp_left)  = temp_evi_low(i);
    trial_evidence(temp_right) = temp_evi_high(i);
end

% Make the true tone cloud value
binary_tone = zeros(length(Tone_cloud),1);
for i = 1:length(Tone_cloud)
    temp_tone = Tone_cloud(i).matrix;
    %Get the data in all sound
    temp1 = find(temp_tone >= 9);
    binary_tone(i) = length(temp1) ./ length(temp_tone);
end
toneType=unique(binary_tone);
Easy=[toneType(1);toneType(end)];
Mode=[toneType(2);toneType(end-1)];
Diff=[toneType(3);toneType(end-2)];

useES_all ={[Easy; Mode;Diff], [Mode;Diff], Diff};
useES = useES_all(useEStype);

%Based on the correct trial, flip the tone cloud
if mean(binary_tone(Correct_side == 1)) < 0.5 %low for right correct(LowRight)
    binary_tone = 1 - binary_tone;% flipping tones 
else
end

%%% get choice data around block change %%%
Choice_trial = find(Outcome == 1 | Outcome == 2);
Outcome  = Outcome(Choice_trial);
Correct_side  = Correct_side(Choice_trial);
LowCorrect  = find(Correct_side ==0 & Outcome == 2);
LowIncorre  = find(Correct_side ==0 & Outcome == 1);
HighCorrect = find(Correct_side ==1 & Outcome == 2);
HighIncorre = find(Correct_side ==1 & Outcome == 1);
binary_tone  = binary_tone(Choice_trial);
StimDuration = StimDuration(Choice_trial);
TrialBlock2 = TrialBlock(Choice_trial);
BlockReward = BlockReward(:,2);
% Find block change
temp = TrialBlock2(2:length(TrialBlock2)) - TrialBlock2(1:length(TrialBlock2)-1);
Blockchange = find(temp ~= 0) + 0.5; %start of new block

%Make long/short performance
StimLength = unique(StimDuration);
StimCategory = ones(1,length(StimDuration));
StimCategory(StimDuration == max(StimLength)) = 2;
short_trial= find(StimCategory == 1);
long_trial = find(StimCategory == 2);

LowCorrect_long  = intersect(LowCorrect, long_trial);
LowIncorre_long  = intersect(LowIncorre, long_trial);
HighCorrect_long = intersect(HighCorrect,long_trial);
HighIncorre_long = intersect(HighIncorre,long_trial);
LowCorrect_short = intersect(LowCorrect, short_trial);
LowIncorre_short = intersect(LowIncorre, short_trial);
HighCorrect_short= intersect(HighCorrect,short_trial);
HighIncorre_short= intersect(HighIncorre,short_trial);

LowC={LowCorrect,LowCorrect_long,LowCorrect_short};
LowW={LowIncorre,LowIncorre_long,LowIncorre_short};
HighC={HighCorrect,HighCorrect_long,HighCorrect_short};
HighW={HighIncorre,HighIncorre_long,HighIncorre_short};

% correct transition around block change (-xx trial ~ end)
if length(Choice_trial)-Blockchange(end) < 20 % at least 30 trials requirement
    num=length(Blockchange)-1;
else
    num=length(Blockchange);
end

% block type
if BlockReward(2)>BlockReward(3)%R->L->R->L->...
    k=1;
elseif BlockReward(3)>BlockReward(2)%L->R->L->R->...
    k=2;
end

preferside =cell(num-1,length(useES),3); nonpreferside =cell(num-1,length(useES),3);
for es=1:length(useES)
    useEStrial = find(sum(binary_tone==useES{es}',2));
    for stim=1:3% whole/ long/ short
        low_c = intersect(useEStrial,LowC{stim});
        low_w = intersect(useEStrial,LowW{stim});
        high_c = intersect(useEStrial,HighC{stim});
        high_w = intersect(useEStrial,HighW{stim});
        
         for i = 2:num% eliminate 1st block change
             RightChoiceOutcome= getChoiceOutcome(low_c, low_w, Blockchange(i),xTrial);%1:correct / 0:incorrect / nan:no data
             LeftChoiceOutcome = getChoiceOutcome(high_c,high_w,Blockchange(i),xTrial);%1:correct / 0:incorrect / nan:no data
           if(mod(k+i,2)==0)
                preferside{i-1,es,stim}   = RightChoiceOutcome;
                nonpreferside{i-1,es,stim}= LeftChoiceOutcome;
            else
                preferside{i-1,es,stim}   = LeftChoiceOutcome;
                nonpreferside{i-1,es,stim}= RightChoiceOutcome;
            end
        end
    end
end

output.preferside =matTransform(preferside(:,:,1));
output.preferside_long =matTransform(preferside(:,:,2));
output.preferside_short=matTransform(preferside(:,:,3));
output.nonpreferside =matTransform(nonpreferside(:,:,1));
output.nonpreferside_long =matTransform(nonpreferside(:,:,2));
output.nonpreferside_short=matTransform(nonpreferside(:,:,3));
end

%-------------------------------------------------------------------------%
function ChoiceOutcome = getChoiceOutcome(cTrial,wTrial,BlockchangeTrial,xTrial)
%-------------------------------------------------------------------------%
% ChoiceOutcome %1:correct / 0:incorrect / nan:no data
tmp_c=cTrial-BlockchangeTrial;
tmp_w=wTrial-BlockchangeTrial;
[B,I]=sort([tmp_c;tmp_w]);

id_one = find(B>0,1);
choiceOutcome1=nan(xTrial(end),1);
for i=1:xTrial(end)
    tmp = id_one+i-1;
    if(ismember(tmp,1:length(I)))
        if(I(tmp)<=length(tmp_c))
            choiceOutcome1(i)=1;
        else
            choiceOutcome1(i)=0;
        end
    else
    end
end
            
id_minusone = id_one-1;
choiceOutcome2=nan(-xTrial(1),1);
for i=1:-xTrial(1)
    tmp = id_minusone-i+1;
    if(ismember(tmp,1:length(I)))
        if(I(tmp)<=length(tmp_c))
            choiceOutcome2(i)=1;
        else
            choiceOutcome2(i)=0;
        end
    else
    end
end
ChoiceOutcome = [flip(choiceOutcome2);choiceOutcome1];
end

%-------------------------------------------------------------------------%
function outcell = matTransform(mat)
%-------------------------------------------------------------------------%
outcell=cell(size(mat,2),1);
for i=1:size(mat,2)
    outcell{i}=cell2mat(mat(:,i)');
end
end

%-------------------------------------------------------------------------%
function out = matTransform2(mat,xTrial,es)
%-------------------------------------------------------------------------%
tmp = (1:length(xTrial))+length(xTrial)*(es-1);
out = mat(tmp,:);
end