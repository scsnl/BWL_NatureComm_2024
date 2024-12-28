%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%------Dopaminergic modulation and dosage effects on brain state dynamics and-------------------------------------------%%%
%%%------working memory component processes in Parkinson's disease, Lee et al, Nature Communications, 2025----%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

%add relevant paths
addpath(genpath(sprintf('%s/functions', pwd)));
dataDir=sprintf('%s/data', pwd);

%% Load relevant files (subject list, timeseries data, BSDS model, covariance matrix, etc)
groupName = 'HCPDcombined'; %'HC', 'PDOFF', 'PDON' 
runName= 'AllRuns'; %'runA', 'runB', 'runC', 'runD', 'AllRuns'
runList={'runA', 'runB', 'runC', 'runD'};
groupList={'HC', 'PDOFF', 'PDON'}; 
taskList={'LL', 'HL', 'DL'}; 
roiSet = '11ROI_both_hemis';
ROI_names = {'lPPC','rPPC','lDLPFC','rDLPFC','DMPFC','lAI','rAI','lCau','rCau','lThal','rThal'};

TR=0.49; %sampling rate
num.ROI = length(ROI_names); % Number of ROIs
num.Task=length(taskList); % Number of task: LL, HL, DL
num.Run=length(runList); % Number of run: runA, runB, runC, runD
num.Vol=778; % length of timeseries
num.Subj=91; %37HC, 27PDOFF, 27

%Load subject information
DM_NPscore_dir=fullfile(dataDir, 'Demographics_n_NPscores', 'HCPD_Subject_Information.csv');
SubjInfo=readtable(DM_NPscore_dir); clear DM_NPscore_dir;

% Load Sternberg task design
taskdesgn_dir = fullfile(dataDir, 'taskdesign_n_taskScores', 'HCPD_taskdesign_n_taskscores.mat');  
load(taskdesgn_dir); clear taskdesgn_dir;

% Load BSDS model (Time varying posterior probabilities of each brain state)
model_dir = fullfile(dataDir, 'model', 'HCPD_BSDS_model.mat');
load(model_dir); clear model_dir;

%Load Covariance matrix estimted from BSDS
Cov_dir=fullfile(dataDir, 'Covariance_Matrix', 'HCPD_covMtx.mat');
load(Cov_dir); clear Cov_dir;


%% Reordering brain states (presentation purpose)
num.State=length(dominantStates);
StateReorder=[4,1,5,2,3];

%% Group index
HC_idx=find(strcmp(SubjInfo.Group, 'HC')==1);
PDOFF_idx=find(strcmp(SubjInfo.Group, 'PDOFF')==1);
PDON_idx=find(strcmp(SubjInfo.Group, 'PDON')==1);
PDOFF_nonMCI_idx=intersect(find(strcmp(SubjInfo.Group, 'PDOFF')==1), find(SubjInfo.CogImpair==0));
PDON_nonMCI_idx=intersect(find(strcmp(SubjInfo.Group, 'PDON')==1), find(SubjInfo.CogImpair==0));


%% extract task onset and temporal evolution of latent brain states in each phase of Sternberg WM from the design matrix
% Encoding phase (2s:stim_onset~maintainence_onsest)
% maintainence phase (4~8s: maintainence_onset~response onset)
% Retreival phase (1s: response onset~1s after response onset)
[ITIonsetNoffset, TaskonsetNoffset] = func_taskblockInfo(task_design_info, num, TR);


%% Supplementary Figure 1: Medication information of the PD participants
%Disease duration of the PD patients
DisDur=SubjInfo.DiseaseDur_y(PDOFF_idx);
% MOCA score
Moca_PD=SubjInfo.MOCA(PDOFF_idx);
% Depression score of PD
GDS_PD=SubjInfo.GDS(PDOFF_idx);
GDS_range = [sum(GDS_PD >= 0 & GDS_PD <= 4), sum(GDS_PD >= 5 & GDS_PD <= 8),... % Normal_range, mild_range
                          sum(GDS_PD >= 9 & GDS_PD <= 11), sum(GDS_PD >= 12 & GDS_PD <= 15)]; %Moderate_range, Severe_range

figure; subplot(1,3,1);
histogram(DisDur);  ylim([0 7]);
xlabel('Disease duration (years)'); ylabel('Number of patients'); title('Disease duration of the patients')
subplot(1,3,2);
histogram(Moca_PD); ylim([0 7]);
xlabel('MOCA'); ylabel('Number of patients'); title('MOCA score of the patients');
subplot(1,3,3);
bar(GDS_range);
% Customize the x-axis labels
set(gca, 'XTickLabel', {'0-4', '5-8', '9-11', '12-15'});
% Add labels and title
xlabel('GDS score'); ylabel('Number of patients'); title('GDS score of the patients');
sgtitle('Medication Information of the PD participants')


%% Figure 2B: Sternberg working memory task performance
% transparent boxplots with no whiskers and jittered datapoints underneath
group_names={'HC', 'PDOFF', 'PDON'};
condition_names={'LL', 'HL', 'DL'};

group_idx = [ones(1,length(HC_idx)), 2.*ones(1,length(PDOFF_idx)), 3.*ones(1,length(PDON_idx))];
Accuracy = [SubjInfo.Sternberg_LL_ACC, SubjInfo.Sternberg_HL_ACC, SubjInfo.Sternberg_DL_ACC];
RT = [SubjInfo.Sternberg_LL_RT, SubjInfo.Sternberg_HL_RT, SubjInfo.Sternberg_DL_RT];

figure;
subplot(1,2,1);
plot_bargraph([Accuracy.*100], group_idx, condition_names, group_names, 'Accuracy(%)', [70, 110], 11, 'Accuracy' );
subplot(1,2,2);
plot_bargraph(RT, group_idx, condition_names, group_names, 'Time(Sec)', [0.7, 2.4], 11, 'Reaction Time' );
sgtitle('Sternberg Working Memory Task Performance');

% % ANOVA analysis
% %%----HC vs PDOFF----%%
% Behavior='RT'; % 'ACC', 'RT
% Group = []; % 1:HC, 2:PDOFF
% Load = [];  % Task type
% % Prepare data for ANOVA
% Behav_Anova = [];
% for task = 1:3
%     if strcmp(Behavior, 'ACC')
%         HC_data = Accuracy(HC_idx, task);
%         PDOFF_data = Accuracy(PDOFF_idx, task);
%     elseif strcmp(Behavior, 'RT')
%         HC_data = RT(HC_idx, task);
%         PDOFF_data =RT(PDOFF_idx, task);
%     end
%     Behav_Anova = [Behav_Anova, HC_data', PDOFF_data']; % Append to OR_Anova
%     Group = [Group, ones(1, length(HC_data)), 2 * ones(1, length(PDOFF_data))];
%     Load = [Load, task * ones(1, length(HC_data) + length(PDOFF_data))];
% end
% % Perform ANOVA
% p = anovan(Behav_Anova, {Group, Load}, 'model', 2, 'varnames', {'Group', 'Load'});
% 
% posthoc=[];
% [~, posthoc.p(1,1), ~, posthoc.stats{1,1}] = ttest(Behav_Anova(find(Load==1)), Behav_Anova(find(Load==2))); % LL vs HL
% [~, posthoc.p(2,1), ~, posthoc.stats{1,2}] = ttest(Behav_Anova(find(Load==1)), Behav_Anova(find(Load==3))); % LL vs DL
% [~, posthoc.p(3,1), ~, posthoc.stats{1,3}] = ttest(Behav_Anova(find(Load==2)), Behav_Anova(find(Load==3))); % HL vs DL

% %%----PDOFF vs PDON----%%
% Behavior='RT'; % 'ACC', 'RT
% F=[]; p=[]; stats=[];
% Behav_Anova = [];    S = [];    Group = [];    Load = [];    count = 1;
% % Iterate through tasks and subjects
% for task = 1:3
%     for group = 1:2
%         if group == 1
%             idx = PDOFF_idx;  group_label = 1; % PDOFF
%         else
%             idx = PDON_idx;  group_label = 2; % PDON
%         end
%         for i = 1:length(idx)
%             if strcmp(Behavior, 'ACC')
%                 Behav_Anova(count) = Accuracy(idx(i),task);
%             elseif strcmp(Behavior,'RT')
%                 Behav_Anova(count) = RT(idx(i), task);
%             end
%             S(count) = i;
%             Group(count) = group_label;
%             Load(count) = task;
%             count = count + 1;
%         end
%     end
% end
% % Perform repeated measures ANOVA
% stats = rm_anova2(Behav_Anova, S, Group, Load, {'Group', 'Load'});
% 
% posthoc=[];
% [~, posthoc.p(1,1), ~, posthoc.stats{1,1}] = ttest(Behav_Anova(find(Load==1)), Behav_Anova(find(Load==2))); % LL vs HL
% [~, posthoc.p(2,1), ~, posthoc.stats{1,2}] = ttest(Behav_Anova(find(Load==1)), Behav_Anova(find(Load==3))); % LL vs DL
% [~, posthoc.p(3,1), ~, posthoc.stats{1,3}] = ttest(Behav_Anova(find(Load==2)), Behav_Anova(find(Load==3))); % HL vs DL


%% Get posterior propability changes of of each state across entire 4 runs in each group
% Figure 3A: Brain state dynamics during Sternberg WM task (run A)
% Fig.3A middle: Posterior probability of latent brain states
figure;
timepoints=0:TR:num.Vol;
subplot(2,1,1);
for run = 1%:numRun
    meanStatePP=squeeze((mean(grpStatePP(run:4:360+run,:,:))));
    stdStatePP=squeeze((std(grpStatePP(run:4:360+run,:,:))));
    colorstrings = {'#515564','#C25B56','#BEB9B5', '#74828F', '#96C0CE' };
    xpts = 0:length(meanStatePP)-1;     xrange = [0, (length(meanStatePP)-1)];     yrange = [0, 1];
    linewidth = 2;
    for i = 1:num.State
        hold on; plot(xpts, meanStatePP(:,i), 'Color', colorstrings{i}, 'LineWidth', linewidth);
    end
    ylim(yrange); xlim(xrange); xlabel('Time(s)', 'FontSize', 12); ylabel('posterior probability', 'FontSize', 12);
    title(sprintf('Posterior probability of latent brain states in Run: %d', run));
end

% Fig. 3A Bottom: Plotting Temporal evolution of latent brain states
Indx=cell(1,num.State);
for i=1:num.State
    Indx{1,i}=find(grpTempEvol==dominantStates(i));
end
grpTempEvol_re=zeros(size(grpTempEvol));
for i=1:num.State
    grpTempEvol_re(Indx{1,i})=i;
end
subplot(2,1,2);
for run =1%:num.Run
    hAxes = gca; imagesc(hAxes, grpTempEvol_re(run:4:360+run,:));
    colormap( hAxes , [0.32 0.33 0.39; 0.76 0.36 0.34; 0.75 0.73 0.71; 0.45 0.51 0.56; 0.59 0.75 0.81]);
    xlabel('Time(s)', 'FontSize', 12); ylabel('Participants', 'FontSize', 12);
    title(sprintf('Temporal Evolution of latent brain states in Run: %d', run));
end

% Figure 3B-D
[StatePP_GroupStat, StatePP_SubjStat, StatePP_everyTrials]=func_statePPproperties(TaskonsetNoffset, grpStatePP, num);
plot_timeVaryingPosteriorProb(StatePP_GroupStat, StatePP_SubjStat, num, groupList, HC_idx, 1)
plot_timeVaryingPosteriorProb(StatePP_GroupStat, StatePP_SubjStat, num, groupList, PDOFF_idx, 2)
plot_timeVaryingPosteriorProb(StatePP_GroupStat, StatePP_SubjStat, num, groupList, PDON_idx, 3)

HC_PP=[];
PDOFF_PP=[];
PDON_PP=[];
for task = 1:num.Task
	for state = 1:num.State
            HC_PP{1,task}(state,:)=mean(StatePP_SubjStat.Avg.wholeTB{1,task}(HC_idx,:,StateReorder(state)));
            PDOFF_PP{1,task}(state,:)=mean(StatePP_SubjStat.Avg.wholeTB{1,task}(PDOFF_idx,:,StateReorder(state)));
            PDON_PP{1,task}(state,:)=mean(StatePP_SubjStat.Avg.wholeTB{1,task}(PDON_idx,:,StateReorder(state)));
    end
end

% Figure S3A-C
[StatePP_GroupStat_ITI, StatePP_SubjStat_ITI, StatePP_everyTrials_ITI]=func_statePPproperties_includingITI(TaskonsetNoffset, grpStatePP, num);
plot_timeVaryingPosteriorProb(StatePP_GroupStat_ITI, StatePP_SubjStat_ITI, num, groupList, HC_idx, 1)
plot_timeVaryingPosteriorProb(StatePP_GroupStat_ITI, StatePP_SubjStat_ITI, num, groupList, PDOFF_idx, 2)
plot_timeVaryingPosteriorProb(StatePP_GroupStat_ITI, StatePP_SubjStat_ITI, num, groupList, PDON_idx, 3)

HC_PP_ITI=[];
PDOFF_PP_ITI=[];
PDON_PP_ITI=[];
for task = 1:num.Task
	for state = 1:num.State
            HC_PP_ITI{1,task}(state,:)=mean(StatePP_SubjStat_ITI.Avg.wholeTB{1,task}(HC_idx,:,StateReorder(state)));
            PDOFF_PP_ITI{1,task}(state,:)=mean(StatePP_SubjStat_ITI.Avg.wholeTB{1,task}(PDOFF_idx,:,StateReorder(state)));
            PDON_PP_ITI{1,task}(state,:)=mean(StatePP_SubjStat_ITI.Avg.wholeTB{1,task}(PDON_idx,:,StateReorder(state)));
    end
end


%% Temporal dynamics of latent brain states during entire task block(E+M+R) in the group --------%%%
targetPhase = 'whole'; %'whole', 'encoding', 'maintenance', 'retrieval'
[~, TemEvol_SubjStat_W]=func_tempEvolProperties(TaskonsetNoffset, grpTempEvol, dominantStates, targetPhase, num);

targetPhase = 'encoding'; %'whole', 'encoding', 'maintenance', 'retrieval'
[~, TemEvol_SubjStat_E]=func_tempEvolProperties(TaskonsetNoffset, grpTempEvol, dominantStates, targetPhase, num);

targetPhase = 'M&R'; %'whole', 'encoding', 'maintenance', 'retrieval'
[~, TemEvol_SubjStat_MR]=func_tempEvolProperties(TaskonsetNoffset, grpTempEvol, dominantStates, targetPhase, num);

%% ANOVA analysis
% %%----HC vs PDOFF----%%
% phase='E'; %'E', 'MR'
% for state = 1: 5 % Specific state to analyze
% Group = []; % 1:HC, 2:PDOFF
% Load = [];  % Task type
% % Prepare data for ANOVA
% Behav_Anova = [];
% for task = 1:3
%     if strcmp(phase, 'E')
%         HC_data = TemEvol_SubjStat_E.OR{1, task}(HC_idx, StateReorder(state));
%         PDOFF_data = TemEvol_SubjStat_E.OR{1, task}(PDOFF_idx, StateReorder(state));    
%     elseif strcmp(phase,'MR')
%         HC_data = TemEvol_SubjStat_MR.OR{1, task}(HC_idx, StateReorder(state));
%         PDOFF_data = TemEvol_SubjStat_MR.OR{1, task}(PDOFF_idx, StateReorder(state));    
%     end
%     Behav_Anova = [Behav_Anova, HC_data', PDOFF_data']; % Append to OR_Anova
%     Group = [Group, ones(1, length(HC_data)), 2 * ones(1, length(PDOFF_data))];
%     Load = [Load, task * ones(1, length(HC_data) + length(PDOFF_data))];
% end
% % Perform ANOVA
% p(:, state) = anovan(Behav_Anova, {Group, Load}, 'model', 2, 'varnames', {'Group', 'Load'});
% end
% 
% %%----PDOFF vs PDON----%%
% phase='E'; ; %'E', 'MR'
% F=[]; p=[]; stats=[];
% for state = 1:5
%     Behav_Anova = [];    S = [];    Group = [];    Load = [];    count = 1;
%     % Iterate through tasks and subjects
%     for task = 1:3
%         for group = 1:2
%             if group == 1
%                 idx = PDOFF_idx;  group_label = 1; % PDOFF
%             else
%                 idx = PDON_idx;  group_label = 2; % PDON
%             end            
%             for i = 1:length(idx)
%                 if strcmp(phase, 'E')
%                     Behav_Anova(count) = TemEvol_SubjStat_E.OR{1, task}(idx(i), StateReorder(state));
%                 elseif strcmp(phase,'MR')
%                     Behav_Anova(count) = TemEvol_SubjStat_MR.OR{1, task}(idx(i), StateReorder(state));
%                 end
%                 S(count) = i;
%                 Group(count) = group_label;
%                 Load(count) = task;
%                 count = count + 1;
%             end
%         end
%     end
%     % Perform repeated measures ANOVA
%     stats{1, state} = rm_anova2(Behav_Anova, S, Group, Load, {'Group', 'Load'});
% end


%%  Plot state switching probability for each group: FIgure S4
StatTransProb.HC = plot_stateSwitchingProb(TemEvol_SubjStat_W, taskList, HC_idx, num, 'HC', StateReorder);
StatTransProb.PDOFF = plot_stateSwitchingProb(TemEvol_SubjStat_W, taskList, PDOFF_idx, num, 'PDOFF', StateReorder);
StatTransProb.PDON = plot_stateSwitchingProb(TemEvol_SubjStat_W, taskList, PDON_idx, num, 'PDON', StateReorder);


%% Plot occupancy rate of latent brain states: Figure S5
group_names={'S1', 'S2', 'S3', 'S4', 'S5'};
condition_names={'LL', 'HL', 'DL'};

group_HC_idx = repelem(1:num.State, length(HC_idx));
group_PD_idx = repelem(1:num.State, length(PDOFF_idx));

% Encoding phase
for task=1:3
    OR_E.HC{1,task}=TemEvol_SubjStat_E.OR{1,task}(HC_idx, StateReorder);
    OR_E.PDOFF{1,task}=TemEvol_SubjStat_E.OR{1,task}(PDOFF_idx, StateReorder);
    OR_E.PDON{1,task}=TemEvol_SubjStat_E.OR{1,task}(PDON_idx, StateReorder);
end

% Plot occupancy rates during Encoding phase
OR_HC_E = cell2mat(arrayfun(@(x) reshape(TemEvol_SubjStat_E.OR{1, x}(HC_idx, StateReorder), [], 1), 1:3, 'UniformOutput', false));
OR_PDOFF_E = cell2mat(arrayfun(@(x) reshape(TemEvol_SubjStat_E.OR{1, x}(PDOFF_idx, StateReorder), [], 1), 1:3, 'UniformOutput', false));
OR_PDON_E = cell2mat(arrayfun(@(x) reshape(TemEvol_SubjStat_E.OR{1, x}(PDON_idx, StateReorder), [], 1), 1:3, 'UniformOutput', false));

figure;
subplot(1,3,1);
plot_bargraph([OR_HC_E.*100], group_HC_idx, condition_names, group_names, 'Occupancy rate(%)', [0, 100], 11, 'HC: Encoding' );
subplot(1,3,2);
plot_bargraph([OR_PDOFF_E.*100], group_PD_idx, condition_names, group_names, 'Occupancy rate(%)', [0, 100], 11, 'PDOFF: Encoding' );
subplot(1,3,3);
plot_bargraph([OR_PDON_E.*100], group_PD_idx, condition_names, group_names, 'Occupancy rate(%)', [0, 100], 11, 'PDON: Encoding' );
sgtitle('Occupancy rate of latent brain states during Encoding phase');

% stats
targetOR_ratio=OR_E.HC;
h=[]; p=[];
for task=1:3
    meanOR_E{1,task}=mean(targetOR_ratio{1,task});
    for state_i=1:5
        for state_j=1:5
            [h{1,task}(state_i, state_j), p{1,task}(state_i,state_j)]=ttest(targetOR_ratio{1,task}(:,state_i), targetOR_ratio{1,task}(:,state_j));                
        end
    end
end

targetOR_ratio=OR_E.PDON;
h=[]; p=[];
for state =1:5
    for task_i=1:3
        for task_j=1:3
            [h{1,state}(task_i, task_j), p{1,state}(task_i,task_j)]=ttest(targetOR_ratio{1,task_i}(:,state), targetOR_ratio{1,task_j}(:,state));                
        end
    end
end


% Maintenance & Retrieval phases
for task=1:3
    OR_MR.HC{1,task}=TemEvol_SubjStat_MR.OR{1,task}(HC_idx, StateReorder);
    OR_MR.PDOFF{1,task}=TemEvol_SubjStat_MR.OR{1,task}(PDOFF_idx, StateReorder);
    OR_MR.PDON{1,task}=TemEvol_SubjStat_MR.OR{1,task}(PDON_idx, StateReorder);
end

OR_HC_MR = cell2mat(arrayfun(@(x) reshape(TemEvol_SubjStat_MR.OR{1, x}(HC_idx, StateReorder), [], 1), 1:3, 'UniformOutput', false));
OR_PDOFF_MR = cell2mat(arrayfun(@(x) reshape(TemEvol_SubjStat_MR.OR{1, x}(PDOFF_idx, StateReorder), [], 1), 1:3, 'UniformOutput', false));
OR_PDON_MR = cell2mat(arrayfun(@(x) reshape(TemEvol_SubjStat_MR.OR{1, x}(PDON_idx, StateReorder), [], 1), 1:3, 'UniformOutput', false));

figure;
subplot(1,3,1);
plot_bargraph([OR_HC_MR.*100], group_HC_idx, condition_names, group_names, 'Occupancy rate(%)', [0, 100], 11, 'HC: Maintenance&Retreival' );
subplot(1,3,2);
plot_bargraph([OR_PDOFF_MR.*100], group_PD_idx, condition_names, group_names, 'Occupancy rate(%)', [0, 100], 11, 'PDOFF: Maintenance&Retreival' );
subplot(1,3,3);
plot_bargraph([OR_PDON_MR.*100], group_PD_idx, condition_names, group_names, 'Occupancy rate(%)', [0, 100], 11, 'PDON: Maintenance&Retreival' );
sgtitle('Occupancy rate of latent brain states during maintenance&retrieval phases');

% stats
targetOR_ratio=OR_MR.PDON;
h=[]; p=[];
for task=1:3
    meanOR_E{1,task}=mean(targetOR_ratio{1,task});
    for state_i=1:5
        for state_j=1:5
            [h{1,task}(state_i, state_j), p{1,task}(state_i,state_j)]=ttest(targetOR_ratio{1,task}(:,state_i), targetOR_ratio{1,task}(:,state_j));                
        end
    end
end

targetOR_ratio=OR_MR.PDON;
h=[]; p=[];
for state =1:5
    for task_i=1:3
        for task_j=1:3
            [h{1,state}(task_i, task_j), p{1,state}(task_i,task_j)]=ttest(targetOR_ratio{1,task_i}(:,state), targetOR_ratio{1,task_j}(:,state));                
        end
    end
end


%% Plot Group level differences in brain measures: OR ratio: FIGURE 3G,H
group_names={'HC', 'PDOFF', 'PDON'};
condition_names={'LL', 'HL', 'DL'};
group_idx = [ones(1,length(HC_idx)), 2.*ones(1,length(PDOFF_idx)), 3.*ones(1,length(PDON_idx))];

% Encoding phase    
for task=1:3
    OR_E_ratios(:,task)=TemEvol_SubjStat_E.OR{1,task}(:,4)./TemEvol_SubjStat_E.OR{1,task}(:,1);
end

% M&R phase
for task=1:3
    OR_MR_ratios(:,task)=TemEvol_SubjStat_MR.OR{1,task}(:,2)./TemEvol_SubjStat_MR.OR{1,task}(:,5);
end

figure;
subplot(1,2,1);
plot_bargraph([OR_E_ratios], group_idx, condition_names, group_names, 'OR Ratio', [0, 4], 11, 'Encoding: OR ratio of State 1 and 2' );
subplot(1,2,2);
plot_bargraph([OR_MR_ratios], group_idx, condition_names, group_names, 'OR Ratio', [0, 3], 11, 'Encoding: OR ratio of State 3 and 4' );
sgtitle('OR ratio')

% stats
targetOR_ratio=OR_E_ratios;
h=[]; p=[];
for task=1:3
    for group_i=1:3
        for group_j=1:3
            if group_i==2 && group_j==3 
                [h{1,task}(group_i, group_j), p{1,task}(group_i, group_j)]=ttest(targetOR_ratio(find(group_idx==group_i), task), targetOR_ratio(find(group_idx==group_j), task));
            elseif group_i==3 && group_j==2
                [h{1,task}(group_i, group_j), p{1,task}(group_i, group_j)]=ttest(targetOR_ratio(find(group_idx==group_i), task), targetOR_ratio(find(group_idx==group_j), task));
            else
                [h{1,task}(group_i, group_j), p{1,task}(group_i, group_j)]=ttest2(targetOR_ratio(find(group_idx==group_i), task), targetOR_ratio(find(group_idx==group_j), task));                
            end
        end
    end
end


%% Plot Group level differences in brain measures for cognitively normal (non-MCI) participants: OR ratio: FIGURE S6A, B
group_names={'HC', 'PDOFF', 'PDON'};
condition_names={'LL', 'HL', 'DL'};
group_idx = [ones(1,length(HC_idx)), 2.*ones(1,length(PDOFF_nonMCI_idx)), 3.*ones(1,length(PDON_nonMCI_idx))];

CogNormal_idx=[HC_idx; PDOFF_nonMCI_idx; PDON_nonMCI_idx];

% Encoding phase    
for task=1:3
    OR_E_CogNormal_ratios(:,task)=TemEvol_SubjStat_E.OR{1,task}(CogNormal_idx, 4)./TemEvol_SubjStat_E.OR{1,task}(CogNormal_idx,1);
end

% M&R phase
for task=1:3
    OR_MR_CogNormal_ratios(:,task)=TemEvol_SubjStat_MR.OR{1,task}(CogNormal_idx,2)./TemEvol_SubjStat_MR.OR{1,task}(CogNormal_idx,5);
end

figure;
subplot(1,2,1);
plot_bargraph([OR_E_CogNormal_ratios], group_idx, condition_names, group_names, 'OR Ratio', [0, 3], 11, 'Encoding: OR ratio of State 1 and 2 in cogitively normal participants' );
subplot(1,2,2);
plot_bargraph([OR_MR_CogNormal_ratios], group_idx, condition_names, group_names, 'OR Ratio', [0, 3], 11, 'Encoding: OR ratio of State 3 and 4 in cognitively normal participants' );
sgtitle('OR ratio in Cognitively normal PD')


% stats
targetOR_ratio=OR_MR_CogNormal_ratios;
h=[]; p=[];
for task=1:3
    for group_i=1:3
        for group_j=1:3
            if group_i==2 && group_j==3 
                [h{1,task}(group_i, group_j), p{1,task}(group_i, group_j)]=ttest(targetOR_ratio(find(group_idx==group_i), task), targetOR_ratio(find(group_idx==group_j), task));
            elseif group_i==3 && group_j==2
                [h{1,task}(group_i, group_j), p{1,task}(group_i, group_j)]=ttest(targetOR_ratio(find(group_idx==group_i), task), targetOR_ratio(find(group_idx==group_j), task));
            else
                [h{1,task}(group_i, group_j), p{1,task}(group_i, group_j)]=ttest2(targetOR_ratio(find(group_idx==group_i), task), targetOR_ratio(find(group_idx==group_j), task));                
            end
        end
    end
end


%% CCA Analysis to investigate relation between difference in OR and difference in behavioral performance of PDOFF&PDON
target_task=3; %Distractor load condition
% Brain feature is composed of difference bewteen OR of PDOFF and PDON.
Brain_Features = [TemEvol_SubjStat_E.OR{1,target_task}(PDON_idx,:)-TemEvol_SubjStat_E.OR{1,target_task}(PDOFF_idx,:),...
    TemEvol_SubjStat_MR.OR{1,target_task}(PDON_idx,:)-TemEvol_SubjStat_MR.OR{1,target_task}(PDOFF_idx,:)];

% Cognitive Feature is composed of diff in ACC and diff in RT between PDOFF and PDON.
CogFlex_Features = [RT(PDON_idx,target_task)-RT(PDOFF_idx,target_task),...
    Accuracy(PDON_idx,target_task)-Accuracy(PDOFF_idx,target_task)];

% Create control variables (age, gender, education)
SubjGender = 2 * strcmp(SubjInfo.Gender, 'M') - 1; % Encode gender as 1 (M), -1 (F)
control_vars = [SubjInfo.Age(PDOFF_idx), SubjGender(PDOFF_idx), SubjInfo.Education_y(PDOFF_idx)];

% Regress out control variables
Brain_Features = func_residual(Brain_Features, control_vars);
CogFlex_Features = func_residual(CogFlex_Features, control_vars);

% Normalize data using vectorized operations
Brain_Features_Norm=zeros(size(Brain_Features));
CogFlex_Features_Norm=zeros(size(CogFlex_Features));
for i=1:size(Brain_Features,2)
    Brain_Features_Norm(:,i)=func_normalizeData(Brain_Features(:,i));
end
for i=1:size(CogFlex_Features,2)
    CogFlex_Features_Norm(:,i)=func_normalizeData(CogFlex_Features(:,i));
end

% Perform canonical correlation analysis (CCA)
[BrainWei, BehaviorWei, CCACorr, BrainCV, BehaviorCV, stats] = canoncorr(Brain_Features_Norm, CogFlex_Features_Norm);

% Permutation test
Nperm = 1000; % Number of permutations
Nobsv = size(Brain_Features, 1); % Number of observations
CCACor_perm = zeros(Nperm, size(CCACorr, 2));

for j = 1:Nperm
    perm_idx = randperm(Nobsv);
    [~, ~, CCACor_perm(j, :), ~, ~, ~] = canoncorr(Brain_Features_Norm, CogFlex_Features_Norm(perm_idx, :));
end

% Calculate corrected p-values
CCAsignificance = mean(CCACor_perm(:, 1) >= CCACorr, 1);
CCAmode = find(CCAsignificance < 0.05); % Significant CCA components


% Check significant contribution of each variable to CC
r_brain=[]; p_brain=[];
for i=1:size(Brain_Features,2)
    [r_brain(i,1), p_brain(i,1)]=corr(BrainCV(:,1), Brain_Features_Norm(:,i));
end

r_behavior=[]; p_behavior=[];
for i=1:size(CogFlex_Features,2)
    [r_behavior(i,1), p_behavior(i,1)]=corr(BehaviorCV(:,1), CogFlex_Features(:,i));
end

% Figure 4B and D: Significant mode of covariation between OR of latent brain states and cognitive-behavioral measures
figure;
subplot(2,2,[1,3]);
scatter(BrainCV(:,CCAmode(1)), BehaviorCV(:,CCAmode(1)), 50, 'filled');
hold on
l = lsline;
set(l,'LineWidth', 2);
ylim([-3, 3]); xlim([-3, 3]);
legend(sprintf('r: %f\nP: %f', CCACorr(CCAmode(1)), CCAsignificance(CCAmode(1))), 'FontSize', 11);
xlabel('Brain Canonical Variate(OR)');
ylabel('Behavioral Canonical Variate');
mycolormap = customcolormap(linspace(0,1,6), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd'});
colorbar('southoutside');
colormap(mycolormap);
% Bar plot of significant brain weights
BrainWei_Sig = (BrainWei([StateReorder, (StateReorder+5)],1).*double(p_brain([StateReorder, (StateReorder+5)])<0.05)).*-1; 
subplot(2,2,2); bar(BrainWei_Sig); title('Weigths Brain','FontSize',13);  % ylim([-3, 5]);
set(gca,'xtick',[1:10],'xticklabel',{'S1E', 'S2E', 'S3E', 'S4E', 'S5E','S1MR', 'S2MR', 'S3MR', 'S4MR', 'S5MR'});
% Bar plot of significant behavioral weights
BehaviorWei_Sig = (BehaviorWei(:,1).*double(p_behavior(:,1)<0.05)).*-1; 
subplot(2,2,4); bar(BehaviorWei_Sig); title('Weights Behavior','FontSize',13);  ylim([-3,5]);
set(gca,'xtick',[1:2],'xticklabel',{'RT', 'ACC'});
sgtitle('CCA analysis: Multivariate occupancy rate-Behavior relationship');

% Figure 4C Predicting brain-behavior relationship using CCA 
[cca_pred_r, null_cca_pred_r, pred_x, pred_y, brainWei_pred] = func_cca_pred_analysis_LOOCV(Brain_Features_Norm,CogFlex_Features_Norm);
p=sum(null_cca_pred_r(:,1) > cca_pred_r(1))/length(null_cca_pred_r(:,1));
if p < 0.05
    fprintf('predicted R: %f, Permutation pval: %f\n' ,cca_pred_r(1), p);
    figure; scatter(pred_x(:,1), pred_y(:,1), 50, CogFlex_Features(:,1), 'filled');    
	hold on; l = lsline; set(l,'LineWidth', 2);

    legend(sprintf('r: %f\nP: %f', cca_pred_r(1), p), 'FontSize', 11);
    xlabel('Brain Canonical Variate');
    ylabel('Behavior Canonical Variate');
    sgtitle('Cross validation analysis: Predicted brain state dynamics and task performance');
else
    fprintf('Prediction result is not significant\n');
end

%% Figure 7. Relationship with LEDD dosage: Testing dopamine overdose hypothesis with brain & behavior CV

% nonlinear relationship
% Testing dopamin overdose hypothesis with brain CV
[p, S]=polyfit(SubjInfo.LEDD(PDOFF_idx), BrainCV(:,1), 2);
Rsquare= 1 - (S.normr/norm(BrainCV(:,1) - mean(BrainCV(:,1))))^2;
x1=linspace(min(SubjInfo.LEDD(PDOFF_idx))-100, max(SubjInfo.LEDD(PDOFF_idx))+100);
y1=polyval(p,x1);

subplot(1,2,1); scatter(SubjInfo.LEDD(PDOFF_idx), BrainCV(:,1), 50, CogFlex_Features(:,1), 'filled');
hold on; plot(x1,y1);
xlim([min(SubjInfo.LEDD(PDOFF_idx))-100, max(SubjInfo.LEDD(PDOFF_idx))+100]);
ylim([-3 3]);
title('Relationship between LEDD and brain CV')
mdl=polyfitn(SubjInfo.LEDD(PDOFF_idx), BrainCV(:,1), 2);
disp(mdl);

% Testing dopamin overdose hypothesis with behavior CV
[p, S]=polyfit(SubjInfo.LEDD(PDOFF_idx), BehaviorCV(:,1), 2);
Rsquare= 1 - (S.normr/norm(BehaviorCV(:,1) - mean(BehaviorCV(:,1))))^2;
x1=linspace(min(SubjInfo.LEDD(PDOFF_idx))-100, max(SubjInfo.LEDD(PDOFF_idx))+100);
y1=polyval(p,x1);

subplot(1,2,2);
scatter(SubjInfo.LEDD(PDOFF_idx), BehaviorCV(:,1), 50, CogFlex_Features(:,1), 'filled');
hold on; plot(x1,y1);
xlim([min(SubjInfo.LEDD(PDOFF_idx))-100, max(SubjInfo.LEDD(PDOFF_idx))+100]);
ylim([-3 3]);
mdl=polyfitn(SubjInfo.LEDD(PDOFF_idx), BehaviorCV(:,1), 2);
title('Relationship between LEDD and behavior CV')
disp(mdl)
sgtitle('Medication dosage effect');


%% Figure S8 Relationship with LEDD dosage: Temporal property & LEDD

%Brain canonical variate: 
figure; 
%linear relationship
[p, S]=polyfit(SubjInfo.LEDD(PDOFF_idx), BrainCV(:,1), 1);
Rsquare= 1 - (S.normr/norm(BrainCV(:,1) - mean(BrainCV(:,1))))^2;
x1=linspace(min(SubjInfo.LEDD(PDOFF_idx))-100, max(SubjInfo.LEDD(PDOFF_idx))+100);
y1=polyval(p,x1);
subplot(2,2,1); scatter(SubjInfo.LEDD(PDOFF_idx), BrainCV(:,1), 50, 'filled');
hold on; plot(x1,y1);
xlim([min(SubjInfo.LEDD(PDOFF_idx))-100, max(SubjInfo.LEDD(PDOFF_idx))+100]);
ylim([-3 3]);
title('Linear relationship between LEDD and brain CV')
mdl=polyfitn(SubjInfo.LEDD(PDOFF_idx), BrainCV(:,1), 1);
disp(mdl);

% nonlinear relationship
[p, S]=polyfit(SubjInfo.LEDD(PDOFF_idx), BrainCV(:,1), 2);
Rsquare= 1 - (S.normr/norm(BrainCV(:,1) - mean(BrainCV(:,1))))^2;
x1=linspace(min(SubjInfo.LEDD(PDOFF_idx))-100, max(SubjInfo.LEDD(PDOFF_idx))+100);
y1=polyval(p,x1);
subplot(2,2,2); scatter(SubjInfo.LEDD(PDOFF_idx), BrainCV(:,1), 50, 'filled');
hold on; plot(x1,y1);
xlim([min(SubjInfo.LEDD(PDOFF_idx))-100, max(SubjInfo.LEDD(PDOFF_idx))+100]);
ylim([-3 3]);
title('Nonlinear relationship between LEDD and brain CV')
mdl=polyfitn(SubjInfo.LEDD(PDOFF_idx), BrainCV(:,1), 2);
disp(mdl);

%Behavioral canonical variate:
%linear relationship
[p, S]=polyfit(SubjInfo.LEDD(PDOFF_idx), BehaviorCV(:,1), 1);
Rsquare= 1 - (S.normr/norm(BehaviorCV(:,1) - mean(BehaviorCV(:,1))))^2;
x1=linspace(min(SubjInfo.LEDD(PDOFF_idx))-100, max(SubjInfo.LEDD(PDOFF_idx))+100);
y1=polyval(p,x1);
subplot(2,2,3); scatter(SubjInfo.LEDD(PDOFF_idx), BehaviorCV(:,1), 50, 'filled');
hold on; plot(x1,y1);
xlim([min(SubjInfo.LEDD(PDOFF_idx))-100, max(SubjInfo.LEDD(PDOFF_idx))+100]);
ylim([-3 3]);
title('Linear relationship between LEDD and behavioral CV')
mdl=polyfitn(SubjInfo.LEDD(PDOFF_idx), BehaviorCV(:,1), 1);
disp(mdl)

%nonlinear relationship
[p, S]=polyfit(SubjInfo.LEDD(PDOFF_idx), BehaviorCV(:,1), 2);
Rsquare= 1 - (S.normr/norm(BehaviorCV(:,1) - mean(BehaviorCV(:,1))))^2;
x1=linspace(min(SubjInfo.LEDD(PDOFF_idx))-100, max(SubjInfo.LEDD(PDOFF_idx))+100);
y1=polyval(p,x1);
subplot(2,2,4); scatter(SubjInfo.LEDD(PDOFF_idx), BehaviorCV(:,1), 50, 'filled');
hold on; plot(x1,y1);
xlim([min(SubjInfo.LEDD(PDOFF_idx))-100, max(SubjInfo.LEDD(PDOFF_idx))+100]);
ylim([-3 3]);
title('Nonlinear relationship between LEDD and behavioral CV')
mdl=polyfitn(SubjInfo.LEDD(PDOFF_idx), BehaviorCV(:,1), 2);
disp(mdl)

sgtitle('Relationship with LEDD dosage: Temporal property & LEDD');


%% Figure S6C&D: CCA Analysis to investigate relation between difference in OR and difference in behavioral performance in COGNITIVELY NORMAL PDOFF&PDON

target_task=3; %Distractor load condition
% Brain feature is composed of difference bewteen OR of PDOFF and PDON.
Brain_Features = [TemEvol_SubjStat_E.OR{1,target_task}(PDON_nonMCI_idx,:)-TemEvol_SubjStat_E.OR{1,target_task}(PDOFF_nonMCI_idx,:),...
    TemEvol_SubjStat_MR.OR{1,target_task}(PDON_nonMCI_idx,:)-TemEvol_SubjStat_MR.OR{1,target_task}(PDOFF_nonMCI_idx,:)];

% Cognitive Feature is composed of diff in ACC and diff in RT between PDOFF and PDON.
CogFlex_Features = [RT(PDON_nonMCI_idx,target_task)-RT(PDOFF_nonMCI_idx,target_task),...
    Accuracy(PDON_nonMCI_idx,target_task)-Accuracy(PDOFF_nonMCI_idx,target_task)];

% Create control variables (age, gender, education)
SubjGender = 2 * strcmp(SubjInfo.Gender, 'M') - 1; % Encode gender as 1 (M), -1 (F)
control_vars = [SubjInfo.Age(PDOFF_nonMCI_idx), SubjGender(PDOFF_nonMCI_idx), SubjInfo.Education_y(PDOFF_nonMCI_idx)];

% Regress out control variables
Brain_Features = func_residual(Brain_Features, control_vars);
CogFlex_Features = func_residual(CogFlex_Features, control_vars);

% Normalize data using vectorized operations
Brain_Features_Norm=zeros(size(Brain_Features));
CogFlex_Features_Norm=zeros(size(CogFlex_Features));
for i=1:size(Brain_Features,2)
    Brain_Features_Norm(:,i)=func_normalizeData(Brain_Features(:,i));
end
for i=1:size(CogFlex_Features,2)
    CogFlex_Features_Norm(:,i)=func_normalizeData(CogFlex_Features(:,i));
end

% Perform canonical correlation analysis (CCA)
[BrainWei, BehaviorWei, CCACorr, BrainCV, BehaviorCV, stats] = canoncorr(Brain_Features_Norm, CogFlex_Features_Norm);

% Permutation test
Nperm = 1000; % Number of permutations
Nobsv = size(Brain_Features, 1); % Number of observations
CCACor_perm = zeros(Nperm, size(CCACorr, 2));

for j = 1:Nperm
    perm_idx = randperm(Nobsv);
    [~, ~, CCACor_perm(j, :), ~, ~, ~] = canoncorr(Brain_Features_Norm, CogFlex_Features_Norm(perm_idx, :));
end

% Calculate corrected p-values
CCAsignificance = mean(CCACor_perm(:, 1) >= CCACorr, 1);
CCAmode = find(CCAsignificance < 0.05); % Significant CCA components

% Check significant contribution of each variable to CC
r_brain=[]; p_brain=[];
for i=1:size(Brain_Features,2)
    [r_brain(i,1), p_brain(i,1)]=corr(BrainCV(:,1), Brain_Features_Norm(:,i));
end

r_behavior=[]; p_behavior=[];
for i=1:size(CogFlex_Features,2)
    [r_behavior(i,1), p_behavior(i,1)]=corr(BehaviorCV(:,1), CogFlex_Features(:,i));
end

% Plotting: Significant mode of covariation between OR of latent brain states and cognitive-behavioral measures
figure;
subplot(2,2,[1,3]);
scatter(BrainCV(:,CCAmode(1)), BehaviorCV(:,CCAmode(1)), 50, 'filled');
hold on
l = lsline;
set(l,'LineWidth', 2);
ylim([-3, 3]); xlim([-3, 3]);
legend(sprintf('r: %f\nP: %f', CCACorr(CCAmode(1)), CCAsignificance(CCAmode(1))), 'FontSize', 11);
xlabel('Brain Canonical Variate(OR)');
ylabel('Behavioral Canonical Variate');
mycolormap = customcolormap(linspace(0,1,6), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd'});
colorbar('southoutside');
colormap(mycolormap);
% Bar plot of significant brain weights
BrainWei_Sig = BrainWei([StateReorder, (StateReorder+5)],1).*double(p_brain([StateReorder, (StateReorder+5)])<0.05); 
subplot(2,2,2); bar(BrainWei_Sig); title('Weigths Brain','FontSize',13);  % ylim([-3, 5]);
set(gca,'xtick',[1:10],'xticklabel',{'S1E', 'S2E', 'S3E', 'S4E', 'S5E','S1MR', 'S2MR', 'S3MR', 'S4MR', 'S5MR'});
% Bar plot of significant behavioral weights
BehaviorWei_Sig = BehaviorWei(:,1).*double(p_behavior(:,1)<0.05); 
subplot(2,2,4); bar(BehaviorWei_Sig); title('Weights Behavior','FontSize',13);  ylim([-3,5]);
set(gca,'xtick',[1:2],'xticklabel',{'RT', 'ACC'});
sgtitle('CCA analysis: Multivariate occupancy rate-Behavior relationship in PD-nonMCI');

%% Examine functional connectivity difference of each state
% convert covariance matrix estimated from BSDS model to Fisher's transformed Pearsons correlation matrix and detrmine functional connectivity each region of each state.

fc_istate=zeros(num.ROI, num.ROI, num.Subj, num.State);
fc_istate_z=zeros(num.ROI, num.ROI, num.Subj, num.State);
for subj=1:num.Subj
    for state= 1:num.State
        k = dominantStates(state);
        if ~isempty(model_subjwise.estimated_covariance{1,subj})
        est_cov = model_subjwise.estimated_covariance{1,subj}{1,k};
        Corr = corrcov(est_cov);
        Corr = Corr - diag(diag(Corr));
        Corr_z = log((1+Corr)./(1-Corr)); % Apply Fisher's z-transformation to the correlation coefficient
        fc_istate(:,:,subj,state) = Corr;
        fc_istate_z(:,:,subj,state) = Corr_z;
        end
    end
end

PDOFF_idx(23)=[];
PDON_idx(23)=[];

%% Figure 5 & S8
StateOfInterest=[1,4]; %State 1: Encoding related state, State 4: Maintenance&retrieval related state
for i=1:length(StateOfInterest)
    state=StateReorder(StateOfInterest(i));
    tscore_range=[-5 5];
    figure;
    subplot(1,3,1);
    pval=0.05;
    [sig_t_mtx_diff, t_mtx, p_mtx] = func_paired_ttest_corr_mtx_nonPair((squeeze(fc_istate_z(:,:,PDOFF_idx,state))), (squeeze(fc_istate_z(:,:,HC_idx,state))), pval);
    plot_connectivity_diff(sig_t_mtx_diff, num.ROI, ROI_names,tscore_range, 'Conn of PDOFF compared to HC');

    subplot(1,3,2);
    [sig_t_mtx_diff, t_mtx, p_mtx] = func_paired_ttest_corr_mtx_nonPair((squeeze(fc_istate_z(:,:,PDON_idx,state))),(squeeze(fc_istate_z(:,:,HC_idx,state))), pval);
    plot_connectivity_diff(sig_t_mtx_diff, num.ROI, ROI_names,tscore_range, 'Conn of PDON compared to HC');

    subplot(1,3,3);
    [sig_t_mtx_diff, t_mtx, p_mtx] = func_paired_ttest_corr_mtx(squeeze((fc_istate_z(:,:,PDOFF_idx,state))),(squeeze(fc_istate_z(:,:,PDON_idx,state))), pval);
    plot_connectivity_diff(sig_t_mtx_diff, num.ROI, ROI_names,tscore_range, 'Conn of PDOFF compared to PDON');
end


%% Figure 6B&D: CCA Analysis to investigate relation between difference in connectivity and in behavioral performance of PDOFF&PDON

target_task=3; %Distractor load condition
% Brain feature is composed of difference bewteen FC of PDOFF and PDON.
StateOfInterest=[1,4];
Brain_Features=[];
pval=0.05;
 for i=1:length(StateOfInterest)
     state=StateReorder(StateOfInterest(i));
     [sig_mtx_diff_PDOFFnON, ~, p_mtx] = func_paired_ttest_corr_mtx(squeeze((fc_istate_z(:,:,PDON_idx,state))),(squeeze(fc_istate_z(:,:,PDOFF_idx,state))), pval);

     [a,b]=find(triu(sig_mtx_diff_PDOFFnON)~=0);
    for i=1:length(a)
        Brain_Features = [Brain_Features, (squeeze(fc_istate_z(a(i),b(i),PDON_idx,state)))-(squeeze(fc_istate_z(a(i),b(i),PDOFF_idx,state)))];
    end
end

% Cognitive Feature is composed of diff in ACC and diff in RT between PDOFF and PDON.
CogFlex_Features = [RT(PDON_idx,target_task)-RT(PDOFF_idx,target_task),...
    Accuracy(PDON_idx,target_task)-Accuracy(PDOFF_idx,target_task)];

% Create control variables (age, gender, education)
SubjGender = 2 * strcmp(SubjInfo.Gender, 'M') - 1; % Encode gender as 1 (M), -1 (F)
control_vars = [SubjInfo.Age(PDOFF_idx), SubjGender(PDOFF_idx), SubjInfo.Education_y(PDOFF_idx)];

% Regress out control variables
Brain_Features = func_residual(Brain_Features, control_vars);
CogFlex_Features = func_residual(CogFlex_Features, control_vars);

% Normalize data using vectorized operations
Brain_Features_Norm=zeros(size(Brain_Features));
CogFlex_Features_Norm=zeros(size(CogFlex_Features));
for i=1:size(Brain_Features,2)
    Brain_Features_Norm(:,i)=func_normalizeData(Brain_Features(:,i));
end
for i=1:size(CogFlex_Features,2)
    CogFlex_Features_Norm(:,i)=func_normalizeData(CogFlex_Features(:,i));
end

% Perform canonical correlation analysis (CCA)
[BrainWei, BehaviorWei, CCACorr, BrainCV, BehaviorCV, stats] = canoncorr(Brain_Features_Norm, CogFlex_Features_Norm);

% Permutation test
Nperm = 1000; % Number of permutations
Nobsv = size(Brain_Features, 1); % Number of observations
CCACor_perm = zeros(Nperm, size(CCACorr, 2));

for j = 1:Nperm
    perm_idx = randperm(Nobsv);
    [~, ~, CCACor_perm(j, :), ~, ~, ~] = canoncorr(Brain_Features_Norm, CogFlex_Features_Norm(perm_idx, :));
end

% Calculate corrected p-values
CCAsignificance = mean(CCACor_perm(:, 1) >= CCACorr, 1);
CCAmode = find(CCAsignificance < 0.05); % Significant CCA components

% Check significant contribution of each variable to CC
r_brain=[]; p_brain=[];
for i=1:size(Brain_Features,2)
    [r_brain(i,1), p_brain(i,1)]=corr(BrainCV(:,1), Brain_Features(:,i));
end
r_behavior=[]; p_behavior=[];
for i=1:size(CogFlex_Features,2)
    [r_behavior(i,1), p_behavior(i,1)]=corr(BehaviorCV(:,1), CogFlex_Features(:,i));
end

% Figure 6B and D: Significant mode of covariation between functional connectivity of latent brain states and cognitive-behavioral measures
figure;
subplot(2,2,[1,3]);
scatter(BrainCV(:,CCAmode(1)), BehaviorCV(:,CCAmode(1)), 50, 'filled');
hold on
l = lsline;
set(l,'LineWidth', 2);
xlim([-3 3]); ylim([-3 3]);
legend(sprintf('r: %f\nP: %f', CCACorr(CCAmode(1)), CCAsignificance(CCAmode(1))), 'FontSize', 11);
xlabel('Brain Canonical Variate(');
ylabel('Behavioral Canonical Variate');
mycolormap = customcolormap(linspace(0,1,6), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd'});
colorbar('southoutside');
colormap(mycolormap);
BrainWei_Sig = BrainWei(:,1).*double(p_brain(:,1)<0.05); 
subplot(2,2,2); bar(BrainWei_Sig); title('Weigths Brain','FontSize',13);  % ylim([-3, 5]);
BehaviorWei_Sig = BehaviorWei(:,1).*double(p_behavior(:,1)<0.05); 
subplot(2,2,4); bar(BehaviorWei_Sig); title('Weights Behavior','FontSize',13);  ylim([-3,5]);
sgtitle('CCA analysis: Multivariate functional connectivity-Behavior relationship');

% Figure 6C: Predicting brain-behavior relationship using CCA
[cca_pred_r, null_cca_pred_r, pred_x, pred_y, brainWei_pred] = func_cca_pred_analysis_LOOCV(Brain_Features_Norm,CogFlex_Features_Norm);
p=sum(null_cca_pred_r(:,1) > cca_pred_r(1))/length(null_cca_pred_r(:,1));
if p < 0.05
    fprintf('predicted R: %f, Permutation pval: %f\n' ,cca_pred_r(1), p);
    figure; scatter(pred_x(:,1), pred_y(:,1), 50, CogFlex_Features(:,1), 'filled');    
	hold on; l = lsline; set(l,'LineWidth', 2);

    legend(sprintf('r: %f\nP: %f', cca_pred_r(1), p), 'FontSize', 11);
    xlabel('Brain Canonical Variate');
    ylabel('Behavior Canonical Variate');
    sgtitle('Cross validation analysis: Predicted brain state dynamics and task performance');
else
    fprintf('Prediction result is not significant\n');
end


%% Figure S9 Relationship with LEDD dosage: Temporal property & LEDD
%Brain canonical variate: 
figure; 
%linear relationship
[p, S]=polyfit(SubjInfo.LEDD(PDOFF_idx), BrainCV(:,1), 1);
Rsquare= 1 - (S.normr/norm(BrainCV(:,1) - mean(BrainCV(:,1))))^2;
x1=linspace(min(SubjInfo.LEDD(PDOFF_idx))-100, max(SubjInfo.LEDD(PDOFF_idx))+100);
y1=polyval(p,x1);
subplot(2,2,1); scatter(SubjInfo.LEDD(PDOFF_idx), BrainCV(:,1), 50, 'filled');
hold on; plot(x1,y1);
xlim([min(SubjInfo.LEDD(PDOFF_idx))-100, max(SubjInfo.LEDD(PDOFF_idx))+100]);
ylim([-3 3]);
title('Linear relationship between LEDD and brain CV')
mdl=polyfitn(SubjInfo.LEDD(PDOFF_idx), BrainCV(:,1), 1);
disp(mdl);

% nonlinear relationship
[p, S]=polyfit(SubjInfo.LEDD(PDOFF_idx), BrainCV(:,1), 2);
Rsquare= 1 - (S.normr/norm(BrainCV(:,1) - mean(BrainCV(:,1))))^2;
x1=linspace(min(SubjInfo.LEDD(PDOFF_idx))-100, max(SubjInfo.LEDD(PDOFF_idx))+100);
y1=polyval(p,x1);
subplot(2,2,2); scatter(SubjInfo.LEDD(PDOFF_idx), BrainCV(:,1), 50, 'filled');
hold on; plot(x1,y1);
xlim([min(SubjInfo.LEDD(PDOFF_idx))-100, max(SubjInfo.LEDD(PDOFF_idx))+100]);
ylim([-3 3]);
title('Nonlinear relationship between LEDD and brain CV')
mdl=polyfitn(SubjInfo.LEDD(PDOFF_idx), BrainCV(:,1), 2);
disp(mdl);

%Behavioral canonical variate:
%linear relationship
[p, S]=polyfit(SubjInfo.LEDD(PDOFF_idx), BehaviorCV(:,1), 1);
Rsquare= 1 - (S.normr/norm(BehaviorCV(:,1) - mean(BehaviorCV(:,1))))^2;
x1=linspace(min(SubjInfo.LEDD(PDOFF_idx))-100, max(SubjInfo.LEDD(PDOFF_idx))+100);
y1=polyval(p,x1);
subplot(2,2,3); scatter(SubjInfo.LEDD(PDOFF_idx), BehaviorCV(:,1), 50, 'filled');
hold on; plot(x1,y1);
xlim([min(SubjInfo.LEDD(PDOFF_idx))-100, max(SubjInfo.LEDD(PDOFF_idx))+100]);
ylim([-3 3]);
title('Linear relationship between LEDD and behavioral CV')
mdl=polyfitn(SubjInfo.LEDD(PDOFF_idx), BehaviorCV(:,1), 1);
disp(mdl)

%nonlinear relationship
[p, S]=polyfit(SubjInfo.LEDD(PDOFF_idx), BehaviorCV(:,1), 2);
Rsquare= 1 - (S.normr/norm(BehaviorCV(:,1) - mean(BehaviorCV(:,1))))^2;
x1=linspace(min(SubjInfo.LEDD(PDOFF_idx))-100, max(SubjInfo.LEDD(PDOFF_idx))+100);
y1=polyval(p,x1);
subplot(2,2,4); scatter(SubjInfo.LEDD(PDOFF_idx), BehaviorCV(:,1), 50, 'filled');
hold on; plot(x1,y1);
xlim([min(SubjInfo.LEDD(PDOFF_idx))-100, max(SubjInfo.LEDD(PDOFF_idx))+100]);
ylim([-3 3]);
title('Nonlinear relationship between LEDD and behavioral CV')
mdl=polyfitn(SubjInfo.LEDD(PDOFF_idx), BehaviorCV(:,1), 2);
disp(mdl)

sgtitle('Relationship with LEDD dosage: Spatial property & LEDD');