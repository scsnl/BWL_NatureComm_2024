function [OverallStat, SubjStat]=func_tempEvolProperties(onsetNoffset, grpTempEvol, dominantStates, targetPhase, num)

% relabel states for ease of computation
grpTempEvol_relabel = zeros(size(grpTempEvol));
for relabel = 1:num.State
    grpTempEvol_relabel(grpTempEvol == dominantStates(relabel)) = relabel;
end

OR = cell(1,num.Task);
MLT = cell(1,num.Task);
OR_whole = cell(1,num.Task);


for task = 1:num.Task
    count = 1;
    idx = 1;
    for subj = 1:num.Subj
        for run = 1:num.Run        
            isubj_onsetNoffset = onsetNoffset{1,task}((onsetNoffset{1,task}.idx_RunxSubj==idx),:);
            for m = 1:size(isubj_onsetNoffset,1)

                %% whole-phase onset & offset
                monset_whole = isubj_onsetNoffset.Encod_On(m);
                moffset_whole = isubj_onsetNoffset.Respons_Off(m);
                % Compute occupancy rate & mean life time of latent brain states
                state_transition_relabel_whole=grpTempEvol_relabel(idx, monset_whole:moffset_whole); 
                [OR_whole{1,task}(count,:), ~]=func_summary_stats_fast(state_transition_relabel_whole, 1:num.State);

                monset_E = isubj_onsetNoffset.Encod_On(m);
                moffset_E = isubj_onsetNoffset.Encod_Off(m);
                % Compute occupancy rate & mean life time of latent brain states
                state_transition_relabel_E=grpTempEvol_relabel(idx, monset_E:moffset_E); 
                [OR_E{1,task}(count,:), ~]=func_summary_stats_fast(state_transition_relabel_E, 1:num.State);

                monset_MR = isubj_onsetNoffset.Maint_On(m);
                moffset_MR = isubj_onsetNoffset.Respons_Off(m);
                % Compute occupancy rate & mean life time of latent brain states
                state_transition_relabel_MR=grpTempEvol_relabel(idx, monset_MR:moffset_MR); 
                [OR_MR{1,task}(count,:), ~]=func_summary_stats_fast(state_transition_relabel_MR, 1:num.State);


                %% phase specific onset & offset
                if strcmp(targetPhase, 'whole')==1
                    monset = isubj_onsetNoffset.Encod_On(m);
                    moffset = isubj_onsetNoffset.Respons_Off(m);                    
                elseif strcmp(targetPhase, 'encoding')==1
                    monset = isubj_onsetNoffset.Encod_On(m);
                    moffset = isubj_onsetNoffset.Encod_Off(m);  
                elseif strcmp(targetPhase, 'maintenance')==1
                    monset = isubj_onsetNoffset.Maint_On(m);
                    moffset = isubj_onsetNoffset.Maint_Off(m);
                elseif strcmp(targetPhase, 'retrieval')==1
                    monset = isubj_onsetNoffset.Respons_On(m);
%                     moffset = isubj_onsetNoffset.Respons_Off(m); % 1sec after Response onset
                    moffset = isubj_onsetNoffset.Respons_Off(m) + 2; % 1sec after Response onset+ add 1 more sec
                elseif strcmp(targetPhase, 'E&M')==1
                    monset = isubj_onsetNoffset.Encod_On(m);
                    moffset = isubj_onsetNoffset.Maint_Off(m);
                elseif strcmp(targetPhase, 'M&R')==1
                    monset = isubj_onsetNoffset.Maint_On(m);
%                     moffset = isubj_onsetNoffset.Respons_Off(m);
                    moffset = isubj_onsetNoffset.Respons_Off(m) + 2; % 1sec after Response onset+ add 1 more sec
                end

                state_transition_relabel=grpTempEvol_relabel(idx, monset:moffset);                 
                % Compute occupancy rate & mean life time of latent brain states
                [OR{1,task}(count,:), MLT{1,task}(count,:)]=func_summary_stats_fast(state_transition_relabel, 1:num.State);
                % compute state switching probability
                % y-axis:state(t), x-axis:state(t+1)
                transitionMatrix = full(sparse(state_transition_relabel(1:end-1),state_transition_relabel(2:end),1,num.State,num.State));
                sum_transitionMatrix = sum(transitionMatrix,2);
                
                % Normalize the transition matrix
                for state=1:num.State
                    transitionMatrix(state, :) = transitionMatrix(state,:)/sum_transitionMatrix(state);                    
                end
                StatTransProb{1,task}(count,:,:) = transitionMatrix;      
                
                count = count+1;
            end
            idx = idx+1;
        end


        % Average the OR, MLT, StatTransProb of ith_subjects' j_th task
        idx_start = find(onsetNoffset{1,task}.idx_RunxSubj==(num.Run*subj-(num.Run-1)), 1, 'first');
        idx_end = find(onsetNoffset{1,task}.idx_RunxSubj==(num.Run*subj), 1, 'last');
        idx_range = idx_start:idx_end;

        % if the OR of whole-phase is dominated by a single state, exclude that trial.

        temp_W=OR_whole{1,task}(idx_range,:);
        temp_E=OR_E{1,task}(idx_range,:);
        temp_MR=OR_MR{1,task}(idx_range,:);
        [a,b]=find(temp_W(:,3)>1);
        idx_range(a)=[];

        mean_OR_subj{1,task}(subj,:)=squeeze(mean(OR{1,task}(idx_range,:)));
        mean_MLT_subj{1,task}(subj,:)=squeeze(mean(MLT{1,task}(idx_range,:)));

        if length(idx_range)==1
            mean_StatSwitchProb_subj{1,task}(subj,:,:)=squeeze((StatTransProb{1,task}(idx_range,:,:)));
        else
            mean_StatSwitchProb_subj{1,task}(subj,:,:)=squeeze(nanmean(StatTransProb{1,task}(idx_range,:,:)));
        end
    end
end

for task = 1:num.Task
    mean_OR{1,task} = squeeze(mean(mean_OR_subj{1,task}));
    mean_MLT{1,task} = squeeze(mean(mean_MLT_subj{1,task}));
    mean_StatSwitchProb{1,task} = squeeze(nanmean(mean_StatSwitchProb_subj{1,task}, 1));
end

SubjStat.OR=mean_OR_subj;
SubjStat.MLT=mean_MLT_subj;
SubjStat.OREachTrial=OR;
SubjStat.MLTEachTrial=MLT;

SubjStat.StatTransProb=mean_StatSwitchProb_subj;

OverallStat.OR=mean_OR;
OverallStat.MLT=mean_MLT;
OverallStat.StatTransProb=mean_StatSwitchProb;



    