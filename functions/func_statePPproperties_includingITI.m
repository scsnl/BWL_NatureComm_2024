function [StatePP_GroupStat, StatePP_SubjStat, StatePP_everyTrial]=func_statePPproperties_includingITI(onsetNoffset, grpStatePP, num)


for task = 1:num.Task
    count = 1; idx = 1;
    for subj = 1:num.Subj
        for run = 1:num.Run
            isubj_onsetNoffset = onsetNoffset{1,task}((onsetNoffset{1,task}.idx_RunxSubj==idx),:);
            for m = 1:size(isubj_onsetNoffset,1)
            	monset_E = isubj_onsetNoffset.Encod_On(m)-1; %Adding ITI: round((1/0.49));
                moffset_E = isubj_onsetNoffset.Encod_Off(m)+1;
                
                monset_M = isubj_onsetNoffset.Maint_On(m);
%                 moffset_M = isubj_onsetNoffset.Maint_Off(m);
                moffset_M = monset_M + 13; %setting to constant value of 6sec wich is an average maintenance phase
                
                monset_R = isubj_onsetNoffset.Respons_On(m);
                moffset_R = isubj_onsetNoffset.Respons_Off(m)+ round((8/0.49)); %Adding ITI:
         
                monset_TB = isubj_onsetNoffset.Encod_On(m); %Task block(TB): From Encod_On to Respons_Off
                moffset_TB = isubj_onsetNoffset.Respons_Off(m);

                %%    
                StatePP.allTrials.Encod{1,task}(count,:,:) = squeeze(grpStatePP(idx,monset_E:moffset_E,:)); %State PP of Encoding phase: [length of the phase x #of States]
                StatePP.allTrials.Maint{1,task}(count,:,:) = squeeze(grpStatePP(idx,monset_M:moffset_M,:)); %State PP of Maintanence phase: [length of the phase x #of States]
                StatePP.allTrials.Respons{1,task}(count,:,:) = squeeze(grpStatePP(idx,monset_R:moffset_R,:)); %State PP of Response phase: [length of the phase x #of States]
                count = count+1;
            end
            idx = idx+1;
        end
        
        %% Gather the StatePP of j_th subjects' s_th trial of k_th task
        idx_start = find(onsetNoffset{1,task}.idx_RunxSubj==(num.Run*subj-(num.Run-1)), 1, 'first');
        idx_end = find(onsetNoffset{1,task}.idx_RunxSubj==(num.Run*subj), 1, 'last');        
                
        StatePP.subj.EachTrial.Encod{1,task}{subj,1}=StatePP.allTrials.Encod{1,task}(idx_start:idx_end,:,:);
        StatePP.subj.EachTrial.Maint{1,task}{subj,1}=StatePP.allTrials.Maint{1,task}(idx_start:idx_end,:,:);
        StatePP.subj.EachTrial.Respons{1,task}{subj,1}=StatePP.allTrials.Respons{1,task}(idx_start:idx_end,:,:);
               
        %% Gather the average StatePPs of j_th subjects averaged over all trials of k_th task
        StatePP.subj.Avg.Encod{1,task}(subj,:,:)=squeeze(nanmean((StatePP.subj.EachTrial.Encod{1,task}{subj,1})));
        StatePP.subj.Avg.Maint{1,task}(subj,:,:)=squeeze(nanmean((StatePP.subj.EachTrial.Maint{1,task}{subj,1})));
        StatePP.subj.Avg.Respons{1,task}(subj,:,:)=squeeze(nanmean((StatePP.subj.EachTrial.Respons{1,task}{subj,1})));        
        StatePP.subj.Avg.wholeTB{1,task}(subj,:,:)=[squeeze(StatePP.subj.Avg.Encod{1,task}(subj,:,:));...
                                                  squeeze(StatePP.subj.Avg.Maint{1,task}(subj,:,:));... 
                                                  squeeze(StatePP.subj.Avg.Respons{1,task}(subj,:,:))];

    end
	%% Gather the average StatePPs of the group averaged over all trials of k_th task
    StatePP.group.Avg.Encod{1,task}=squeeze(nanmean(StatePP.allTrials.Encod{1,task}));
    StatePP.group.Avg.Maint{1,task}=squeeze(nanmean(StatePP.allTrials.Maint{1,task}));
    StatePP.group.Avg.Respons{1,task}=squeeze(nanmean(StatePP.allTrials.Respons{1,task}));
end

%%Combining Encod+Maint+Respons
for task=1:num.Task
    StatePP.group.Avg.wholeTB{1,task}=[StatePP.group.Avg.Encod{1,task};StatePP.group.Avg.Maint{1,task};StatePP.group.Avg.Respons{1,task}];
end

StatePP_GroupStat=StatePP.group;
StatePP_SubjStat=StatePP.subj;
StatePP_everyTrial=StatePP.allTrials;
