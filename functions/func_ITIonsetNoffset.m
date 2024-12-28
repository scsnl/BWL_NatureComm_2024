function [isubj_onsetNoffset]=func_ITIonsetNoffset(group_taskSeq, idx, PID, runType)
          

isubj_taskSeq=group_taskSeq(idx,:);
indx=find(isubj_taskSeq==0); %tasklist=0: ITI
finding_onset=[(find((indx(2:end)-indx(1:end-1))>1)+1)];
ITI_Onset = indx(finding_onset);
ITI_Onset(end) = [];
ITI_length = finding_onset(2:end)-finding_onset(1:end-1);

ITI_Offset = ITI_Onset + ITI_length-1; % Determined fixed (shared) duration
isubj_onsetNoffset=[ITI_Onset', ITI_Offset', idx*ones(length(ITI_Onset),1), PID*ones(length(ITI_Onset),1), runType*ones(length(ITI_Onset),1)];    
isubj_onsetNoffset(find(isubj_onsetNoffset(:,2)==1),:)=[]; % find ITI with duration 1
