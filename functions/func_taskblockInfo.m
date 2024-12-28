function [ITIonsetNoffset, TaskonsetNoffset] = func_taskblockInfo(task_design_info, num, TR)


% Get onset & duration of ITI
group_taskSeq=zeros(num.Subj*num.Run, num.Vol);
ITIonsetNoffset = cell(1,1);
onsetNoffset_merge=[];
idx = 1;
for subj =1:num.Subj
    PID = task_design_info(subj).SubjInfo.PID;
    for run = 1:num.Run
        isubj_task_dsgn = task_design_info(subj).taskdesign(run).taskgn;
       for taskLoad=1:num.Task %ONLY the correct trials. tasklist=1:LL, tasklist=2:HL, tasklist=3:DL 
            Onset=ceil(isubj_task_dsgn(4+taskLoad).Onset./TR); %Stimulus(Encoding) onsent
            Dur=ceil(isubj_task_dsgn(4+taskLoad).Blocklength./TR); %Dur: Block-length
            for j=1:length(Onset)
                group_taskSeq(idx, Onset(j):Onset(j)+Dur(j)-1)=taskLoad;
            end
        end
        isubj_onsetNoffset = func_ITIonsetNoffset(group_taskSeq, idx, PID, run);
        onsetNoffset_merge = [onsetNoffset_merge; isubj_onsetNoffset];

        idx = idx+1;
    end
end
ITIonsetNoffset = array2table(onsetNoffset_merge,...
'VariableNames',{'ITI_On','ITI_Off','idx_RunxSubj','PID', 'run'});    


% Get onset & duration of LL, HL, DL
TaskonsetNoffset = cell(1,num.Task);
for taskLoad=1:num.Task
    idx = 1;
    onsetNoffset_merge=[];
    for subj =1:num.Subj
        PID = task_design_info(subj).SubjInfo.PID;
        for run = 1:num.Run
            isubj_task_dsgn = task_design_info(subj).taskdesign(run).taskgn;
            isubj_onsetNoffset = func_taskblockOnsetNoffset(isubj_task_dsgn, taskLoad, idx, PID, run, TR);
            onsetNoffset_merge = [onsetNoffset_merge; isubj_onsetNoffset];
            idx = idx+1;
        end
    end
    TaskonsetNoffset{1,taskLoad} = array2table(onsetNoffset_merge,...
	'VariableNames',{'Encod_On','Encod_Off','Maint_On','Maint_Off','Respons_On','Respons_Off', 'idx_RunxSubj', 'PID', 'run'});    
end
