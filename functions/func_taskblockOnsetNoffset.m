function [onsetNoffset]=func_taskblockOnsetNoffset(isubj_task_dsgn, taskLoad, idx, PID, run, TR)

%% block length of the task: Onset of the encoding phase (stimulus onset) to 1 seconds after the onset of the probe (response onset) (Poston et al., 2015). 

%% stimulation onset & offset
encoding_Onset_vec = ceil(isubj_task_dsgn(4+taskLoad).Onset./TR); 
encoding_Offset_vec = encoding_Onset_vec + ceil(isubj_task_dsgn(4+taskLoad).Duration./TR)-1; % Duration of stimulus encoding is 2s (i.e., 5 frame considering TR)

%% maintenance onset & offset
maintenance_Onset_vec = encoding_Offset_vec + 1;
maintenance_Offset_vec = maintenance_Onset_vec + ceil(isubj_task_dsgn(8+taskLoad).Duration./TR)-1; % Duration of maintainence

%% response onset & offset 
response_Onset_vec = maintenance_Offset_vec + 1; 
response_Offset_vec = response_Onset_vec + ceil(1/TR)-1; % 1sec after the response onset

%% Merge
isubj_onsetNoffset = zeros(length(encoding_Onset_vec), 9);
itrial_numb=0;
for ions = 1:length(encoding_Onset_vec)
	itrial_numb = itrial_numb + 1;
    
    isubj_onsetNoffset(itrial_numb,1) = encoding_Onset_vec(ions);
    isubj_onsetNoffset(itrial_numb,2) = encoding_Offset_vec(ions);
    isubj_onsetNoffset(itrial_numb,3) = maintenance_Onset_vec(ions);
    isubj_onsetNoffset(itrial_numb,4) = maintenance_Offset_vec(ions);    
	isubj_onsetNoffset(itrial_numb,5) = response_Onset_vec(ions);
    isubj_onsetNoffset(itrial_numb,6) = response_Offset_vec(ions);
    isubj_onsetNoffset(itrial_numb,7) = idx;
    isubj_onsetNoffset(itrial_numb,8) = PID;
    isubj_onsetNoffset(itrial_numb,9) = run;
end
onsetNoffset = isubj_onsetNoffset;
