function[CSD,w,sinks,sources,LFP_sinks] = Current_Source_Density(Traces,CSD_window);
%% CSD,w = Current_Source_Density(Traces,CSD_window)
%  extracts current source density for LFP
% 
% REQUIRED INPUTS:
% 
% Traces  : N x M x C matrix, N= number of channels, M = Number of trials,
% C = time points 
% containing LFT recordings(mV)
% 
% OPTIONAL:
% 
% OUTPUT: 
% CSD : N x C array, N number of channels, C timepoints
% containing CSD (uA/mm^3) averaged by trials.
% 
% w:  NxC array of input LFPS averaged by trials. 
% 
%sinks : Nx1 array with the minimum CSD for each channel in the CSD_window
%time period
%
%sources:Nx1 array with the maximum CSD for each channel in the CSD_window
%time period

%phi = LFP at z
% t = time; 
%CSD at time t = sigma x ( LFP at 1 point shallower - 2xLFP + sigma at 1
%point deeper at time t)/ channel depth; 
% z = 25um
%CSD = uA/mm3
%
%Alex von klemperer - 2020
 Artifact_delay = 3;%delay to ignore as likely due to artifact
delta_z = 0.025; %um
sigma = -0.3; %mS/mm


for k = 1 : size(Traces,2)
    temp_traces = squeeze(Traces(:,k,:));
        t = 1:size(temp_traces,2); 
    z = 2:size(temp_traces,1)-1;
    CSD_temp(:,:) = sigma*(temp_traces(z-1,t) - 2*temp_traces(z,t) + temp_traces(z+1,t))/(delta_z^2);
    CSD(:,:,k) = CSD_temp;
end;
    
CSD = nanmean(CSD,3);
CSD_window = CSD_window+Artifact_delay;
w = squeeze(nanmean(Traces(:,:,:),2));

  for k = 1 : size(CSD,1) 
    sinks(k) = min(CSD(k,CSD_window(1):CSD_window(2)));
    sources(k) = max(CSD(k,CSD_window(1):CSD_window(2)));
    LFP_sinks(k) = min(w(k,CSD_window(1):CSD_window(2)));
   end;
end    