function Stim_response = Unit_Plot(input_cond,Stim_Type,fn, units, resp_win, artifact_win)

%% function which tests whether all units for a given experiment type are responsive to either whisker stimulus or Light stimulus
% INPUT : 
% input_cond  : condition of interest from ephys structure. eg
%              ephys_data.condition(1)
% Stim_Type   : Type of stimulation , 'Whsk','Opto','Both';
% fn          : file name and directory which figures are saved to
% Units       : Which units to test, by default all
% Resp window : window of time to test responsiveness to stimulus (window
%               inwhich spike rates and probablities are tested)
% Artifact_win: window of time to NaN out as artifact. 
%
% OUTPUT :
% StimResponse : structure with responsiveness statistics for stimulus. 
%
%
% Alexander von klemperer 2020
%%
% Default to all units
if nargin < 4 || isempty(units)
    units        = ':';
end

% Default to resp win from 6ms (after any artifacts) to 30ms (should capture
% most of the direct stimulus-driven activity
if nargin < 5 || isempty(resp_win)
    resp_win        = [0.005 0.015];
    control_win(1) = 0-resp_win(2);
    control_win(2) = 0-resp_win(1);
end

% Set any spikes during this window to NaN; -0.001 to 0.006 is where any piezo artifacts
% may occur
if nargin < 7
    artifact_win    = [0.00 0.000];
end

% Hardcoded for now:
rate_kernel_size    = 0.01;
opto_resp_win       = [0.000 0.05];
    opto_control_win(1) = 0-opto_resp_win(2);
    opto_control_win(2) = 0-opto_resp_win(1);


%% Code execution starts here

    opto_onsets             = [input_cond.LED_onset];
    n_delta_ts              = length(opto_onsets);
   
    % Fetch data for this condition
    this_cond                       = input_cond;
    
    % Find stimulus data for this condition
    this_t_whisk                    = this_cond.whisk_onset;
    this_t_opto                     = this_cond.LED_onset;
    if isfield('whisk_stimulator','this_cond')
        this_whisker_nr                 = this_cond.whisk_stimulator;
    else
        this_whisker_nr                 = this_cond.whisk_stim_nr;
    end
    
    if (Stim_Type == 'Opto') 
        this_t_whisk = 0;
    end;
        
    n_trials            	= this_cond.n_trials;
    
    delta_t               = this_t_opto - this_t_whisk;
    disp(delta_t);
    whisker_nr             = this_whisker_nr;
    opto_power            = this_cond.LED_power;
    
    % Get spike data and remove artifact spikes
    spikes                          = this_cond.spikes(units, :, :) - this_t_whisk;
    
    q_artifact                      = spikes > artifact_win(1) & spikes < artifact_win(2);
    spikes(q_artifact)              = NaN;
    
    if (Stim_Type == 'Whsk') | (Stim_Type == ['Both']) % if whisker stim or both 
    Whisk_Resp = Stim_Responsive(spikes,resp_win,control_win,n_trials,delta_t,true,'Whisker Stim Response',fn); % tests whether spike rate or probablity is signifactly larger following stimulus
    else
        Whisk_Resp = [];
    end;
    %% opto response
     if (Stim_Type == 'Opto') | (Stim_Type == 'Both')
    opto_spikes                  	= (spikes + this_t_whisk) - this_t_opto;
    opto_spike_rates(:,:) 	= spike_rates_individual(opto_spikes, opto_resp_win);
    Opto_Resp = Stim_Responsive(opto_spikes,opto_resp_win,opto_control_win,n_trials,delta_t,true,'Opto Stim Response',fn);
     else
         Opto_Resp = [];
     end;
     
     Stim_response.Opto_Resp = Opto_Resp;
     Stim_response.Whisk_Resp = Whisk_Resp;
end
