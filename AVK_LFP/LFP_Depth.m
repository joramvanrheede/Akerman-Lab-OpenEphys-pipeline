function LFP_out = LFP_Depth(experiment_type,MUE_data,OutputFn,plot_figs);

window_edges = 0.05; % time in seconds before and after stimulus onset/offset to keep;
CSD_window = 0.03; % window in which CSD to be measured following stimulus

channels_All = [1:32];
channels_L4 = [6:16];
channels_L5 = [15:32];

LFP_out = struct;
LFP_out.exp_type =experiment_type;
LFP_out.Chans_All = channels_All;    

        si = 0.001; % sampling interval in s
         window_edges = window_edges/si; 
        CSD_window = CSD_window/si;
        CSD_window = [window_edges+1 window_edges+1+CSD_window];

  switch experiment_type
    case 'Drive' % if drive experiment
         disp('Experiment Drive');
         cond_data = MUE_data.conditions(2); % first gets control drive condition.
        case 'Timing'
         disp('Experiment TIMING');
         cond_data = MUE_data.conditions(numel(MUE_data.conditions)); % takes last timing conditon (control condition)
    end;
        LFP_out.Cond_data = cond_data;
     
       % Performs CSD analysis on Layer 4 for whisk stimulus   
        [CSD_Depth,LFP_Depth,Sinks,Sources,LFP_Sinks] = CSD_analysis(window_edges,CSD_window,channels_All,cond_data,'Whisk',[OutputFn '\Figs\'],plot_figs);
        
       %% outputs
       LFP_out.CSD_Depth = CSD_Depth;
       LFP_out.LFP_Depth = LFP_Depth;
       LFP_out.Sinks = Sinks;
       LFP_out.Sources = Sources;
       LFP_out.LFP_Sinks = LFP_Sinks;
       
 end