function[CSD_out,LFP_Avg_out,sinks,sources] = CSD_Analysis(window_edges,CSD_window,si,channels,cond_data,stim_type,OutputFn,Figout);

%%CSD_Anaylsis(window_edges,CSD_window,si,channels,cond_data,stim_type)
%
%Performs Current Source Analysis on ephys condition Cond_data
% Inputs : 
%           Window edges : time period in S before and after stimulus to
%           keep.
%           CSD_window : time windown in S after stimulus to analyse for
%           sinks
%           si : sampling interval in s
%           cond_data : condition data from an ephys data
%           stim type : either 'LED', 'Whisk' or 'LED_WWindow'

%% Set timing from seconds to sample points
    whisk_time = cond_data.whisk_onset/si; % whisk onset in sampled points; %1
    LED_onset_time = cond_data.LED_onset/si; % 
    LED_duration_time = cond_data.LED_duration/si;
    LED_offset_time = LED_onset_time+LED_duration_time;
    Chan_string = ['Chans ' num2str(channels(1)) '_' num2str(channels(end))];
    
    %% Set traces for whisk and LFP
    whisk_traces = cond_data.LFP_trace(channels,:,whisk_time-window_edges:whisk_time+window_edges); 
    LED_traces = cond_data.LFP_trace(channels,:,LED_onset_time-window_edges:LED_offset_time+window_edges);
    
    LED_Whisk_window = cond_data.LFP_trace(channels,:,LED_onset_time-window_edges+500:LED_onset_time+window_edges+500);
    whisk_traces = whisk_traces/1000; % converts from uV to mV
    LED_traces = LED_traces/1000; % converts from uV to MV
    LED_Whisk_window =   LED_Whisk_window/1000;   
clear w CSD sinks sources;

%% calculate CSD
switch stim_type
    case 'LED'
         LFP_traces = LED_traces;
         
           [CSD_out,LFP_Avg_out,sinks,sources,a_window] = Current_Source_Density(LFP_traces,CSD_window);
    case 'Whisk'
           LFP_traces = whisk_traces;
           [CSD_out,LFP_Avg_out,sinks,sources,a_window] = Current_Source_Density(LFP_traces,CSD_window);
    case 'LED_WWindow'
        disp('LED whisk window');
         LFP_traces = LED_Whisk_window;
            [CSD_out,LFP_Avg_out,sinks,sources,a_window] = Current_Source_Density(LFP_traces,CSD_window);
  
end;

%% plot Figures
if Figout == true;
    ha = figure('Name',[stim_type ' ' Chan_string],'NumberTitle','off');
    set(gcf,'Units','normalized','Position',[.2 .1 .4 .8])% left bottom width height
    
    colormap('hot');   % set colormap
    subplot(2,2,1)
    imagesc(LFP_Avg_out);        % draw image and scale colormap to values range
    title('LFP averaged by trial');
    xlabel('time(ms)');
    ylabel('Channel');
    c = colorbar;          % show color scale
    c.Label.String = 'LFP (mV)';
    
    subplot(2,2,2)
    imagesc(CSD_out);        % draw image and scale colormap to values range
    title('CSD');
    c = colorbar;          % show color scale
    xlabel('time(ms)');
    ylabel('Channel');
    c.Label.String = 'CSD uA/m^3';
    hold on
    line_1_x = [a_window(1) a_window(2); a_window(1) a_window(2)];
    line_1_y = [channels(1) channels(end); channels(1) channels(end)];
    plot(line_1_x,line_1_y','b--');
    
    subplot(2,2,3)
    plot(sinks);
    title('Sinks');
    xlabel('Channel');
    ylabel('CSD uA/m^3');
    
    
    subplot(2,2,4)
    plot(sources);
    title('Sources');
    xlabel('Channel');
    ylabel('CSD uA/m^3');

disp('saving figure output...')  
disp([OutputFn stim_type Chan_string '.fig'])
saveas(ha,[OutputFn stim_type Chan_string '.fig']);
saveas(ha,[OutputFn stim_type Chan_string '.png']);
disp('saved');
end;
end


 
 
 
 