% LTP_RWS_vs_control
clearvars
close all

data_dir        = '/Volumes/Akermanlab/Joram/Cortex_20x100Hz_Plasticity';

expt_folders    = dir(data_dir);


channels        = [15:26];


% Response assessment windows
spont_win       = [-1 -0.100]; % window at the start of the trial for assessing spontaneous activity
whisk_win       = [0 0.025]; % Window for assessing post-whisk spiking response
LFP_win         = [0 0.200]; % Window for assessing peak-to-peak


plotwin         = [-0.01 0.1]; % time window for PSTHs to plot

psth_smoothing  = 5; % smoothing window for PSTHs when plotted


qprint          = 1;
print_folder    = '/Users/Joram/Dropbox/Akerman Postdoc/Figures/2019_Lab_meeting';
print_postfix  	= '100Hz_L5';

%% Manually set these variables for now... --> hacks!

% 8Hz vals: 75, 2400, 100
% 100Hz vals: 30, 800, 75

max_n_trials          	= 30; % limit to N trials
post_start_time        	= 800;
tdiff                   = 75;


%%

% To do:
% Responsiveness criterion
%   - Relative to spontaneous
%   - Absolute
% Subtract spontaneous activity
% Remove non-responsive channels
% Deal with control files that have extra data with a new induction
% protocol
% Work out when to take the sum of spike numbers and when to take mean
% Do exact timing of data points


%% 
expt_names                  = {expt_folders.name};

qremove                     = ismember(expt_names,{'.','..','.DS_Store'}); 
expt_names(qremove)         = [];

test_pre_spikecount         = [];
test_post_spikecount        = [];
control_pre_spikecount      = [];
control_post_spikecount     = [];

test_pre_LFP_P2P            = [];
test_post_LFP_P2P           = [];
control_pre_LFP_P2P         = [];
control_post_LFP_P2P        = [];

control_pre_trial_spikes    = {};
control_post_trial_spikes   = {};
test_pre_trial_spikes       = {};
test_post_trial_spikes      = {};

control_pre_trial_LFP       = {};
control_post_trial_LFP      = {};
test_pre_trial_LFP          = {};
test_post_trial_LFP         = {};

test_pre_time               = datetime;
test_post_time              = datetime;
control_pre_time            = datetime;
control_post_time           = datetime;

figure(1)
set(gcf,'Visible','off')
figure(2)
set(gcf,'Visible','off')
for a = 1:length(expt_names)
    this_expt   = expt_names{a};
    
    if ~ismember(this_expt,{'Test_1_pre_0' 'Test_1_pre_1' 'Test_1_post_0' 'Test_1_post_1'})
        continue
    end
    
    this_expt_folder = fullfile(data_dir, this_expt);
    
    date_folders    = dir(this_expt_folder);
    dates           = {date_folders.name};
    
    qremove         = ismember(dates,{'.','..','.DS_Store'});
    dates(qremove)  = [];
    
    for b = 1:length(dates)
        this_date   = dates{b};
        
        this_date_folder = [this_expt_folder filesep this_date];
        
        data_files  = dir([this_date_folder filesep '*.mat']);
        data_files  = {data_files.name};
        
        plot_psth_stack             = []; % add PSTHS for experiments of the same type to make mean later
        plot_LFP_stack              = []; % stack LFPs for experiments of the same type to make mean later
        
        for c = 1:length(data_files)
            
            %% Data loading
            this_data_file = data_files{c};
            full_data_file = [this_date_folder filesep this_data_file];
            
            disp(['Loading ' full_data_file '...']);
            load(full_data_file)
            
            %  add psth info to ephys_data
            ephys_data      = ephys_data_psths(ephys_data,0.0001, 0.010);
           
            % get protocol time
            protocol_time   = get_protocol_time(ephys_data.data_folder);
            
            %% Whisker response
            
            time_bins       = ephys_data.conditions(1).psth_bins;
            whisk_time      = ephys_data.conditions(1).whisk_onset;
            whisk_bins      = time_bins - whisk_time;
            
            q_whisk_win     = whisk_bins >= whisk_win(1) & whisk_bins < whisk_win(2);
            q_spont_win     = whisk_bins >= spont_win(1) & whisk_bins < spont_win(2);
            
            cond_psths      = ephys_data.conditions(1).trial_psths;
            
            cond_psth_mean_by_channel   = squeeze(mean(cond_psths(channels,:,:),1)); % number of spikes, mean over channels
            cond_psth_mean_by_trial     = mean(cond_psth_mean_by_channel,1); % number of spikes, mean over channels and trials
            
            % Determine whisker response (summed spikes in whisk bin, spikes/channel/trial) 
            whisk_spike          = sum(cond_psth_mean_by_trial(q_whisk_win));
            
            % Whisker response for each individual trial
            whisk_spike_by_trial = sum(cond_psth_mean_by_channel(:,q_whisk_win),2);
            
            % Spontaneous spiking response
            spont_spikes        = sum(cond_psth_mean_by_trial(q_spont_win));
            
            whisk_win_size      = whisk_win(2) - whisk_win(1);
            spont_win_size      = spont_win(2) - spont_win(1);
            win_ratio           = spont_win_size / whisk_win_size;
            
            corr_spont_spikes   = spont_spikes / win_ratio;
            
            % add psth to stack for this experiment type
            plot_psth_stack     = [plot_psth_stack; cond_psth_mean_by_trial];
            
            % Spont spiking per trial
            
            spont_spike_by_trial = sum(cond_psth_mean_by_channel(:,q_spont_win),2);
            spont_spike_corr     = spont_spike_by_trial / win_ratio;
            spont_spike_std      = std(spont_spike_corr) * sqrt(win_ratio);
            
            % Whisk spiking per trial
            
            whisk_spike_by_trial = sum(cond_psth_mean_by_channel(:,q_whisk_win),2);
            whisk_spike_corr     = whisk_spike_by_trial / win_ratio;
            whisk_spike_std      = std(whisk_spike_corr);
            
            % Correction for baseline activity happens here
            whisk_spike       	= whisk_spike - corr_spont_spikes;
            
            %% LFP response
            
            LFP_traces          = ephys_data.conditions.expt_LFPs(channels,1:3000);
            LFP_timebins        = [1:size(LFP_traces,2)]/1000;
            LFP_timebins        = LFP_timebins - whisk_time;
            mean_LFP_trace      = mean(LFP_traces,1);
            
            mean_LFP_trace      = notch_filt(mean_LFP_trace,1000,50);
            
            q_LFP_resp          = LFP_timebins >= LFP_win(1) & LFP_timebins < LFP_win(2);
            
            LFP_N1              = min(mean_LFP_trace(q_LFP_resp));
            LFP_P1              = max(mean_LFP_trace(q_LFP_resp));
            LFP_P2P             = LFP_P1 - LFP_N1;
            
            LFP_traces_by_trial = squeeze(mean(ephys_data.conditions(1).LFP_trace,1));
            
            LFP_P2P_by_trial    = NaN(size(LFP_traces_by_trial,1),1);
            for d = 1:size(LFP_traces_by_trial,1)
                this_LFP_trace          = notch_filt(LFP_traces_by_trial(d,:),1000,50);
                this_LFP_max            = max(this_LFP_trace(q_LFP_resp));
                this_LFP_min            = min(this_LFP_trace(q_LFP_resp));
                LFP_P2P_by_trial(d)     = this_LFP_max - this_LFP_min;
            end
            
            plot_LFP_stack      = [plot_LFP_stack; mean_LFP_trace];
            
            whisk_spike_by_trial    = whisk_spike_by_trial(1:max_n_trials);
            LFP_P2P_by_trial        = LFP_P2P_by_trial(1:max_n_trials);
            
            if any(strfind(data_files{c},'pre_0'))
                figure(1)
                subplot(2,2,3)
                title('Control condition - before blank stimulus')
                figure(2)
                subplot(2,2,3)
                title('Control condition - before blank stimulus')
                
                control_pre_spikecount(b,c)  	= whisk_spike;
                control_pre_LFP_P2P(b,c)       	= LFP_P2P;
                control_pre_time(b,c)           = protocol_time;
                
                control_pre_trial_spikes(b,c)  	= {whisk_spike_by_trial};
                control_pre_trial_LFP(b,c)  	= {LFP_P2P_by_trial};
                
                
            elseif any(strfind(data_files{c},'post_0'))
                figure(1)
                subplot(2,2,4)
                title('Control condition - after blank stimulus')
                figure(2)
                subplot(2,2,4)
                title('Control condition - after blank stimulus')
                
                control_post_spikecount(b,c)	= whisk_spike;
                control_post_LFP_P2P(b,c)     	= LFP_P2P;
                control_post_time(b,c)      	= protocol_time;
                
                control_post_trial_spikes(b,c) 	= {whisk_spike_by_trial};
                control_post_trial_LFP(b,c) 	= {LFP_P2P_by_trial};
                
                
            elseif any(strfind(data_files{c},'pre_1'))
                figure(1)
                subplot(2,2,1)
                title('Test condition - before induction')
                figure(2)
                subplot(2,2,1)
                title('Test condition - before induction')
                
                test_pre_spikecount(b,c)        = whisk_spike;
                test_pre_LFP_P2P(b,c)           = LFP_P2P;
                test_pre_time(b,c)              = protocol_time;
                
                test_pre_trial_spikes(b,c)    	= {whisk_spike_by_trial};
                test_pre_trial_LFP(b,c)         = {LFP_P2P_by_trial};
                
                
            elseif any(strfind(data_files{c},'post_1'))
                figure(1)
                subplot(2,2,2)
                title('Test condition - after induction')
                figure(2)
                subplot(2,2,2)
                title('Test condition - after induction')
                
                test_post_spikecount(b,c)       = whisk_spike;
                test_post_LFP_P2P(b,c)          = LFP_P2P;
                test_post_time(b,c)             = protocol_time;
                
                test_post_trial_spikes(b,c)   	= {whisk_spike_by_trial};
                test_post_trial_LFP(b,c)        = {LFP_P2P_by_trial};
                
            else
                warning('Unsure of file type - no ''pre_1'' or ''post_0'' designation found in filename')
            end
            
            
            
        end
        
        
        figure(1)
        plot(whisk_bins,smooth(mean(plot_psth_stack,1),psth_smoothing))
        xlim(plotwin)
        hold on
        
        figure(2)
        plot(LFP_timebins,mean(plot_LFP_stack,1))
        xlim(plotwin)
        hold on
        
        
        
    end
    
end

% some figure things
for a = 1:2
    figure(a)
    set(gcf,'Visible','on')
    set(gcf,'Units','Normalized')
    set(gcf,'Position',[.2 .2 .6 .7])
    % Set all y axes to the same range (based on the largest range)
    plotaxes    = get(gcf,'Children');
    maxy        = cellfun(@max,get(plotaxes,'Ylim'));
    miny        = cellfun(@min,get(plotaxes,'Ylim'));
    
    maxy        = max(maxy);
    miny        = min(miny);
    
    set(plotaxes,'YLim',[miny maxy]); % add margin 
    
    hold off

end


%% Combine pre and post responses in matrix for plotting purposes
test_spike_counts       = [test_pre_spikecount test_post_spikecount];
control_spike_counts    = [control_pre_spikecount control_post_spikecount];

test_LFP_P2Ps           = [test_pre_LFP_P2P test_post_LFP_P2P];
control_LFP_P2Ps        = [control_pre_LFP_P2P control_post_LFP_P2P];

test_times              = [test_pre_time test_post_time];
control_times           = [control_pre_time control_post_time];


test_mean_spike_resp    = mean(test_spike_counts,2);
q_test_spike_resp       = test_mean_spike_resp > 0.5;

control_mean_spike_resp = mean(control_spike_counts,2);
q_control_spike_resp  	= control_mean_spike_resp > 0.5;


test_spike_counts       = test_spike_counts(q_test_spike_resp,:);
test_times              = test_times(q_test_spike_resp,:);
test_LFP_P2Ps           = test_LFP_P2Ps(q_test_spike_resp,:);

control_spike_counts    = control_spike_counts(q_control_spike_resp,:);
control_times           = control_times(q_control_spike_resp,:);
control_LFP_P2Ps       	= control_LFP_P2Ps(q_control_spike_resp,:);

% 
% time loop
test_minutes            = NaN(size(test_times));
for a = 1:size(test_minutes,2)
	if a == 1
        test_minutes(:,a) = 0;
        continue
    end
    test_minutes(:,a) = minutes(test_times(:,a) - test_times(:,1));
end

% time loop
control_minutes            = NaN(size(control_times));
for a = 1:size(control_minutes,2)
	if a == 1
        control_minutes(:,a) = 0;
        continue
    end
    control_minutes(:,a) = minutes(control_times(:,a) - control_times(:,1));
end



%% Spiking figures
figure
subplot(2,1,1)
plot(test_minutes',test_spike_counts','k.-','LineWidth',2,'MarkerSize',20)
title('Induction spike counts')
xlimz = xlim;
xlim([-5 xlimz(2)])

subplot(2,1,2)
plot(control_minutes',control_spike_counts','k.-','LineWidth',2,'MarkerSize',20)
title('Control spike counts')
xlimz = xlim;
xlim([-5 xlimz(2)])

%% long term plot figure


%% Spikes
% control pre
control_pre_spikes      = cell2mat(control_pre_trial_spikes')';
control_pre_means       = mean(control_pre_spikes,2);
control_pre_points      = [];
for a = 1:(size(control_pre_spikes,2)/5)
    inds    = (a-1)*5+[1:5];
    control_pre_points(:,a)     = mean(control_pre_spikes(:,inds),2)./control_pre_means;
end

control_pre_serr            = serr(control_pre_points);
control_pre_time_points     = [1:size(control_pre_points,2)]*tdiff;

% control post
control_post_spikes         = cell2mat(control_post_trial_spikes')';
control_post_points         = [];
for a = 1:(size(control_post_spikes,2)/5)
    inds    = (a-1)*5+[1:5];
    control_post_points(:,a)    = mean(control_post_spikes(:,inds),2)./control_pre_means;
end

control_post_serr           = serr(control_post_points);
control_post_time_points  	= post_start_time + [1:size(control_post_points,2)]*tdiff;

% merge pre and post
control_all_points          = [control_pre_points control_post_points];
control_all_serrs           = [control_pre_serr control_post_serr];
control_time_vect           = [control_pre_time_points control_post_time_points]/60;

figure
set(gcf,'Units','Normalized')
set(gcf,'Position',[0.1 0.2 0.8 0.7])
subplot(2,1,1)
errorbar(control_time_vect,mean(control_all_points), control_all_serrs,'k.','MarkerSize',20);
title('Spiking response - Control')
xlabel('Time (minutes)')
ylabel('Normalised spike count')
set(gca,'FontSize',14,'FontName','Arial','LineWidth',2)

%
% test pre
test_pre_spikes      = cell2mat(test_pre_trial_spikes')';
test_pre_means       = mean(test_pre_spikes,2);  
test_pre_points      = [];
for a = 1:(size(test_pre_spikes,2)/5)
    inds    = (a-1)*5+[1:5];
    test_pre_points(:,a)     = mean(test_pre_spikes(:,inds),2)./test_pre_means;
end

test_pre_serr            = serr(test_pre_points);
test_pre_time_points     = [1:size(test_pre_points,2)]*tdiff;

% test post
test_post_spikes         = cell2mat(test_post_trial_spikes')';
test_post_points         = [];
for a = 1:(size(test_post_spikes,2)/5)
    inds    = (a-1)*5+[1:5];
    test_post_points(:,a)    = mean(test_post_spikes(:,inds),2)./test_pre_means;
end

test_post_serr        	= serr(test_post_points);
test_post_time_points 	= post_start_time + [1:size(test_post_points,2)]*tdiff;

% merge pre and post
test_all_points      	= [test_pre_points test_post_points];
test_all_serrs      	= [test_pre_serr test_post_serr];
test_time_vect      	= [test_pre_time_points test_post_time_points]/60;

subplot(2,1,2)
errorbar(test_time_vect,mean(test_all_points), test_all_serrs,'k.','MarkerSize',20);
title('Spiking response - Induction')
xlabel('Time (minutes)')
ylabel('Normalised spike count')
set(gca,'FontSize',14,'FontName','Arial','LineWidth',2)

% Both in the same figure

figure
set(gcf,'Units','Normalized')
set(gcf,'Position',[0.1 0.3 0.8 0.6])
errorbar(control_time_vect,mean(control_all_points), control_all_serrs,'k.','MarkerSize',40,'MarkerFaceColor',[.6 .6 .6],'MarkerEdgeColor',[.6 .6 .6],'LineWidth',2,'Color',[.6 .6 .6]);
hold on
errorbar(test_time_vect,mean(test_all_points), test_all_serrs,'k.','MarkerSize',40,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'LineWidth',2);
% xlim([0 95])
title('Spiking response - Control vs LTP induction')
xlabel('Time (minutes)')
ylabel('Normalised spike count')
legend({'Control' 'Induction'},'FontSize',18)
set(gca,'FontSize',18,'FontName','Arial','LineWidth',2)


%% LFP

% pre_data, post_data, trial_interval, protocol_interval

% control pre
control_pre_LFP         = cell2mat(control_pre_trial_LFP')';
control_pre_means       = mean(control_pre_LFP,2);
control_pre_points      = [];
for a = 1:(size(control_pre_LFP,2)/5)
    inds    = (a-1)*5+[1:5];
    control_pre_points(:,a)     = mean(control_pre_LFP(:,inds),2)./control_pre_means;
end

control_pre_serr            = serr(control_pre_points);
control_pre_time_points     = [1:size(control_pre_points,2)]*tdiff;

% control post
control_post_LFP            = cell2mat(control_post_trial_LFP')';
control_post_points         = [];
for a = 1:(size(control_post_LFP,2)/5)
    inds    = (a-1)*5+[1:5];
    control_post_points(:,a)	= mean(control_post_LFP(:,inds),2)./control_pre_means;
end

control_post_serr           = serr(control_post_points);
control_post_time_points  	= post_start_time + [1:size(control_post_points,2)]*tdiff;

% merge pre and post
control_all_points          = [control_pre_points control_post_points];
control_all_serrs           = [control_pre_serr control_post_serr];
control_time_vect           = [control_pre_time_points control_post_time_points]/60;

figure
set(gcf,'Units','Normalized')
set(gcf,'Position',[0.1 0.2 0.8 0.7])
subplot(2,1,1)
errorbar(control_time_vect,mean(control_all_points), control_all_serrs,'k.','MarkerSize',20);
title('LFP size - Control')
xlabel('Time (minutes)')
ylabel('Normalised LFP peak-to-peak')
set(gca,'FontSize',14,'FontName','Arial','LineWidth',2)

%
% test pre
test_pre_LFP      = cell2mat(test_pre_trial_LFP')';
test_pre_means       = mean(test_pre_LFP,2);  
test_pre_points      = [];
for a = 1:(size(test_pre_LFP,2)/5)
    inds    = (a-1)*5+[1:5];
    test_pre_points(:,a)     = mean(test_pre_LFP(:,inds),2)./test_pre_means;
end

test_pre_serr            = serr(test_pre_points);
test_pre_time_points     = [1:size(test_pre_points,2)]*tdiff;

% test post
test_post_LFP         = cell2mat(test_post_trial_LFP')';
test_post_points         = [];
for a = 1:(size(test_post_LFP,2)/5)
    inds    = (a-1)*5+[1:5];
    test_post_points(:,a)    = mean(test_post_LFP(:,inds),2)./test_pre_means;
end

test_post_serr        	= serr(test_post_points);
test_post_time_points 	= post_start_time + [1:size(test_post_points,2)]*tdiff;

% merge pre and post
test_all_points      	= [test_pre_points test_post_points];
test_all_serrs      	= [test_pre_serr test_post_serr];
test_time_vect      	= [test_pre_time_points test_post_time_points]/60;

subplot(2,1,2)
errorbar(test_time_vect,mean(test_all_points), test_all_serrs,'k.','MarkerSize',20);
title('LFP size - Induction')
xlabel('Time (minutes)')
ylabel('Normalised LFP peak-to-peak')
set(gca,'FontSize',14,'FontName','Arial','LineWidth',2)

% Both in the same figure

figure
set(gcf,'Units','Normalized')
set(gcf,'Position',[0.1 0.3 0.8 0.6])
errorbar(control_time_vect,mean(control_all_points), control_all_serrs,'k.','MarkerSize',40,'MarkerFaceColor',[.6 .6 .6],'MarkerEdgeColor',[.6 .6 .6],'LineWidth',2,'Color',[.6 .6 .6]);
hold on
errorbar(test_time_vect,mean(test_all_points), test_all_serrs,'k.','MarkerSize',40,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'LineWidth',2);
% xlim([0 95])
title('LFP size - Control vs LTP induction')
xlabel('Time (minutes)')
ylabel('Normalised LFP peak-to-peak')
legend({'Control' 'Induction'},'FontSize',18)
set(gca,'FontSize',18,'FontName','Arial','LineWidth',2)
ylimz = ylim;
ylim([0 ylimz(2)])


%% LFP figures

figure
subplot(2,1,1)
plot(test_minutes',test_LFP_P2Ps','k.-','LineWidth',2,'MarkerSize',20)
title('Induction LFP P2P')
xlimz = xlim;
xlim([-5 xlimz(2)])

subplot(2,1,2)
plot(control_minutes',control_LFP_P2Ps','k.-','LineWidth',2,'MarkerSize',20)
title('Control LFP P2P')
xlimz = xlim;
xlim([-5 xlimz(2)])

for a = [3 4 6 8]
    figure(a)
    %set(gcf,'Visible','on')
    %set(gcf,'Units','Normalized')
    %set(gcf,'Position',[.2 .2 .6 .7])
    % Set all y axes to the same range (based on the largest range)
    plotaxes    = get(gcf,'Children');
    maxy        = cellfun(@max,get(plotaxes,'Ylim'));
    miny        = cellfun(@min,get(plotaxes,'Ylim'));
    
    maxy        = max(maxy);
    miny        = 0;
    
    set(plotaxes,'YLim',[miny maxy]); % add margin 
    
    hold off

end
% 
if qprint
    fig1_title  = ['PSTHs_2x2_all_recs' '_' print_postfix];
    print(1,[print_folder filesep fig1_title '.png'],'-dpng','-r300')
    print(1,[print_folder filesep fig1_title '.eps'],'-depsc','-r300')
    
    fig2_title  = ['LFPs_2x2_all_recs' '_' print_postfix];
    print(2,[print_folder filesep fig2_title '.png'],'-dpng','-r300')
    print(2,[print_folder filesep fig2_title '.eps'],'-depsc','-r300')
    
    fig3_title  = ['Spikes_ind_vs_control_all_' print_postfix];
    print(3,[print_folder filesep fig3_title '.png'],'-dpng','-r300')
    print(3,[print_folder filesep fig3_title '.eps'],'-depsc','-r300')
    
    fig4_title  = ['Spike_timeseries_separate_' print_postfix];
    print(4,[print_folder filesep fig4_title '.png'],'-dpng','-r300')
    print(4,[print_folder filesep fig4_title '.eps'],'-depsc','-r300')
    
    fig5_title  = ['Spike_timeseries_joint_' print_postfix];
    print(5,[print_folder filesep fig5_title '.png'],'-dpng','-r300')
    print(5,[print_folder filesep fig5_title '.eps'],'-depsc','-r300')
    
    fig6_title  = ['LFP_timeseries_separate_' print_postfix];
    print(6,[print_folder filesep fig6_title '.png'],'-dpng','-r300')
    print(6,[print_folder filesep fig6_title '.eps'],'-depsc','-r300')
   
    fig7_title  = ['LFP_timeseries_joint_' print_postfix];
    print(7,[print_folder filesep fig7_title '.png'],'-dpng','-r300')
    print(7,[print_folder filesep fig7_title '.eps'],'-depsc','-r300')
    
    fig8_title  = ['LFP_ind_vs_control_all_' print_postfix];
    print(8,[print_folder filesep fig8_title '.png'],'-dpng','-r300')
    print(8,[print_folder filesep fig8_title '.eps'],'-depsc','-r300')
    
end

