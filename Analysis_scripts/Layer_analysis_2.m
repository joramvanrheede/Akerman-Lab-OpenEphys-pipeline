close all

%%
reload  = 0;
%%
if reload == 1
    clear all
load_folder = 'F:\Synced\Synced_individual'
group_folders     	= {'POM';'Cortical'}; 
reload =1;
end;


%
bin_size = 0.005;
spont_win = [0 0.03]
resp_win = [0.006 0.03]; 
% Max and min depth for units to include
max_depth           = -50;  % Depth relative to L4 sink (negative numbers = towards pia)
min_depth           = -300; % Depth relative to L4 sink (negative numbers = towards pia)

L5_max_depth        =  400;   % Depth relative to L4 sink (negative numbers = towards pia)
L5_min_depth        =  100;  % Depth relative to L4 sink (negative numbers = towards pia)

L5_chans            = 5:12; % Layer 5 channels relative to the L4 sink (sink channel + [L5_chans])

%% Script starts : 
if reload == 1;
    for k = 1: numel(group_folders)
    load_fn = [load_folder '\' group_folders{k} '_selected_exps.mat']
    disp(['Loading...' load_fn]);
    load(load_fn);
    groups(k) = selected_exps;
    disp('Done');
    end;
else 
    disp('No reload');    
end;

clear selected_exps;

%% sanity checks on unit behaviour 

group_i = 2;

for k = 1: 7;%numel(groups(group_i).ephys_data(1).conditions)
Unit_behave(k) = All_unit_behaviour(groups,spont_win,resp_win,k);
end;

disp(['Saving' group_folders{group_i} 'unit_behaviour.mat']) 
save([load_folder '\' group_folders{group_i} 'unit_behaviour.mat'],'Unit_behave');
disp('Saved');

 ha = figure('Name','Unit Behaviour','NumberTitle','off');
 set(gcf,'Units','normalized','Position',[.05 .1 0.9 .8]) %left bottom width height

 %{
subplot(3,4,1);
for a = 1 : numel(Unit_behave.groups)
for j = 1 : numel(Unit_behave.groups(a).ephys_data)
for k = 1:size (Unit_behave.groups(a).ephys_data(j).Spont,1);
    spont = Unit_behave.groups(a).ephys_data(j).Spont(k,:);
    whisk = Unit_behave.groups(a).ephys_data(j).Whisk(k,:);
   hold on;
    plot(spont,whisk,'o');
    xlabel('spontaneous spike rate');
    ylabel('whisker evoked spike rate');
end;
end;
end;
%}
subplot(3,3,1)
hold on;
scatter(Unit_behave.all_unit_rate,Unit_behave.all_unit_whisk_rate);
eq = [0 max(Unit_behave.all_unit_whisk_rate)]
plot(eq,eq,'b--','LineWidth', 1);
xlabel('spontaneous spike rate');
ylabel('whisker evoked spike rate');

subplot(3,3,2)
hold on;
a = Unit_behave.all_unit_whisk_rate'./(Unit_behave.all_unit_whisk_rate'+Unit_behave.all_unit_rate);
plot(1,a,'bo');
ylabel('whisk/(whisk+spont rate)');
refline(0,0.5);

subplot(3,3,3)
hold on;
scatter(Unit_behave.all_unit_depth,a);
h = refline(0,0.5);
eq = [400 0;400 1];
plot(eq(:,1),eq(:,2),'b--','LineWidth', 1);
xlabel('Unit "depth" (800-kilosort_depth)');
ylabel('whisk /(whisk+spont rate)');

%{
subplot(3,4,5);
for a = 1 : numel(Unit_behave.groups)
for j = 1 : numel(Unit_behave.groups(a).ephys_data)
for k = 1:size (Unit_behave.groups(a).ephys_data(j).Spont,1);
    spont = Unit_behave.groups(a).ephys_data(j).Spont(k,:);
    opto = Unit_behave.groups(a).ephys_data(j).Opto(k,:);
   hold on;
    plot(spont,opto,'o');
    xlabel('spontaneous spike rate');
    ylabel('LED evoked spike rate');
end;
end;
end;
%}
subplot(3,3,4)
hold on;
scatter(Unit_behave.all_unit_rate,Unit_behave.all_unit_opto_rate);
eq = [0 max(Unit_behave.all_unit_opto_rate)];
plot(eq,eq,'b--','LineWidth', 1);

xlabel('spontaneous spike rate');
ylabel('LED evoked spike rate');

subplot(3,3,5)
hold on;
a = Unit_behave.all_unit_opto_rate'./(Unit_behave.all_unit_opto_rate'+Unit_behave.all_unit_rate);
plot(1,a,'bo');
ylabel('LED /(LED+spont rate)');
refline(0,0.5);

subplot(3,3,6)
hold on;
scatter(Unit_behave.all_unit_depth,a);
h = refline(0,0.5);
eq = [400 0; 400 1];
plot(eq(:,1),eq(:,2),'b--','LineWidth', 1);
xlabel('Unit "depth" (800-kilosort_depth)');
ylabel('LED rate/(LED spike rate)');

%{
subplot(3,4,9);
for a = 1 : numel(Unit_behave.groups)
for j = 1 : numel(Unit_behave.groups(a).ephys_data)
for k = 1:size (Unit_behave.groups(a).ephys_data(j).Spont,1);
    whisk = Unit_behave.groups(a).ephys_data(j).Whisk(k,:);
    opto = Unit_behave.groups(a).ephys_data(j).Opto(k,:);
   hold on;
    plot(whisk,opto,'o');
    xlabel('Whisk spike rate');
    ylabel('LED evoked spike rate');
end;
end;
end;
%}
subplot(3,3,7)
hold on;
scatter(Unit_behave.all_unit_whisk_rate,Unit_behave.all_unit_opto_rate);
eq = [0 max(Unit_behave.all_unit_whisk_rate)];
plot(eq,eq,'b--','LineWidth', 1);

xlabel('Whisk spike rate');
ylabel('LED evoked spike rate');

subplot(3,3,8)
hold on;
a = Unit_behave.all_unit_opto_rate./(Unit_behave.all_unit_opto_rate+Unit_behave.all_unit_whisk_rate);
plot(1,a,'bo');
ylabel('LED/(LED + Whisk rate)');
refline(0,0.5);

subplot(3,3,9)
hold on;
scatter(Unit_behave.all_unit_depth,a);
h = refline(0,0.5);
eq = [400 0; 400 1];
plot(eq(:,1),eq(:,2),'b--','LineWidth', 1);
xlabel('Unit "depth" (800-kilosort_depth)');
ylabel('LED/(LED+Whisk rate)');
%}


n_bins = numel(Unit_behave.depth_edges)-1;
figure ();
hold on;
scatter(Unit_behave.depth_bins,Unit_behave.all_unit_rate);
plot(1:n_bins,Unit_behave.mean_unit_rate,'b-');

scatter(Unit_behave.depth_bins,Unit_behave.all_unit_opto_rate,'r*');
plot(1:n_bins,Unit_behave.mean_opto_rate,'r-');

scatter(Unit_behave.depth_bins,Unit_behave.all_unit_whisk_rate,'g*');
plot(1:n_bins,Unit_behave.mean_whisk_rate,'g-');


figure ();
hold on;
scatter(Unit_behave.depth_bins,Unit_behave.all_unit_isi);
plot([1:n_bins],Unit_behave.mean_unit_isi);

%% Checks channel correlations with unit activity. 
clear chan_corr;
chan_corr = [];
unit_depths = [];
channel_depths = [];
for z = 1 : numel(groups);
for i = 1 : numel(groups(z).ephys_data)

unit_trial_spikes = [];
MUA_trial_spikes = [];

unit_depths = [unit_depths;groups(z).ephys_data(i).unit_depths];

for j = 1 : numel(groups(z).ephys_data(i).conditions);
    
unit_spikes = groups(z).ephys_data(i).conditions(j).unit_spikes; 
MUA_spikes = groups(z).ephys_data(i).conditions(j).multiunit_spikes;
binned_unit_spikes = squeeze(sum(unit_trial_bins(unit_spikes,bin_size,resp_win),3));
binned_MUA_spikes = squeeze(sum(unit_trial_bins(MUA_spikes,bin_size,resp_win),3));
unit_trial_spikes = [unit_trial_spikes binned_unit_spikes];
MUA_trial_spikes = [MUA_trial_spikes binned_MUA_spikes];
end;

clear rho pval u_d m_d;
if ~ isempty(unit_trial_spikes);
[rho,pval] = corr(unit_trial_spikes',MUA_trial_spikes','type','Spearman');
[~,cc] = max(rho,[],2);
m_d = cc*25;
else
    rho = NaN;
    m_d = [];
end;
channel_depths = [channel_depths;m_d];
end;
end;

figure();
plot(unit_depths,channel_depths,'bo');
title('Channel Depth in which activity best correlates to unit activity')
xlabel('unit depths(um)');
ylabel('channel depths(um)');



%}

%{
for k = 1 : numel(selected_exps(1).ephys_data);
    if ~isempty(selected_exps.ephys_data(k).conditions)
          group(k).cond = selected_exps.ephys_data(k).conditions(end);
          group(k).max_sink_chan = selected_exps.ephys_data(k).max_sink_chan;
          group(k).unit_depths  = selected_exps.ephys_data(k).unit_depths;
     end;
end;

%cleans data removing all empty experiments
for k = 1 : numel(group);
    rem(k) = isempty(selected_exps.ephys_data(k).conditions);
     end;
group(rem) = [];
clear rem;



%%
% Get L5 multi_unit activity and L23 units. 
for k  = 1:numel(group)
     these_L5_chans                          = group(k).max_sink_chan + L5_chans; % Finds layer 5 chans
     these_L5_chans(these_L5_chans>32)       = []; 
     group(k).L5_chans                       = these_L5_chans; 
     group(k).max_sink_depth                 = group(k).max_sink_chan*25;
     unit_depths                             = group(k).unit_depths;
     unit_depths                             =  unit_depths(:) - group(k).max_sink_depth;
     group(k).q_L23_depth                    =  unit_depths <= max_depth & unit_depths >= min_depth;
     group(k).q_L5_depth                     =  unit_depths <= L5_max_depth & unit_depths >= L5_min_depth;
     group(k).L23_units                      = sum(group(k).q_L23_depth);
     clear these_L5_chans unit_depth q_L23_depth q_L5 depth;
end;

exp_prep = 4;
L23_spikes = group(exp_prep).cond.unit_spikes(group(exp_prep).q_L23_depth,:,:);
L23_units_binned = unit_trial_bins(L23_spikes,bin_size,resp_win);

L5_MUA_spikes = group(exp_prep).cond.multiunit_spikes(group(exp_prep).L5_chans,:,:);
L5_spikes_binned = unit_trial_bins(L5_MUA_spikes,bin_size,resp_win);
L5_spikes_binned = sum(L5_spikes_binned,1);

a =squeeze(sum(L23_units_binned,3));
b = squeeze(sum(L5_spikes_binned,3));


figure(1);
subplot(2,2,1)
raster_plot(L23_spikes,2,0.2);
subplot(2,2,2)
bar(a);
subplot(2,2,3);
raster_plot(L5_MUA_spikes,2,0.2);
subplot(2,2,4)
bar(b);

figure(); 
hold on;
for k = 1: size(a,1)
scatter(b,a(k,:));
end;
[h,p] = corr(b',a','type','Spearman','rows','pairwise');
%}