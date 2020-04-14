clear all; 
close all;
str_Directory = ['D:\Multi_unit Coalated\2_Cortical\L2_3_evoked\_max']; %['D:\Multi_unit Coalated\1_POM\LED Powers\2019_11_26_A1'];
matfiles = dir(fullfile(str_Directory, '*.mat'));
nfiles = length(matfiles);
x = (-50:199)/1000;
ignore_list = []
PSTH_fig        = figure;
set(gcf,'Units','normalized','Position',[.3 .1 .4 .5],'PaperPositionMode','auto')

for i = 1 : nfiles
    if ~ismember(i,ignore_list);
    disp(['Loading...' matfiles(i).name]);
    load(fullfile(matfiles(i).folder, matfiles(i).name));
    L2_3_Data(i) = L2_3_Evoked;
    spike_delta(i) = nanmean(L2_3_Evoked.spike_rate-L2_3_Evoked.spont_rate);
    spont_rate(i) = nanmean(L2_3_Evoked.spont_rate);
    
    mid_spike_delta(i) = nanmean(L2_3_Evoked.mid_spike_rate-L2_3_Evoked.mid_spont_rate);
    mid_spont_rate(i) = nanmean(L2_3_Evoked.mid_spont_rate);
    
    late_spike_delta(i) = nanmean(L2_3_Evoked.late_spike_rate-L2_3_Evoked.late_spont_rate);
    late_spont_rate(i) = nanmean(L2_3_Evoked.late_spont_rate);
    
    L5(i) = L2_3_Evoked.Layer5_Response_Calculated;
    PSTH(i,:) = L2_3_Evoked.psth_L2_3;%- L2_3_Evoked.psth_Spont;
    PSTH(i,50:54) = NaN*ones(1,5);
    hold on;
    plot(x,PSTH(i,:),'LineWidth',1)
    else
        PSTH(i,:) = NaN*x;
    end;
    refline(0,0);
end;
xlabel('Time (s)')
ylabel('Spike count in time bin')
fixplot

PSTH_mean          = figure;
set(gcf,'Units','normalized','Position',[.3 .1 .4 .5],'PaperPositionMode','auto')
PSTH = nanmean(PSTH,1);
baseline = nanmean(PSTH(1:40));
PSTH = PSTH-baseline;
plot(x,PSTH,'k','LineWidth',2);
refline(0,0);

figure();
subplot(2,2,1);
scatter(L5,spike_delta);

subplot(2,2,2);
scatter(L5,mid_spike_delta);

subplot(2,2,3);
scatter(L5,late_spike_delta);

legend([{'POM'};{'Cortical'}])
output = [spike_delta;mid_spike_delta;late_spike_delta]

disp(['Saving figures to...' str_Directory '\']);
    saveas(PSTH_fig,[str_Directory '\' 'Cortical_PSTHs.fig']);
saveas(PSTH_fig,[str_Directory '\' 'Cortical_PSTHs.fig']);

save([str_Directory '\' 'Cortical_PSTH.mat'],'PSTH');
