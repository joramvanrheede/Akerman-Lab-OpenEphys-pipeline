clear all; close all; clc;

testfiledir = 'D:\Multi_unit Coalated\1_POM\Timing_analysed';

matfiles = dir(fullfile(testfiledir, '*.mat'));
 
Global_Delta_T = [];
Layer_5_Target = 25;
nfiles = length(matfiles);
for i = 1 : nfiles
    disp(['Loading...' matfiles(i).name]);
   load(fullfile(matfiles(i).folder, matfiles(i).name));
   index = [Timing_Exps(:).Mean_L5_Opto_delta];
   [~, index] = min(abs(index(:)-Layer_5_Target));
   Global_Delta_T = [Global_Delta_T Timing_Exps(index).delta_t];
   Data_Set(i) = Timing_Exps(index);
   L5_out(i) = Timing_Exps(index).Mean_L5_Opto_delta;
  end;

  Global_Delta_T = unique(Global_Delta_T);
  
  Binned_Rates = NaN*ones(numel(Global_Delta_T),nfiles);
  Peak_Rates = NaN*ones(numel(Global_Delta_T),nfiles);
  Peak_times = NaN*ones(numel(Global_Delta_T),nfiles);
  
  for k = 1: numel(Data_Set)
      for j = 1 : numel(Data_Set(k).delta_t)
      t_index(j) = find(Data_Set(k).delta_t(j) == Global_Delta_T);
      Binned_Rates(t_index(j),k)= Data_Set(k).mean_binned_spike_rate(j);
      Peak_Rates(t_index(j),k)= Data_Set(k).mean_peak_spike_rate(j);
      Peak_times(t_index(j),k)= Data_Set(k).mean_peak_spike_time(j);    
      end;
  Norm_Binned_Rates(:,k) = Binned_Rates(:,k)./Binned_Rates(13,k)
  Norm_Peak_Rates(:,k) = Peak_Rates(:,k)./Peak_Rates(13,k)
  Norm_Peak_times(:,k) = Peak_times(:,k)./Peak_times(13,k) 
  end;
 
  
  figure();
  scatter(Global_Delta_T,Binned_Rates(:,4));
  disp('saving');
  save(['D:\Multi_unit Coalated\1_POM\Timing_max'],'Data_Set','Binned_Rates','Peak_Rates','Peak_times','Global_Delta_T');