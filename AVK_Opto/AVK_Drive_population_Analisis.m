clear all; close all; clc;

testfiledir = 'D:\Multi_unit Coalated\1_POM\Drive_analysed';

matfiles = dir(fullfile(testfiledir, '*.mat'));
 
Global_Delta_T = [];
Layer_5_Target = 1000;
nfiles = length(matfiles);
for i = 1 : nfiles
    disp(['Loading...' matfiles(i).name]);
   load(fullfile(matfiles(i).folder, matfiles(i).name));
   L5Resp = [Drive_Exps(:).h_L5_resp];
   index = [Drive_Exps(:).Mean_L5_Opto];
   
   for j = 1 : numel(L5Resp)
       if L5Resp(j) ==0;
         index(j) = NaN;
       end;
   end;
   [~, index] = min(abs(index(:)-Layer_5_Target));
   if L5Resp(index)==0 
   disp('Warning not signifcant')
   end;
   Data_Set(i) = Drive_Exps(index);
   L5_out(i) = Drive_Exps(index).Mean_L5_Opto;
   Binned_rate(i,:) = nanmean(Drive_Exps(index).binned_rate_by_trial,1);
   Peak_rate(i,:) = nanmean(Drive_Exps(index).peak_rate_by_trial,1);
   Peak_time_by_trial(i,:)=nanmean(Drive_Exps(index).peak_time_by_trial,1);
end;

Binned_rate = Binned_rate(:,1)-Binned_rate(:,2);
Peak_rate = Peak_rate(:,1)-Peak_rate(:,2);
Peak_time_by_trial = Peak_time_by_trial(:,1)-Peak_time_by_trial(:,2);
  