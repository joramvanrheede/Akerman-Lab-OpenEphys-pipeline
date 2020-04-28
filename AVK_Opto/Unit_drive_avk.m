%clear all
save_dir = 'E:\Sorted\POm\Drive\'
depths = ephys_data.unit_depths;
ix = find(depths >300,1,'first');

for k = 1:ix-1
  a(k) =  drive_plots(ephys_data,k)
end;

%save([save_dir '2019_02_01_6' '.mat'],'a');
close all;
Binned_light = [];
Binned_Whisk = [];
for k = 1: size(a,2)
    Binned_light = [Binned_light a(k).binned_rate_by_trial(:,1)];
    Binned_Whisk = [Binned_Whisk a(k).binned_rate_by_trial(:,2)];
end;