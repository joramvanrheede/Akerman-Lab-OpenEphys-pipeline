clear all;
testfiledir = 'C:\Users\Alex Vk\Documents\2017\DPHIL\In Vivo\PreProcessed Data\Drive\2020_01_27_a2';
matfiles = dir(fullfile(testfiledir, '*.mat'));
nfiles = length(matfiles);
for i = 1 : nfiles
    disp(['Loading...' matfiles(i).name]);
   load(fullfile(matfiles(i).folder, matfiles(i).name));
    b(i) = AVK_drive_opto_only(ephys_data);
    a(i) = AVK_drive_analysis(ephys_data,(1:8));
    binned_rates(i,:) = nanmean(a(i).binned_rate_by_trial,1);
    close all;
end;


