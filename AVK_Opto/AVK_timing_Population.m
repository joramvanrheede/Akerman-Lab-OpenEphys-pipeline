clear all;
testfiledir = 'D:\Multi_unit Coalated\2_Cortical\Timing';
exp_dir = ['2019_05_21']
output_dir = 'D:\Multi_unit Coalated\2_Cortical\Timing_analysed'

matfiles = dir(fullfile(testfiledir, exp_dir, '*.mat'));

nfiles = length(matfiles);
for i = 1 : nfiles
    disp(['Loading...' matfiles(i).name]);
   load(fullfile(matfiles(i).folder, matfiles(i).name));
   Timing_Exps(i) = AVK_timing_Analyse(ephys_data);
   close all;
 end;

if ~exist([output_dir], 'dir')
       mkdir([output_dir])
end
disp('Saving outpus');
save([output_dir '\' exp_dir 'timing_analysed.mat'],'Timing_Exps');
disp('Saved');