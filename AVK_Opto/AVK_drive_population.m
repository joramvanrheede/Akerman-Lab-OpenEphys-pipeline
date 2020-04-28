clear all;
testfiledir = 'D:\Multi_unit Coalated\2_Cortical\Drive';
exp_dir = ['2019_11_25_A1']
output_dir = 'D:\Multi_unit Coalated\2_Cortical\Drive_analysed'

matfiles = dir(fullfile(testfiledir, exp_dir, '*.mat'));

nfiles = length(matfiles);
for i = 1 : nfiles
    disp(['Loading...' matfiles(i).name]);
   load(fullfile(matfiles(i).folder, matfiles(i).name));
   Drive_Exps(i) = AVK_drive_analysis(ephys_data);
 
 end;

if ~exist([output_dir], 'dir')
       mkdir([output_dir])
end
disp('Saving outpus');
save([output_dir '\' exp_dir 'Drive_analysed.mat'],'Drive_Exps');
disp('Saved');