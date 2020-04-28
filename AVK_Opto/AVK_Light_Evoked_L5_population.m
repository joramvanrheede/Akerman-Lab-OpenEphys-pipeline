clear all; 
close all;
str_Directory = ['D:\Multi_unit Coalated\1_POM\Layer_5_pulse_response\']; %['D:\Multi_unit Coalated\1_POM\LED Powers\2019_11_26_A1'];
matfiles = dir(fullfile(str_Directory, '*.mat'));
nfiles = length(matfiles);

PoP_Median_responses = [];
PoP_Mean_responses = [];

for i = 1 : nfiles
    disp(['Loading...' matfiles(i).name]);
    load(fullfile(matfiles(i).folder, matfiles(i).name));
    Animal_Median_responses = [];
    Animal_Mean_responses = [];
    for k = 1 : numel(File_results)
        Animal_Median_responses = [Animal_Median_responses File_results(k).Median_Light_Delta(File_results(k).Diff_From_Spont(:,:)==1)];
        Animal_Mean_responses = [Animal_Mean_responses File_results(k).Mean_Light_Delta(File_results(k).Diff_From_Spont(:,:)==1)];
    end;
    if ~isempty(Animal_Mean_responses)
    PoP_Max_Mean_response(i) = max(Animal_Mean_responses);
    PoP_Max_Median_response(i) = max(Animal_Median_responses);
    PoP_Min_Mean_response(i) = min(Animal_Mean_responses);
    PoP_Min_Median_response(i) = min(Animal_Median_responses);
    else
     PoP_Max_Mean_response(i) = NaN;
    PoP_Max_Median_response(i) = NaN;
    PoP_Min_Mean_response(i) = NaN;
    PoP_Min_Median_response(i) = NaN;
    
    end;
    
    PoP_Median_responses = [PoP_Median_responses Animal_Median_responses];
    PoP_Mean_responses = [PoP_Mean_responses Animal_Mean_responses];
end;

