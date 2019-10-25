function protocol_time = get_protocol_time(folder_name)
%% Takes the name of the experiment folder and extracts
% information about the date and time of the recordings, which is used
% to plot the timeline of the experiment.

% Extract date and time information from the folder name:
% 1 - find the string that has the date and time
date_ind            = strfind(folder_name,'_20');
date_start          = date_ind(1) + 1;
date_str            = folder_name(date_start:date_start+18);
time_str            = date_str(12:19);

% 2 - break down the date and time into individual values
this_year           = str2num(date_str(1:4));
this_month          = str2num(date_str(6:7));
this_day          	= str2num(date_str(9:10));

this_hour           = str2num(time_str(1:2));
this_min            = str2num(time_str(4:5));
this_sec            = str2num(time_str(7:8));

% Make a time variable that can be compared to the start time of the
% first protocol that was run
protocol_time       = datetime(this_year,this_month,this_day,this_hour,this_min,this_sec);
