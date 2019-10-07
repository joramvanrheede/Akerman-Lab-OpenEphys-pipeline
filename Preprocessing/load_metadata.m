function [metadata, headers]    = load_metadata(metadata_file)
% [METADATA, HEADERS]    = load_metadata(METADATA_FILE)
% 
% Load metadata, for further processing / reading by READ_METADATA function.
% 
% METADATA_FILE: Full path to metadata file.
% 
% METADATA: The data contained in the metadata spreadsheet as a large cell.
%
% HEADERS: The row with the headers of the metadata spreadsheet.
% 

% Read excel file containing metadata
[~, ~, metadata]    = xlsread(metadata_file,1); % use XLSread function to read the excel spreadsheet
headers             = metadata(1,:); % top line contains the headers for the columns
metadata            = metadata(2:end,:); % the rest is data
