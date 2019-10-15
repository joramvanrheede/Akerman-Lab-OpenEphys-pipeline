function LFP_colors = LFP_colormap
% function LFP_colors = LFP_colormap
% 
% Colormap for LFP visualisation; ranges from blue (low) through white (middle)
% to red (high).
% 

fade_streak             = [0:63] / 63;
fade_columns            = [fade_streak' fade_streak'];

LFP_colors              = ones(128,3);

LFP_colors(1:64, 1:2)   = fade_columns;

LFP_colors(65:128, 2:3) = flipud(fade_columns);
