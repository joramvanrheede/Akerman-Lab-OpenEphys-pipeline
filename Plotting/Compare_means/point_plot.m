function point_plot(in_data, in_groups,x_vals)
% point_plot_h    = point_plot(in_data, in_groups, x_vals)
% or:
% point_plot_h    = point_plot(in_data,x_vals)
% 
% Make point plot as alternative to bar graph
% Adds individual data points as points with alpha level, plus mean and 
% standard error lines
% 
% if in_data is a vector, the function requires 'in_groups' to specify for 
% each point in in_data what group it belongs to.
% 
% if in_data is a matrix, each column in the matrix is treated as a separate
% group.
% 
% function will plot the points at 1:n_groups unless x_vals are specified;
% if X-vals are specified these will determine the X position for each group.
% 


marker_alpha    = [0.3];

if ismatrix(in_data)
    if exist('in_groups','var')
        x_vals      = in_groups;
    end
    in_groups       = repmat(1:size(in_data,2),size(in_data,1),1);
    in_data         = in_data(:);
    in_groups       = in_groups(:);
    group_nrs       = unique(in_groups);
    
    if nargin < 2
        x_vals  = 1:length(group_nrs);
    end
else
    if nargin < 2
        in_groups = ones(size(in_data));
    end
    group_nrs       = unique(in_groups);
    if nargin < 3
        x_vals   = 1:length(group_nrs);
    end
end

% How wide are mean and standard error line?
mean_width      = 0.5 * median(diff(x_vals));
serr_width      = 0.25 * median(diff(x_vals));

for a = 1:length(group_nrs)
    data_in_group   = in_data(in_groups == group_nrs(a)); % Select data for this group
    
    this_x          = x_vals(a);
    
    scatter(repmat(this_x,1,length(data_in_group)),data_in_group,50,'o','filled','MarkerFaceColor',[1 0 0],'MarkerEdgeAlpha',marker_alpha,'MarkerFaceAlpha', marker_alpha);
    hold on
    
    line([this_x-mean_width this_x+mean_width],[mean(data_in_group) mean(data_in_group)],'LineWidth',3,'Color',[0 0 0])
    
    line([this_x-serr_width this_x+serr_width],[mean(data_in_group)+serr(data_in_group) mean(data_in_group)+serr(data_in_group)],'LineWidth',3,'Color',[0 0 0])
    line([this_x-serr_width this_x+serr_width],[mean(data_in_group)-serr(data_in_group) mean(data_in_group)-serr(data_in_group)],'LineWidth',3,'Color',[0 0 0])
    
    line([this_x this_x],[mean(data_in_group)-serr(data_in_group) mean(data_in_group)+serr(data_in_group)],'LineWidth',3,'Color',[0 0 0])
    
end

