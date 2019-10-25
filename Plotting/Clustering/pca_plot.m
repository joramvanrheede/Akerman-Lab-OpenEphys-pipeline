function [pca_res, plot_h] = pca_plot(indata,n_dims, group_var)
% function [pca_res, plot_h] = pca_plot(indata,n_dims, group_var)
% 
% Plot principal components coeffs of indata against each other as a 2d or 3d 
% scatter plot
% 

if nargin == 1
    n_dims = 2;
elseif n_dims < 2 || n_dims > 3
    error('Only 2 or 3 dimensions are supported')
end

if nargin < 3
    group_var   = ones(size(indata,2),1);
end

% Create colour variable, black by default
colour_var       = zeros(length(group_var),3);
size(colour_var)
% Hard-code some sensible colours
% 1 = black
colour_var(group_var == 1,1)        = 1;    % red
colour_var(group_var == 2,3)        = 1;    % blue
colour_var(group_var == 3,2)        = .7;   % darkish green
colour_var(group_var == 4,[1 3])    = 1;    % Violet?
colour_var(group_var == 5,[2 3])  	= .5;   % Green-blue
colour_var(group_var > 5,:)         = .5;   % The rest is grey

colour_var

% Do PCA
[coeff, score, latent, tsquared, explained] = pca(indata);

if n_dims == 2
    
    plot_h = scatter(coeff(:,1),coeff(:,2),50,colour_var,'filled');
    
    xlim([-1 1])
    ylim([-1 1])
elseif n_dims == 3
    
    plot_h = scatter3(coeff(:,1),coeff(:,2),coeff(:,3),50,colour_var,'filled');
    
    xlim([-1 1])
    ylim([-1 1])
    zlim([-1 1])
end

fixplot

% output
pca_res.coeff       = coeff;
pca_res.score       = score;
pca_res.latent      = latent;
pca_res.tsquared   	= tsquared;
pca_res.explained   = explained;



