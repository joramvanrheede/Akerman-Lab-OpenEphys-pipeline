function corr_data  = unit_correlations(merged_data,corr_type,mua_channels,psth_edges,resp_win)


if nargin < 2 || isempty(corr_type)
    corr_type       = 'spont';
end

if nargin < 3 || isempty(mua_channels)
    mua_channels       = ':'; 
end

if nargin < 4 || isempty(psth_edges)
	psth_edges      = [0:0.01:3];
end

if nargin < 5
    if strcmpi(corr_type,'spont')
        resp_win    = [0 0.5];
    else
        resp_win    = [];
    end
end

sta_psth_bins   = [-0.1:0.001:0.1]; % spike triggered average PSTH bins
smooth_gauss    = normpdf(-3:0.25:3);

unit_depths     = merged_data.unit_depths;

n_units         = length(unit_depths);
n_channels      = size(merged_data.conditions(1).multiunit_spikes,1);
n_conditions    = length(merged_data.conditions);

spont_corr_win  = [0 0.5];

% spike-triggered LFP?
% spike-triggered population rate
% Population coupling index from Okun paper
% Population percentile index
% correlate across trials, then avg over time bins?

cond_unit_psths         = NaN(n_units, n_conditions, length(psth_edges)-1);
triggered_percentile    = NaN(n_units, n_conditions);
triggered_psths         = NaN(n_units, n_conditions, length(sta_psth_bins)-1);
triggered_norm_psths    = NaN(n_units, n_conditions, length(sta_psth_bins)-1);
for a = 1:n_conditions
    
    this_cond   = merged_data.conditions(a);
    unit_spikes = this_cond.spikes;
    mua_spikes  = this_cond.multiunit_spikes(mua_channels,:,:);
    
    switch lower(corr_type)
        case 'whisk'
            unit_spikes     = unit_spikes - this_cond.whisk_onset;
            mua_spikes      = mua_spikes - this_cond.whisk_onset;
        case 'opto'
            unit_spikes  	= unit_spikes - this_cond.LED_onset;
            mua_spikes   	= mua_spikes - this_cond.LED_onset;
        case 'spont'
            % do nothing; case included to avoid error in 'otherwise' statement below
        otherwise
            error('Unrecognised correlation option')
    end
	
    for i = 1:n_units
        these_spikes  = squeeze(unit_spikes(i,a,:));
        
        [~, cond_unit_psths(i,a,:)]  = psth(these_spikes(i,:,:),psth_edges,false);
        
    end
    
    % USE
    % spike_triggered_mua function
    
    [triggered_mua_psths, norm_triggered_psths, n_valid_spikes]     = spike_triggered_mua(unit_spikes, mua_spikes, resp_win,sta_psth_bins);
    triggered_nbins             = size(triggered_mua_psths,2);
    half_nbins                  = floor(triggered_nbins / 2);
    triggered_percentile(:,a)   = sum(norm_triggered_psths(:,1:half_nbins),2);
    
    triggered_psths(:,a,:)      = triggered_mua_psths;
    triggered_norm_psths(:,a,:) = norm_triggered_psths;
%     triggered_percentile    = sum(norm_triggered_psths)
%     
end

corr_data.cond_unit_psths       = cond_unit_psths;
corr_data.triggered_psths       = triggered_psths;
corr_data.triggered_norm_psths  = triggered_norm_psths;
corr_data.triggered_percentile  = triggered_percentile;
corr_data.smooth_STPR           = smoothdata(triggered_psths,3,'gaussian',12);
corr_data.smooth_norm_STPR     	= smoothdata(triggered_norm_psths,3,'gaussian',12);

