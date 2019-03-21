

% run process_channels_generic first
probe_traces = permute(squeeze(sdata(1).expt(1).LFP_traces(66,:,1:2000)),[2 1]);

diff_table  = [];
figure
for a = 1:size(probe_traces,2)
    
	thistrace           = probe_traces(:,a);
    
    plot(thistrace,'Color',[0.2 0.2 0.2]+a*0.025,'LineWidth',2)
    hold on
    
    testmat             = repmat(thistrace,[1,size(probe_traces,2)]);
    
    diff_mat            = abs(probe_traces - testmat);
    diff_sum            = sum(diff_mat);
    diff_sum(a)         = NaN;
    
    diff_table          = [diff_table; diff_sum];
    [mindiff, minind]   = min(diff_sum);
    
    disp(['For channel ' num2str(a) ' the most similar channel is ' num2str(minind)])
end

%% Note - can also do something like this with 'corrcoef'


figure
imagesc(diff_table)
figure
imagesc(corrcoef(probe_traces))

colormap gray
ylabel('Channel number')
xlabel('Channel number')
title('Differences between channels')



