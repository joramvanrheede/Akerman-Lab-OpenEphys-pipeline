clear h i j k r ref spikes t target test this_cond corr_window bin_edges stim_spikes_ref stim_spikes_target jpsth j_sum jpsth_windows;
close all
this_cond = ephys_data.conditions(1);
stim_time = ephys_data.conditions(1).whisk_onset;
spikes = this_cond.spikes;
corr_window = 25;
bin_size = 1;
psth_bins       = [0:1:100];

stim_spikes = spikes;
stim_spikes = stim_spikes;

target = squeeze(stim_spikes(6,:,:));
ref = squeeze(stim_spikes(8,:,:));    

target_shift = target(2:end,:);
target_shift = [target_shift;target(1,:)];

[N,Edges] = Cross_Correlogram(ref,target,corr_window,bin_size);
[Period_shift,Edges_2] = Cross_Correlogram(ref,target_shift,corr_window,bin_size);
Corrected = N-Period_shift;

clear x y;


Lamda = Period_shift;
%variance = (std(Period_shift))^2;
for x = 1 : numel(N);
y(x) = Poisson_prob(N(x),Lamda(x));
end;
figure()
scatter(N,y);


figure();
subplot(3,1,1);
bar(Edges(1:end-1),N);
title('Cross Correlogram')

subplot(3,1,2);
bar(Edges_2(1:end-1),Period_shift);
title('Period shifted correlogram');

subplot(3,1,3);
bar(Edges(1:end-1),Corrected);
title('Corrected Correlogram');


%{
for z = 1: size(stim_spikes,2)
stim_spikes_ref = ref(z,:);
stim_spikes_target = target(z,:);

q = jpsth_window(1) < stim_spikes_ref & stim_spikes_ref <jpsth_window(2);
stim_spikes_ref = stim_spikes_ref(q);

q = jpsth_window(1) < stim_spikes_target & stim_spikes_target <jpsth_window(2);
stim_spikes_target = stim_spikes_target(q);

q = histcounts(stim_spikes_ref,bin_edges);
r = histcounts(stim_spikes_target,bin_edges);
disp(['q :' num2str(max(q)) 'r:' num2str(max(r))]);
for k = 1:numel(q)
for j = 1:numel(r)
    a(k,j,z) = q(k) & r(j);
end;
end;
clear q r k j;
end;
jpsth = sum(a,3);


bin_edges = 1:bin_size/0.1:size(jpsth,1);
for j = 1:numel(bin_edges)-1
for k = 1:numel(bin_edges)-1
    j_sum(k,j) = sum(sum(jpsth(bin_edges(k):bin_edges(k+1),bin_edges(j):bin_edges(j+1))));
end;
end;
clear a;

colormap('hot');   % set colormap
im = imagesc(j_sum);        % draw image and scale colormap to values range
title('LFP Whisk Alone');
im.YData = [size(j_sum,1),1]
xlabel('unit x');
ylabel('unit y');
c = colorbar;          % show color scale
c.Label.String = 'Spike Count';

%}

