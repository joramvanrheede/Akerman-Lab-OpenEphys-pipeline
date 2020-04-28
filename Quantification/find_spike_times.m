function[first_spike,second_spike] = find_spike_times(spikes,window);
 
for k = 1: size(spikes,1);
    a = find(spikes(k,:) > window(1),1,'first');
    if~isempty(a)
    first_spike(k) = spikes(k,a);
    second_spike(k) = spikes(k,a+1);
    else
    first_spike(k) = NaN;  
    second_spike(k) = NaN;
    end;
end;

    first_spike(first_spike > window(2)) = NaN;
    second_spike = second_spike(~isnan(first_spike));
    second_spike(second_spike > 0.5) = NaN;


end