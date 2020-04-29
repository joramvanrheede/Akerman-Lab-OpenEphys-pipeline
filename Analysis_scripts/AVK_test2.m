close all;

for a = 1 %: numel(Analysed_group)
    figure();
for k = 50 : numel(Analysed_group(a).all_units);
    unit = Analysed_group(a).all_units(k);
    disp(k)
    hold on;
    quick_plot_unit(unit,'grid');
    
end;
end;    
    
disp('done');
