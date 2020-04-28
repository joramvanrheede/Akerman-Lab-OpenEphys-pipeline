clear all; close all; clc;
load_folder = 'F:\Synced'
group_folders     	= {'POM';'Cortical'}; 

for k = 1: numel(group_folders);
disp(['Loading' group_folders{k} 'unit_behaviour.mat']) 
load([load_folder '\' group_folders{k} 'unit_behaviour.mat']);
Group(k).behaviour = Unit_behave;
disp('loaded');
clear Unit_behave;
end

ha = figure('Name','Unit Behaviour','NumberTitle','off');
  set(gcf,'Units','normalized','Position',[.05 .1 0.9 .8]) %left bottom width height   
for j = 1: numel(Group);
    Group(j).All_whisk = [];
 for k = 1:numel(Group(j).behaviour)
    a = Group(j).behaviour(k).all_unit_whisk_rate';
    b = Group(j).behaviour(k).all_unit_rate;
    Group(j).behaviour(k).whisk_spont_ratio = a;%./(a+b);
 
end;

 for k = 1:numel(Group(j).behaviour)
    Group(j).behaviour(k).whisk_spont_ratio = Group(j).behaviour(k).whisk_spont_ratio - Group(j).behaviour(end).whisk_spont_ratio;
    subplot(2,2,(2*j)-2+1)
    hold on;
    plot(k,Group(j).behaviour(k).whisk_spont_ratio,'bo');
    title (['Delta Whisk Responsiveness ratio over time :' group_folders{j}]);
    ylabel('whisk spike rate/(whisk+spont spike rate)');
    xlabel('Time point');
    refline(0,0);

    subplot(2,2,(2*j)-2+2)
    hold on;
    title (['Delta Whisk Responsiveness against depth :' group_folders{j}]);
    scatter(Group(j).behaviour(k).all_unit_depth,Group(j).behaviour(k).whisk_spont_ratio);
    legend({'t1','t2','t3','t4','t5','t6','control'});
    xlabel('Unit "depth" (800-kilosort_depth)');
    ylabel('whisk spike rate/(whisk+spont spike rate)');
   

end;

for k = 1:numel(Group(j).behaviour)
[h,p] =  ranksum(Group(j).behaviour(k).whisk_spont_ratio,Group(j).behaviour(end).whisk_spont_ratio);
Group(j).behaviour(k).h = h;
Group(j).behaviour(k).p = p;
Group(j).behaviour(k).mean_behaviour = nanmean(Group(j).behaviour(k).whisk_spont_ratio);
Group(j).All_whisk = [Group(j).All_whisk Group(j).behaviour(k).whisk_spont_ratio]
end

effect_test = Group(j).All_whisk(:,3:6);
[m,index] = max(abs(effect_test),[],2);
disp(['Group : ' num2str(j)]);
disp([num2str(index)]);

for q = 1 : numel(index)
Group(j).Max_effect(q) = effect_test(q,index(q));
end;
end;

ha = figure('Name','Unit Behaviour','NumberTitle','off');
  set(gcf,'Units','normalized','Position',[.05 .1 0.9 .8]) %left bottom width height
  hold on;
scatter(Group(1).behaviour(1).all_unit_depth,Group(1).Max_effect,'bo');
scatter(Group(2).behaviour(1).all_unit_depth,Group(2).Max_effect,'ro');
xlabel('Unit Depth (800-Kilosort uM)');
ylabel('Unit whisk response/spont response');
legend({'POm group';'Cortical group'});



