clear all;
close all;
clc;

sorted_folder       = 'F:\Synced';
group_folders     	= {'Pom';'Cortical'}; %

for k = 1 : numel(group_folders);
   load([sorted_folder '\' group_folders{k} '_exps.mat']);
    Group_results(k).exps = exps;
    Group_results(k).L23_Opto_all = L23_Opto_all;
    Group_results(k).L23_Whisk_all = L23_Whisk_all;
    Group_results(k).L23_Spont_all = L23_Spont_all;
    Group_results(k).L5_Opto_all = L5_Opto_all;
    Group_results(k).L5_Spont_all = L5_Spont_all;
    Group_results(k).L5_Whisk_all = L5_Whisk_all;
    Group_results(k).Group = group_folders{k};
    
    
    
    %[rho(k),pval(k)] = corr(Group_results(k).L23_Opto_all',Group_results(k).L5_Opto_all','Type','Spearman','Rows','Pairwise');
    %[rho(k),pval(k)] = corr(Group_results(k).L23_Opto_all',Group_results(k).L5_Opto_pos','Type','Spearman','Rows','Pairwise');
    
    
    clear exps L23_Opto_all L23_Whisk_all L23_Spont_all L5_Opto_all L5_Whisk_all L5_Spont_all;
end;
clear k

for b = 1 : numel(Group_results);
rhotot = []
    for a = 1 : numel(Group_results(b).exps);
        Exp_of_interest = Group_results(b).exps(a);
        L5_Opto = Exp_of_interest.L5_whisk_rates' - Exp_of_interest.L5_spont_rates';
        L23_Opto = Exp_of_interest.L23_whisk_rates' - Exp_of_interest.L23_spont_rates';
 
        toremove = [];
        for k = 1:size(L5_Opto,2);
            if sum(L5_Opto(:,k)) == 0;
                toremove =[toremove k];
            end;
        end;
       L5_Opto(:,toremove) = [];
  
       toremove = [];
       for k = 1:size(L23_Opto,2)
           if sum(L23_Opto(:,k)) == 0;
                toremove =[toremove k];  
           end;
       end;
       L23_Opto(:,toremove) = [];
 
      clear rho pval;
       if ~isempty(L23_Opto) & ~isempty(L5_Opto)
           [rho,pval] = corr(L23_Opto,L5_Opto,'Type','Spearman','Rows','pairwise');
         %{  
           ha = figure('Name',[group_folders{b} ' Opto Correlation'],'NumberTitle','off');
    set(gcf,'Units','normalized','Position',[.2 .1 .4 .8])% left bottom width height

                    ax1 = subplot(2,2,1);
                    im = imagesc(rho);
                    title('Correlation between L23 and L5 units');
                    xlabel('Layer 5 units');
                    ylabel('Layer 2/3 units'); 
                    colormap(ax1,hot)
                    c = colorbar;
                    c.Label.String = 'Correlation Coefficient (\rho)';
                    rhohist = reshape(rho,numel(rho),1);
                     subplot(2,2,2)
                     histogram(rhohist);
                     xlabel('Correlation Coefficient (\rho)');
                     ylabel('count');
 
                     bin_edges = 0:0.01:1;
                     
                     ax2 = subplot(2,2,3);
                     title('Correlation between L23 and L5 units');
                     clims = [0.049 0.051];
                     im = imagesc(pval,clims);
                     xlabel('Layer 5 units');
                     ylabel('Layer 2/3 units'); 
                      colormap(ax2,gray);
                     c = colorbar('Ticks',[0.01 0.05 0.1]);
                      c.Label.String = 'Significance of Correlation (p)';
                     phist = reshape(pval,numel(pval),1);
                     subplot(2,2,4);
                     histogram(phist,bin_edges);
                     xlabel('Significance of Correlation (p)');
                     ylabel('count');
 
 
 
           %}
           
           rho(pval >1) = NaN;
           rho = rho(~isnan(rho));
           rho = reshape(rho,numel(rho),1);
           rhotot = [rhotot;rho];
         end;
end;
correlated(b).rhotot = rhotot;
end;

POM_Count = correlated(1).rhotot;
Cort_Count = correlated(2).rhotot;

POM_Count(POM_Count>0) = 0;
POM_Count(POM_Count<0) = 1;

Cort_Count(Cort_Count>0) =0;
Cort_Count(Cort_Count<0) =1;


ha = figure('Name','Correlations','NumberTitle','off');
    set(gcf,'Units','normalized','Position',[.2 .1 .4 .8])% left bottom width height
    subplot(2,1,1);
    bin_edges = [-0.8:0.2:0.8];
    h_1 = histogram(correlated(1).rhotot,bin_edges);
    
    title(group_folders(1));
    subplot(2,1,2);
    histogram(correlated(2).rhotot,bin_edges);
    [p,h] = ranksum(correlated(1).rhotot,correlated(2).rhotot);    
    title(group_folders(2));
    
hb = figure('Name','Correlations','NumberTitle','off');
    set(gcf,'Units','normalized','Position',[.2 .1 .4 .8])% left bottom width height
    
    subplot(2,1,1)
    histogram(L23_Opto);
    title('Evoked Activity in Layer 2/3 unit');
    xlabel('Delta Spike (Hz)');
    ylabel('Trials');
    subplot(2,1,2)
    histogram(L5_Opto(:,9));
    title('Evoked Activity in Layer 5 unit');
    xlabel('Delta Spike (Hz)');
    ylabel('Trials');
    

 
 %{   
 
 
 

ha = figure('Name','Whisker rates','NumberTitle','off');
    set(gcf,'Units','normalized','Position',[.2 .1 .4 .8])% left bottom width height
    subplot(3,1,1);
    title('Relationship between layers');
    hold on;
    scatter(Group_results(1).L5_Whisk_all,Group_results(1).L23_Whisk_all,'r*');
    scatter(Group_results(2).L5_Whisk_all,Group_results(2).L23_Whisk_all,'b');
    legend(group_folders);
     ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    xlabel('Delta firing rate in L5(Hz)');
    ylabel('Delta firing rate in L2/3(Hz)');
    subplot(3,1,3);
    hold on;
    histogram(Group_results(1).L5_Whisk_all);
    histogram(Group_results(2).L5_Whisk_all);
    legend(group_folders);
    title('Layer 5 rates'); 
    subplot(3,1,2);
    hold on;
    histogram(Group_results(1).L23_Whisk_all);
    histogram(Group_results(2).L23_Whisk_all);
    legend(group_folders);
    title('Layer 23 rates');
   
    hb = figure('Name','Spontaneous rates','NumberTitle','off');
    set(gcf,'Units','normalized','Position',[.2 .1 .4 .8])% left bottom width height
    subplot(3,1,1);
    hold on;
    title('Relationship between layers');
    scatter(Group_results(1).L5_Spont_all,Group_results(1).L23_Spont_all,'r*');
    scatter(Group_results(2).L5_Spont_all,Group_results(2).L23_Spont_all,'b');
    legend(group_folders);
    xlabel('Mean firing rate in L5(Hz)');
    ylabel('Mean firing rate in L2/3(Hz)');
    subplot(3,1,3);
    hold on;
    histogram(Group_results(1).L5_Spont_all);
    histogram(Group_results(2).L5_Spont_all);
    legend(group_folders);
    
    title('Layer 5 rates'); 
    subplot(3,1,2);
    hold on;
    histogram(Group_results(1).L23_Spont_all);
    histogram(Group_results(2).L23_Spont_all);
    legend(group_folders);
    title('Layer 23 rates');
    
    hc = figure('Name','Opto evoked rates','NumberTitle','off');
    set(gcf,'Units','normalized','Position',[.2 .1 .4 .8])% left bottom width height
    subplot(2,1,1);
    hold on;
    title('Relationship between layers');
    scatter(Group_results(1).L5_Opto_all,Group_results(1).L23_Opto_all,'r*');
    scatter(Group_results(2).L5_Opto_all,Group_results(2).L23_Opto_all,'b');
    legend(group_folders);
     ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    xlabel('Delta firing rate in L5(Hz)');
    ylabel('Delta firing rate in L2/3(Hz)');
    
    subplot(2,1,2);
    hold on;
    title('Relationship between layers');
    scatter(Group_results(1).L5_Opto_pos,Group_results(1).L23_Opto_all,'r*');
    scatter(Group_results(2).L5_Opto_pos,Group_results(2).L23_Opto_all,'b');
    legend(group_folders);
     ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    xlabel('Delta firing rate in L5(Hz)');
    ylabel('Delta firing rate in L2/3(Hz)');
    
    
     hd = figure('Name','Opto evoked rates','NumberTitle','off');
    set(gcf,'Units','normalized','Position',[.2 .1 .4 .8])% left bottom width height
    subplot(2,1,2);
    hold on;
    histogram(Group_results(1).L5_Opto_all,'Normalization','probability');
    histogram(Group_results(2).L5_Opto_all,'Normalization','probability');
    legend(group_folders);
    title('Layer 5 rates'); 
    subplot(2,1,1);
    hold on;
    histogram(Group_results(1).L23_Opto_all,'Normalization','probability');
    histogram(Group_results(2).L23_Opto_all,'Normalization','probability');
    legend(group_folders);
    title('Layer 23 rates');
   %}
  