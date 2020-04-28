
clear opto_effect opto_effect_pom L2_3 unit med neg;
hold on;
%close all;
clear opto_effect opto_effect_pom L2_3 unit med neg spont_corr_unit ;
clear spont_corr_local spont_corr_l5 opto_corr_l5 opto_corr_local opto_corr_unit rho rh_5 o_5 o diff diff_5; 

all_times(1:2,1:7,1:200) = NaN;
sig_times(1:2,1:7,1:200) = NaN;
corr_bins = (3:18)     
b =1;

diff(1:2,1:200) = NaN;
    diff_5(1:2,1:200) = NaN;
   
    
rho_local(1:2,1:200) = NaN; 
pval_local(1:2,1:200) = NaN; 
rho_l5(1:2,1:200) = NaN; 

for b = 1:2
    figure('Name',['group' num2str(b)]);
    clear spont_corr_unit spont_corr local opto_corr_unit opto_corr_local
    count = 0
for k = 1:numel(Analysed_group(b).all_units)
unit = Analysed_group(b).all_units(k);
unit_delta = [];
l5_delta = [];
spont_corr_unit = [];
spont_corr_local = [];
spont_corr_l5 = [];
opto_corr_unit = [];
opto_corr_local = [];
opto_corr_l5 = [];


if unit.unit_depth < 225; %& ~unit.whisk_responsive;
    for j = [1 7]  
        spont_corr_unit = [spont_corr_unit;sum(unit.conditions(j).spont_behaviour.psth_100(:,corr_bins),2)];
        spont_corr_local = [spont_corr_local sum(unit.conditions(j).spont_behaviour.multi_local_psth_100(corr_bins,:),1)];
        spont_corr_l5 =  [spont_corr_l5 sum(unit.conditions(j).spont_behaviour.multi_l5_psth_100(corr_bins,:),1)];
        opto_corr_unit = [opto_corr_unit;sum(unit.conditions(j).opto_behaviour.psth_100(:,corr_bins),2)];
        opto_corr_local = [opto_corr_local sum(unit.conditions(j).opto_behaviour.multi_local_psth_100(corr_bins,:),1)];
        opto_corr_l5 =  [opto_corr_l5 sum(unit.conditions(j).opto_behaviour.multi_l5_psth_100(corr_bins,:),1)];
        end;
        delta_unit_spike = opto_corr_unit-spont_corr_unit;
        delta_local_spike = opto_corr_local-spont_corr_local;
        delta_l5_spike = opto_corr_l5-spont_corr_l5;
        hold on;
     [r,p] = corr(spont_corr_unit,spont_corr_local','type','spearman');
     [r_2,p_2] = corr(opto_corr_unit,opto_corr_local','type','spearman');
     
rho_local(b,k) = r_2-r;
pval_local(b,k) = p;

subplot(2,1,1);
hold on;
        if pval_local(b,k) <0.05
            hold on;    
            plot(1,rho_local(b,k),'g*');
        else 
        hold on;
        plot(1,rho_local(b,k),'ro');
        end;
 title('Unit Activity correlation with Local multiunit activity');
 ylabel('delta spearman rho');
        
     [r,p] = corr(spont_corr_unit,spont_corr_l5','type','spearman');
     [r_2,p_2] = corr(opto_corr_unit,opto_corr_l5','type','spearman');
     
rho_l5(b,k) = r_2-r;
pval_l5(b,k) = p;
subplot(2,1,2);
hold on;
        if pval_l5(b,k) <0.05
            hold on;    
            plot(1,rho_l5(b,k),'g*');
        else 
        hold on;
        plot(1,rho_l5(b,k),'ro');
        end;
title('Unit Activity correlation with L5 multiunit activity');
ylabel('delta spearman rho');
  %}
end;
  
end;
end;
 
[h,p] = ranksum(rho_local(1,:),rho_local(2,:));
[h_5,p_5] = ranksum(rho_l5(1,:),rho_l5(2,:));


