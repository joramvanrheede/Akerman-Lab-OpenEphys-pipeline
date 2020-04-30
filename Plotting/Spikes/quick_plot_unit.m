function quick_plot_unit(unit,stim_type);

 if (stim_type == 'opto') | (stim_type == 'both')    
        %figure('Name','Optoeffect'); %left bottom width height
        %set(gcf,'Units','normalized','Position',[0.3 .3 .4 .4]);
      
        opto_effect = sum((unit.conditions(end).opto_behaviour.psth_100-unit.conditions(end).spont_behaviour.psth_100),1);
        bin_no = size(opto_effect,2); 
        x = 0.005:0.005:bin_no*0.005;
        
        %opto_linear = sum(unit.conditions(line_test).whisk_behaviour.psth_100)./sum(unit.conditions(line_test).opto_behaviour.psth_100+unit.conditions(end).whisk_behaviour.psth_100);
        
        %bar(x,opto_effect);
        ylabel('spike_count');
        xlabel('time bin (s)');
               
      end;
      
      if (stim_type == 'grid') | (stim_type == 'both')
        figure('Name',['Unit depth :' num2str(unit.unit_depth)]); %left bottom width height
        set(gcf,'Units','normalized','Position',[0.05 .1 .9 .9]);
             bin_no = 20; 
             x = 0.005:0.005:bin_no*0.005;
              no_conds = numel(unit.conditions);
     
            for k = 1: no_conds;
                raster_control = unit.conditions(k).spont_behaviour.raster_100;
                raster_opto = unit.conditions(k).opto_behaviour.raster_100;
                raster_data = unit.conditions(k).whisk_behaviour.raster_100;
                psth_plot = sum(unit.conditions(k).whisk_behaviour.psth_100,1);
                psth_control = unit.conditions(end).whisk_behaviour.raster_100;
                
                subplot(4,no_conds,k);
                raster_plot(raster_control,1,0.2);
                title('spontaneous');
                ylabel('Trial no');
                xlabel('Time(s) post stimulus');

                subplot(4,no_conds,no_conds+k);
                raster_plot(raster_opto,1,0.2);
                title('opto');
                ylabel('Trial no');
                xlabel('Time(s) post stimulus');

                
                subplot(4,no_conds,2*no_conds+k);
                raster_plot(raster_data,1,0.2);
                ylabel('Trial no');
                xlabel('Time(s) post stimulus');

                
                subplot(4,no_conds,3*no_conds+k);
                hold on;
                b = bar(x,psth_plot,'hist','b');
                h = histogram(psth_control,x);
                h.FaceAlpha = 0.2;
                h.FaceColor = 'r';
                ylabel('spike count');
                xlabel('bin');
                
            end;
            
        end;

end