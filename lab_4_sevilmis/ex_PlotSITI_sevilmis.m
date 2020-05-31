function j = ex_PlotSITI_sevilmis(si, ti,tmp)
    
    %PLOT SI and TI - BAR PLOTS - FOR ALL CONTENTS
    figure('Name','SI & TI Metrics','Color',[1 1 1]);
    title('SI and TI');
    xlabel('Contents');

    names = {'Campfire','Runners','Sintel2','TrafficFlow','TreeShades'};
    set(gca,'xtick',[1:5],'xticklabel',names);
    ylabel('SI & TI');
    hold on; grid on;
    x = categorical({'Campfire', 'Runners', 'Sintel2', 'TrafficFlow', 'TreeShades'});

    data = zeros(size(ti,1),2);
    for i = 1 : size(ti,1)
        data(i,1) = ti(i);
        data(i,2) = si(i);
    end
    b = bar(data);
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels1 = string(b(1).YData);
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
    
    xtips2 = b(2).XEndPoints;
    ytips2 = b(2).YEndPoints;
    labels2 = string(b(2).YData);
    text(xtips2,ytips2,labels2,'HorizontalAlignment','center', 'VerticalAlignment','bottom')   
    
    legend({'TI', 'SI'}, 'Location', 'NorthEast');
    saveas(gcf, ['res/siti',num2str(tmp),'.png']);
    
end