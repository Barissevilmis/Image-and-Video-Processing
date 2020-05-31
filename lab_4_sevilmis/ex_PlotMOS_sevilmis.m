function j = ex_PlotMOS_sevilmis(mos, ci, bitrates, cnts, cdc, contents, codecs,tmp)

    %MOS VALUES VS BITRATE - PLOT - SPECIFIC PER CONTENT
    figure('Name','MOS plots','Color',[1 1 1]);
    title(['MOS vs BPP for ', contents{tmp}]);
    xlabel('BPP');
    ylabel('MOS');
    hold on; grid on;
    
    %colors = ['r', 'b', 'g', 'm'];
    shapes = ["--", "-.+", ":*", "-h"];
    for i = 1 : 4
        data = intersect(cnts, cdc(i,:));
        errorbar(bitrates(data), mos(data), ci(data), shapes(i),'LineWidth',1.5);
    end
    legend(codecs, 'Location', 'SouthEast');
    %saveas(gcf, ['res/mos',num2str(tmp),'.png']);
    
end