function j = ex_PlotNR_sevilmis(score, bitrates, cnts, cdc, contents, codecs, n,tmp,name)

    % Specific content: Plot for each codec seperately. NR vs BPP
    figure('Name',['NR(',name,') Metrics vs Bitrates'],'Color',[1 1 1]);
    title([name,' vs BPP for ', contents{n}]);
    xlabel('BPP');
    ylabel(['NR: ',name]);
    hold on; grid on;
    
    for i = 1 : 4
        data = intersect(cnts, cdc(i,:));
        plot(bitrates(data), score(data));
    end
    legend(codecs, 'Location', 'SouthEast');
    %saveas(gcf, ['res/nr',num2str(tmp),'.png']);

end