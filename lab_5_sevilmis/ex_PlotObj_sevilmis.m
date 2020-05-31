function ex_PlotObj_sevilmis(score, bitrates, cnts, cdc, contents, codecs,n, tmp, name)
    
    % Specific content: Plot for each codec seperately. FR vs BPP
    figure('Name',['Obj(',name,') Metrics vs Bitrates'],'Color',[1 1 1]);
    title([name ,' vs BPP for ', contents{n}]);
    xlabel('BPP');
    ylabel(name);
    hold on; grid on;
    
    for i = 1 : 4
        data = intersect(cnts, cdc(i,:));
        plot(bitrates(data), score(data));
    end
    legend(codecs, 'Location', 'SouthEast');
    %saveas(gcf, ['results/obj',num2str(tmp),'.png']);
    
end