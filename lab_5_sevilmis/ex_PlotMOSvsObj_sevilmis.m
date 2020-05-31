function ex_PlotMOSvsObj_sevilmis(mos, ci, scores, cnts, contents,tmp,name)

    % SUBJECTIVE VS OBJECTIVE SCORES - ALL CONTENTS
    figure('Name','MOS plots','Color',[1 1 1]);
    title(['MOS vs ',name,' for all contents']);
    xlabel(name);
    ylabel('MOS');
    hold on; grid on;
    
    %colors = ['r', 'b', 'g', 'm'];
    shapes = ["--", "-.+", ":*", "-h","-.p"];
    for i = 1 : 4
        errorbar(scores(cnts(i,:)), mos(cnts(i,:)), ci(cnts(i,:)), shapes(i),'LineWidth',1.5);
    end
    %plot(scores(data), mos(data));
    legend(contents, 'Location', 'SouthEast');
    %saveas(gcf, ['results/mosvsobj',num2str(tmp),'.png']);
    
end