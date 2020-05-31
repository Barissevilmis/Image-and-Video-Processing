function [J, r1, r2] = ex_Outlier_sevilmis(data)

    % FIND R1 AND R2 USING PEARSON
    r1 = zeros(size(data,2),1);
    r2 = zeros(size(data,2),1);
    for i = 1 : size(data,2)
        r1(i) = corr(ex_MOS_sevilmis(data), data(:,i),'Type','Pearson');
        r2(i) = corr(ex_MOS_sevilmis(data), ex_MOS_sevilmis(data(:,i)),'Type','Pearson');
    end
    
    % REMOVE OUTLIERS
    outliers = (r1 < 0.75) & (r2 < 0.8);
    J = data(:,~outliers);
    
end
