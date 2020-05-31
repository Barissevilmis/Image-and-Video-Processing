function [p, s, rmse] = ex_SpearmanPearsonRMSE_sevilmis(mos, data, content)

    % COMPUTE PEARSON, SPEARMAN COEEFICIENTS AND RMSE
    p=corr(data(content),mos(content),'Type','Pearson');
    s=corr(data(content),mos(content),'Type','Spearman');
    rmse=sqrt(sum((data(content)-mos(content)).^2)/length(mos(content)));
    
end