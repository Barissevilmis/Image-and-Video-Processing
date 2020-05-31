function out = ex_Outlier_sevilmis(data)

    Pi = zeros(size(data, 1), size(data, 2));
    Qi = zeros(size(data, 1), size(data, 2));
    
    % COMPUTE KURTOSIS AND APPLY THRESHOLDS
    for i = 1 : size(data, 1)
        k = kurtosis(data(i, :));
        mu = mean(data(i, :));
        sigma = std(data(i, :));
        for j = 1 : size(data ,2)
            if k >= 2 && k <= 4        
                Pi(i, j) = data(i, j) < mu - (2*sigma);
                Qi(i, j) = data(i, j) > mu + (2*sigma);
            else
                Pi(i, j) = data(i, j) < mu - (sqrt(20)*sigma);
                Qi(i, j) = data(i, j) > mu + (sqrt(20)*sigma);
            end
        end
    end
    
    P_sub = sum(Pi, 1);
    Q_sub = sum(Qi, 1);
    
    
    % FIND INDICES OF OUTLIER SUBJECTS
    outliers = ((P_sub+Q_sub)/ size(data, 1) > 0.05) & ((P_sub-Q_sub)/(P_sub+Q_sub) < 0.3);
    
    % REMOVE OUTLIERS
    out = data(:,~outliers);

end
