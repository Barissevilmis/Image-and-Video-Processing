function [MOS, CI] = ex_MOS_sevilmis(data)
    
    
    N = size(data, 2);
    
    %ALPHA VALUES AND TWO TAILED T DISTRIBUTION
    %FIND STD, USE TO CALCULATE CI
    %MEAN FOR MOS
    alpha = 0.05;
    alphaup = 1-alpha/2;
    t = tinv(alphaup,N-1);
    
    CI = zeros(size(data, 1), 1);
    MOS = zeros(size(data, 1), 1);
    
    for i = 1 : size(data, 1)
        sigma = std(data(i,:));
        CI(i) = t * (sigma/sqrt(N));
        MOS(i) = sum(data(i,:))/ N;
    end
    

end