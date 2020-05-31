function [MSE, HAU] = ex_PointToPoint_sevilmis(model, ref)
    
    dist = model - ref;
    val = sqrt(dist(:,1).^2 + dist(:,2).^2 + dist(:,3).^2);
    MSE = mean(val.^2);
    HAU = max(val);
    
end