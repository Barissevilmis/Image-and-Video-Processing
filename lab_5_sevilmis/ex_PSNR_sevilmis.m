function J = ex_PSNR_sevilmis(model, ref)
    
    dist = abs(double(model - ref));
    val = dist(:,1).^2 + dist(:,2).^2 + dist(:,3).^2;
    MSE = mean(val.^2);
    
    J = 10*log10((255^2)/MSE);
end