function PSNR = ex_PSNRAVG_sevilmis(model, ref)
    
    R = abs(double(model(:,1) - ref(:,1)));
    G = abs(double(model(:,2) - ref(:,2)));
    B = abs(double(model(:,3) - ref(:,3)));
    
    MSE_R = mean(R.^2);
    MSE_G = mean(G.^2);
    MSE_B = mean(B.^2);
    
    PSNR_R = 10*log10((255^2)/MSE_R);
    PSNR_G = 10*log10((255^2)/MSE_G);
    PSNR_B = 10*log10((255^2)/MSE_B);
    
    PSNR = (PSNR_R + PSNR_G + PSNR_B) / 3;
    
end