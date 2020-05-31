function PSNR = ex_PSNRYUV_sevilmis(model, ref)
    
    [MY, MU, MV] = rgb2yuv(model(:,1), model(:,2), model(:,3));
    [RY, RU, RV] = rgb2yuv(ref(:,1), ref(:,2), ref(:,3));
    model_yuv = [MY MU MV];
    ref_yuv = [RY RU RV];
    Y = abs(double(model_yuv(:,1) - ref_yuv(:,1)));
    U = abs(double(model_yuv(:,2) - ref_yuv(:,2)));
    V = abs(double(model_yuv(:,3) - ref_yuv(:,3)));
    
    MSE_Y = mean(Y.^2);
    MSE_U = mean(U.^2);
    MSE_V = mean(V.^2);
    
    PSNR_Y = 10*log10((255^2)/MSE_Y);
    PSNR_U = 10*log10((255^2)/MSE_U);
    PSNR_V = 10*log10((255^2)/MSE_V);
    
    PSNR = (6*PSNR_Y + PSNR_U + PSNR_V) / 8;
    
end