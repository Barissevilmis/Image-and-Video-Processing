function J = ex2_3_sevilmis(lena, trees, gamma, map_t)

% Normalize Lena
    lena_d = double(lena);
    lena_d = lena_d/255.0;
    
    figure('Name','Ex-2.1.3: Lena - Normalized Image Gamma Correlation');
    subplot(1, size(gamma,2)+1,1);imshow(lena_d);
    title('Normalized Lena');
    
    %Apply gamma correction for Lena
    for j = 1 : size(gamma,2) 
        lena_temp = lena_d;       
        lena_temp = lena_temp.^ gamma(j);
        subplot(1,size(gamma,2)+1,j+1);imshow(lena_temp);
        title(['Gamma Correlated Lena with gamma:=',num2str(gamma(j))]);
    end
    
    % Apply gamma correction for trees using the map
    figure('Name','Ex-2.1.3: Trees - Normalized Image Gamma Correlation');
    subplot(1, size(gamma,2)+1,1);imshow(uint8(trees), map_t);
    title('Normalized Trees');
    for j = 1 : size(gamma,2) 
        map_temp = map_t;
        map_temp = map_temp.^gamma(j);   
        subplot(1,size(gamma,2)+1,j+1);imshow(trees, map_temp);
        title(['Gamma Correlated Trees with gamma:=',num2str(gamma(j))]);
    end
    
end