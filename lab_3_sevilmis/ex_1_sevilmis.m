function [J_l2, J_l1] = ex_1_sevilmis(img, filt_1, filt_2, thr, noise, tmp)
   
    
    % Choose filter type
    filt_name = ""; 
    if filt_1(1,1) == 1/4
        filt_name = "Sobel";
    elseif filt_1(1,1) == 1/3
        filt_name = "Prewitt";
    else
        filt_name = "Roberts";
    end

    % Filter with given template filters in requested direction
    hor = imfilter(img,filt_1,'symmetric','same');
    ver = imfilter(img,filt_2,'symmetric','same');
    
    % L2 and L1 norms
    J_l2 = sqrt(hor.^2 + ver.^2);
    J_l1 = abs(hor) + abs(ver);
    
    % Threshold L2 norm
    J_l2 = uint8(J_l2);
    J_l2 = J_l2 > thr;
    
    %Threshold L1 norm
    J_l1 = uint8(J_l1);
    J_l1 = J_l1 > thr;
    
    if tmp < 3
        figure('Name', 'Ex-1 Figures:');
        subplot(1,2,1);imshow(uint8(ver));title([filt_name,' Horizontal Edge Detection']);
        subplot(1,2,2);imshow(uint8(hor));title([filt_name,' Vertical Edge Detection']);
        %saveas(gcf, ['res/VH-',num2str(tmp),'.png'])
    end
    
    figure('Name', 'Ex-1 Results:');
    subplot(1,2,1);imshow(J_l2);title([filt_name,' L2 with Noise = ', num2str(noise)]);
    subplot(1,2,2);imshow(J_l1);title([filt_name,' L1 with Noise = ', num2str(noise)]);
    %saveas(gcf, ['res/1-',num2str(tmp),'.png'])
end