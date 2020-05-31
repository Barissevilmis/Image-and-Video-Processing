function J = ex_3_sevilmis(img, sigma, thr, tmp)
    [row, col, ch] = size(img);
    
    % Gaussian filter with specified sigma
    log_filter = fspecial('log', 2*ceil(3*sigma)+1, sigma);
    
    % Gaussian filter for noise reduction
    filtered_img = imfilter(img,log_filter,'symmetric','same');
    J = zeros(size(filtered_img));
    
    
    % Zero crossing: Opposite signs + Check if difference with a specific neighbor > Threshold
    % If yes, then negative pixel is an edge point
    for i = 2 : row
        for j = 2 : col
            tmp_row_up = (filtered_img(i-1,j) < 0) & (filtered_img(i,j) >= 0);
            tmp_row_down = (filtered_img(i-1,j) >= 0) & (filtered_img(i,j) < 0);
            tmp_col_right = (filtered_img(i,j-1) >= 0) & (filtered_img(i,j) < 0);
            tmp_col_left = (filtered_img(i,j-1) < 0) & (filtered_img(i,j) >= 0);

            if tmp_row_up && ((filtered_img(i,j) - filtered_img(i-1,j))>= thr)
                J(i-1,j) = 255.0;
            elseif  tmp_row_down && ((filtered_img(i-1,j) - filtered_img(i,j))>= thr)
                J(i,j) = 255.0;
            end     
            if tmp_col_right && ((filtered_img(i,j-1) - filtered_img(i,j)) >= thr)
                J(i,j) = 255.0;
            elseif tmp_col_left && ((filtered_img(i,j) - filtered_img(i,j-1)) >= thr)
                J(i,j-1) = 255.0;
            end
        end
    end
    
    J = uint8(J);
    figure('Name', 'Laplace Operator');
    imshow(J);title(['Laplace with Thr. = ', num2str(thr)]);
    %saveas(gcf, ['res/3-',num2str(tmp),'.png'])
end