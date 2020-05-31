function img = ex_3_sevilmis(img, threshold_mat)

    %Create a filter in requested format:
    %If/Else to create a filter same size as the image
    N = double(size(threshold_mat,1));
    N2 = double(size(threshold_mat,2));
    if (size(img,1)/N ~= floor(size(img,1)/N))||(size(img,2)/N2 ~= floor(size(img,2)/N2))
        S = repmat(threshold_mat , floor(size(img,1)/N), floor(size(img,2)/N2));
        row = ceil((size(img,1)/N - floor(size(img,1)/N))*N);
        col = ceil((size(img,2)/N2 - floor(size(img,2)/N2))*N2);
        
        if (col + (floor(size(img,2)/N2)*N2) ~= size(img,2))&&(row + (floor(size(img,1)/N)*N) ~= size(img,1))
            S = [S S(:,1:col-1); S(1:row-1,:) S(1:row-1,1:col-1)]; 
        elseif (col + (floor(size(img,2)/N2)*N2) ~= size(img,2))
            S = [S S(:,1:col-1); S(1:row,:) S(1:row,1:col-1)];
        elseif (row + (floor(size(img,1)/N)*N) ~= size(img,1))
            S = [S S(:,1:col); S(1:row-1,:) S(1:row-1,1:col)];
        else
            S = [S S(:,1:col); S(1:row,:) S(1:row,1:col)]; 
        end
        
    else
        S = repmat(threshold_mat , size(img,1)/N, size(img,2)/N2);
    end
    img(img > S) = 255.0;
    img(img < S) = 0.0;
    
end