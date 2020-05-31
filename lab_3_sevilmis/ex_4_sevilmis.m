function J = ex_4_sevilmis(img, thr, tmp)

    % All the filters entered manually
    f0=(1/(2*sqrt(2)))*[1 sqrt(2) 1; 0 0 0; -1 -sqrt(2) -1];
    f1=(1/(2*sqrt(2)))*[1 0 -1; sqrt(2) 0 -sqrt(2); 1 0 -1];
    f2=(1/(2*sqrt(2)))*[0 -1 sqrt(2); 1 0 -1;-sqrt(2) 1 0];
    f3=(1/(2*sqrt(2)))*[sqrt(2) -1 0; -1 0 1;0 1 -sqrt(2)];
    f4=(1/2)*[0 1 0;-1 0 -1;0 1 0];
    f5=(1/2)*[-1 0 1; 0 0 0; 1 0 -1];
    f6=(1/6)*[1 -2 1;-2 4 -2;1 -2 1];
    f7=(1/6)*[-2 1 -2;1 4 1;-2 1 -2];
    f8=(1/3)*ones(3,3);
    
    % Project image onto filter subspaces
    imgf0 = imfilter(img, f0 , 'symmetric', 'same');
    imgf1 = imfilter(img, f1 , 'symmetric', 'same');
    imgf2 = imfilter(img, f2 , 'symmetric', 'same');
    imgf3 = imfilter(img, f3 , 'symmetric', 'same');
    imgf4 = imfilter(img, f4 , 'symmetric', 'same');
    imgf5 = imfilter(img, f5 , 'symmetric', 'same');
    imgf6 = imfilter(img, f6 , 'symmetric', 'same');
    imgf7 = imfilter(img, f7 , 'symmetric', 'same');
    imgf8 = imfilter(img, f8 , 'symmetric', 'same');
    
    img_arr = {imgf0 imgf1 imgf2 imgf3 imgf4 imgf5 imgf6 imgf7 imgf8};
    m = zeros(size(img));
    s = zeros(size(img));
    J = zeros(size(img));
    
    % Find m and s
    for i = 1 : 9
        if i < 3
            m = m + img_arr{i}.^2; 
        end
        s = s + img_arr{i}.^2;
        
    end
    
    % Find magnitude
    filtered_img = sqrt(m./s);
  
    % Mark edge points
    indices_w = (filtered_img >= thr);
    J(indices_w) = 255.0;
    
    J = uint8(J);
    figure('Name','Frei-Chen Method');
    imshow(J);title(['Frei-Chen with Thr.',num2str(thr)]);
    %saveas(gcf, ['res/4-',num2str(tmp),'.png'])


end