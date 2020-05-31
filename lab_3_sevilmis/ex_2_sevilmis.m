function J = ex_2_sevilmis(img ,thr, tmp)
    [row, col, ch] = size(img);

    % Create an empty kirsch operator matrices
    kirsch = zeros(3,3,7);
    kirsch(:,:,1) = [-3 -3 5; -3 0 5; -3 -3 5];
    kirsch(:,:,2) = [-3 5 5; -3 0 5; -3 -3 -3];
    
    filtered_img = zeros(row, col, 7);

    % Rotate counterclockwise by 90 degrees according to the i - 2: Same with
    % 45 degrees i - 1
    for i=3:7
        kirsch(:,:,i) = rot90(kirsch(:,:,i-2));
    end
        
    % Apply kirsch
    for m = 1 : 7
        filtered_img(:,:,m) = imfilter(img,kirsch(:,:,m),'symmetric','same');
    end
 
    % Choose maximum cooridnates from every dimension after the absolute
    J = max(abs(filtered_img), [], 3);
    
    % Threshold
    J = uint8(J);
    J = J > thr;
    
    figure('Name', 'Kirsch Operator');
    imshow(J);title(['Kirsch with Thr. = ',num2str(thr)]);
    
    %saveas(gcf, ['res/2-',num2str(tmp),'.png'])

end