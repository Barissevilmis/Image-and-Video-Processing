function SI = ex_SI_sevilmis(data)
    
    s1 = floor(data.FrameRate*data.Duration);
    Fn = zeros(s1,1);
    
    % READ VIDEO FROM 0
    % APPLY SOBEL TO LUMINANCE PLANE
    % TAKE STD
    % CHOOSE MAX STD
    i = 1;
    data.CurrentTime = 0;
    while hasFrame(data)
        curr = readFrame(data);
        R = curr(:,:,1);
        G = curr(:,:,2);
        B = curr(:,:,3);
        [Y,~,~] = rgb2yuv(R,G,B);
        sobel = fspecial('sobel');
        res = imfilter(double(Y), sobel);
        Fn(i) = std2(res);
        i = i + 1;
    end
    SI = max(Fn);
   
end
