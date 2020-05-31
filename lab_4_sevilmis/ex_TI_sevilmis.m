function TI = ex_TI_sevilmis(data)

    s1 = floor(data.FrameRate*data.Duration);
    Mn = zeros(s1,1);
    
    i = 1;
    
    % READ DATA FROM 0
    % FIND DIFFERENCE BETWEEN CONSECUTIVE LUMINANCE PLANES
    % TAKE STD
    % TAKE MAX STD
    data.CurrentTime = 0;
    prev = readFrame(data);
    prevR = prev(:,:,1);
    prevG = prev(:,:,2);
    prevB = prev(:,:,3);
    [prevY,~,~] = rgb2yuv(prevR,prevG,prevB);
    while hasFrame(data)
        curr = readFrame(data);
        R = curr(:,:,1);
        G = curr(:,:,2);
        B = curr(:,:,3);
        [Y,~,~] = rgb2yuv(R,G,B);
        res = double(Y) - double(prevY);
        Mn(i) = std2(res);
        prevY = Y;
        i = i + 1;        
    end
    
    TI = max(Mn);
    
end
