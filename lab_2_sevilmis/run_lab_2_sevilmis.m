%Lab 2 - Baris Sevilmis 23/10/2019

%% Exercise 1 - Fixed Threshold Method
%Read and show original images
clear all; close all; clc;
lena_old = imread('lena-y.png');
wool_old = imread('wool.png');
lena = double(lena_old);
wool = double(wool_old);
figure('Name','Ex-1: Lena-y & Wool Image');
subplot(1,2,1);imshow(uint8(lena));
title('Lena');
subplot(1,2,2);imshow(uint8(wool));
title('Wool');

%Thresholding
lena_range_bl = lena < 128.0;
lena_range_wh = lena >= 128.0;
lena(lena_range_bl) = 0.0;
lena(lena_range_wh) = 255.0;

wool_range_bl = wool < 128.0;
wool_range_wh = wool >= 128.0;
wool(wool_range_bl) = 0.0;
wool(wool_range_wh) = 255.0;

figure('Name','Ex-1: Lena-y & Wool Fixed Threshold with Th. = 128 Image');
subplot(1,2,1);imshow(uint8(lena));
title('Fixed Thresholded Lena');
subplot(1,2,2);imshow(uint8(wool));
title('Fixed Thresholded Wool');

%MSE
N1 = size(lena_old,1) * size(lena_old,2);
N2 = size(wool_old,1) * size(wool_old,2);
mse_lena = sum((lena - double(lena_old)).^2, 'all')/N1;
mse_wool= sum((wool - double(wool_old)).^2, 'all')/N2;
%% Exercise 2 - Random Threshold Method
%Read and show images
clear all; close all; clc;
lena_old = imread('lena-y.png');
wool_old = imread('wool.png');
lena = double(lena_old);
wool = double(wool_old);
figure('Name','Ex-2: Lena-y & Wool Image');
subplot(1,2,1);imshow(uint8(lena));
title('Lena');
subplot(1,2,2);imshow(uint8(wool));
title('Wool');

%Noise amplitude variation
for i = 10 : 40 : 250
    %noise_lena = -i + (i+i)*rand(size(lena,1),size(lena,2));
    %noise_wool = -i + (i+i)*rand(size(wool,1),size(wool,2));
    noise_lena = unidrnd(i,size(lena)) - i/2;
    noise_wool = unidrnd(i,size(wool)) - i/2;
    lena_tmp = noise_lena + lena;
    wool_tmp = noise_wool + wool;
    
    %Thresholding
    lena_range_bl = lena_tmp < 128.0;
    lena_range_wh = lena_tmp >= 128.0;
    lena_tmp(lena_range_bl) = 0.0;
    lena_tmp(lena_range_wh) = 255.0;

    wool_range_bl = wool_tmp < 128.0;
    wool_range_wh = wool_tmp >= 128.0;
    wool_tmp(wool_range_bl) = 0.0;
    wool_tmp(wool_range_wh) = 255.0;

    figure('Name','Ex-2: Lena-y & Wool Random Threshold with Th. = 128 Image');
    subplot(1,2,1);imshow(uint8(lena_tmp));
    title(['Fixed Thresholded Lena, amplitude of noise: ',num2str(i)]);
    subplot(1,2,2);imshow(uint8(wool_tmp));
    title(['Fixed Thresholded Wool, amplitude of noise: ',num2str(i)]);
    %MSE for a noise amplitude of 10
    if i == 10
        N1 = size(lena_old,1) * size(lena_old,2);
        N2 = size(wool_old,1) * size(wool_old,2);
        mse_lena = sum((lena_tmp - double(lena_old)).^2, 'all')/N1;
        mse_wool= sum((wool_tmp - double(wool_old)).^2, 'all')/N2;
    end
end


%% Exercise 3: Ordered Threshold Method
%Read and show
clear all; close all; clc;
lena_old = imread('lena-y.png');
wool_old = imread('wool.png');
lena = double(lena_old);
wool = double(wool_old);
figure('Name','Ex-3: Lena-y & Wool Image');
subplot(1,2,1);imshow(uint8(lena));
title('Lena');
subplot(1,2,2);imshow(uint8(wool));
title('Wool');

S = [34 29 17 21 30 35; 28 14 9 16 20 31; 13 8 4 5 15 19; 12 3 0 1 10 18; 27 7 2 6 23 24; 33 26 11 22 25 32];

%Quantize
N = max(S(:))+2;
S = S + 0.5;
S_map = S * (256/N);

%Threshold
res_lena = ex_3_sevilmis(lena, S_map);
res_wool = ex_3_sevilmis(wool, S_map);

figure('Name','Ex-3: Lena-y & Wool Ordered Threshold');
subplot(1,2,1);imshow(uint8(res_lena));
title('Ordered Thresholded Lena');
subplot(1,2,2);imshow(uint8(res_wool));
title('Ordered Thresholded Wool');

%MSE
N1 = size(lena_old,1) * size(lena_old,2);
N2 = size(wool_old,1) * size(wool_old,2);
mse_lena = sum((res_lena - double(lena_old)).^2, 'all')/N1;
mse_wool= sum((res_wool - double(wool_old)).^2, 'all')/N2;

%% Exercise 4: Ordered Matrix With Centered Points
%Read and show
clear all; close all; clc;
lena_old = imread('lena-y.png');
wool_old = imread('wool.png');
lena = double(lena_old);
wool = double(wool_old);
figure('Name','Ex-4: Lena-y & Wool Image');
subplot(1,2,1);imshow(uint8(lena));
title('Lena');
subplot(1,2,2);imshow(uint8(wool));
title('Wool');

%Filters
C6 = [34 25 21 17 29 33; 30 13 9 5 12 24; 18 6 1 0 8 20; 22 10 2 3 4 16; 26 14 7 11 15 28; 35 31 19 23 27 32];
E6 = [30 22 16 21 33 35; 24 11 7 9 26 28; 13 5 0 2 14 19; 15 3 1 4 12 18; 27 8 6 10 25 29; 32 20 17 23 31 34];

%Quantize
N1 = max(C6(:))+2;
C6 = C6 + 0.5;
C6_map = C6 * (256/N1);

N2 = max(E6(:))+2;
E6 = E6 + 0.5;
E6_map = E6 * (256/N2);

%Threshold
c6_lena = ex_3_sevilmis(lena, C6_map);
c6_wool = ex_3_sevilmis(wool, C6_map);

e6_lena = ex_3_sevilmis(lena, E6_map);
e6_wool = ex_3_sevilmis(wool, E6_map);

figure('Name','Ex-4: Lena-y & Wool Ordered Threshold with C6');
subplot(1,2,1);imshow(uint8(c6_lena));
title('Ordered Thresholded Lena - C6');
subplot(1,2,2);imshow(uint8(e6_lena));
title('Ordered Thresholded Lena - E6');

figure('Name','Ex-4: Lena-y & Wool Ordered Threshold with E6');
subplot(1,2,2);imshow(uint8(c6_wool));
title('Ordered Thresholded Wool - C6');
subplot(1,2,1);imshow(uint8(e6_wool));
title('Ordered Thresholded Wool - E6');

N1 = size(lena_old,1) * size(lena_old,2);
N2 = size(wool_old,1) * size(wool_old,2);

mse_lena = sum((c6_lena - double(lena_old)).^2, 'all')/N1;
mse_wool= sum((c6_wool - double(wool_old)).^2, 'all')/N2;

mse2_lena = sum((e6_lena - double(lena_old)).^2, 'all')/N1;
mse2_wool= sum((e6_wool - double(wool_old)).^2, 'all')/N2;

%% Exercise 5 - Diagonal Ordered Matrix with Balanced Centered Points
%Read and Show
clear all; close all; clc;
lena_old = imread('lena-y.png');
wool_old = imread('wool.png');
lena = double(lena_old);
wool = double(wool_old);
figure('Name','Ex-5: Lena-y & Wool Image');
subplot(1,2,1);imshow(uint8(lena));
title('Lena');
subplot(1,2,2);imshow(uint8(wool));
title('Wool');

%Filters
O81 = [13 9 5 12; 6 1 0 8; 10 2 3 4; 14 7 11 15];
O82 = [18 22 26 19; 25 30 31 23; 21 29 28 27; 17 24 20 16];

O8 = [O81 O82; O82 O81];

N = max(O8(:))+2;
O8 = O8 + 0.5;
O8_map = O8 * (256/N);

O8_diag = O8(:,2:8);
M = size(O8_diag,1)/2 -1;
for i = 1 : M
    O8_diag(1:i,M-i+1) = 0;
    O8_diag(1:i,M+i+1) = 0;
    O8_diag(2*M+i-1:2*(M+1),i) = 0; 
    O8_diag(2*M+i-1:2*(M+1),2*(M+1)-i) = 0; 
end

O8_diag = O8_diag * (256/N);

%Threshold
o8_lena = ex_3_sevilmis(lena, O8_map);
o8_wool = ex_3_sevilmis(wool, O8_map);

o8_lena2 = ex_3_sevilmis(lena, O8_diag);
o8_wool2 = ex_3_sevilmis(wool, O8_diag);

figure('Name','Ex-5: Lena-y & Wool Diagonal Ordered Threshold with O8');
subplot(1,2,1);imshow(uint8(o8_lena));
title('Ordered Thresholded Lena - O8');
subplot(1,2,2);imshow(uint8(o8_wool));
title('Ordered Thresholded Lena - O8');

figure('Name','Ex-5: Lena-y & Wool Diagonal Ordered Threshold with O8 - Diagonal');
subplot(1,2,1);imshow(uint8(o8_lena2));
title('Ordered Thresholded Lena - O8-Diag');
subplot(1,2,2);imshow(uint8(o8_wool2));
title('Ordered Thresholded Lena - O8-Diag');

%MSE
N1 = size(lena_old,1) * size(lena_old,2);
N2 = size(wool_old,1) * size(wool_old,2);

mse_lena = sum((o8_lena - double(lena_old)).^2, 'all')/N1;
mse_wool = sum((o8_wool - double(wool_old)).^2, 'all')/N2;

%% Exercise 6 - Ordered Matrix with Dispersed Dots
%Read and Show
clear all; close all; clc;
lena_old = imread('lena-y.png');
wool_old = imread('wool.png');
lena = double(lena_old);
wool = double(wool_old);
figure('Name','Ex-6: Lena-y & Wool Image');
subplot(1,2,1);imshow(uint8(lena));
title('Lena');
subplot(1,2,2);imshow(uint8(wool));
title('Wool');

%Filters
D2 = [0 2; 3 1];
D3 = [8 4 5; 3 0 1; 7 2 6];
U2 = ones(2,2);
U3 = ones(3,3);

D6 = [4*D3 (4*D3 + 2*U3); (4*D3 + 3*U3) (4*D3 + U3)];
D4 = [4*D2 (4*D2 + 2*U2); (4*D2 + 3*U2) (4*D2 + U2)];

N = max(D6(:))+2;
D6 = D6 + 0.5;
D6_map = D6 * (256/N);

N = max(D4(:))+2;
D4 = D4 + 0.5;
D4_map = D4 * (256/N);

%Threshold
d6_lena = ex_3_sevilmis(lena, D6_map);
d6_wool = ex_3_sevilmis(wool, D6_map);

d4_lena = ex_3_sevilmis(lena, D4_map);
d4_wool = ex_3_sevilmis(wool, D4_map);

figure('Name','Ex-6: Lena-y & Wool Ordered Matrix with Dispersed Dots D6');
subplot(1,2,1);imshow(uint8(d6_lena));
title('Ordered Sparse Thresholded Lena - D6');
subplot(1,2,2);imshow(uint8(d6_wool));
title('Ordered Sparse Thresholded Lena - D6');

figure('Name','Ex-6: Lena-y & Wool Ordered Matrix with Dispersed Dots D4');
subplot(1,2,1);imshow(uint8(d4_lena));
title('Ordered Sparse Thresholded Lena - D4');
subplot(1,2,2);imshow(uint8(d4_wool));
title('Ordered Sparse Thresholded Lena - D4');

%MSE
N1 = size(lena_old,1) * size(lena_old,2);
N2 = size(wool_old,1) * size(wool_old,2);

mse_lena = sum((d6_lena - double(lena_old)).^2, 'all')/N1;
mse_wool= sum((d6_wool - double(wool_old)).^2, 'all')/N2;

mse2_lena = sum((d4_lena - double(lena_old)).^2, 'all')/N1;
mse2_wool= sum((d4_wool - double(wool_old)).^2, 'all')/N2;

%% Exercise 7 - Error Diffusion Method
%Read and Show
clear all; close all; clc;
lena_old = imread('lena-y.png');
wool_old = imread('wool.png');
lena = double(lena_old);
wool = double(wool_old);
lena = double(lena)/255;
wool = double(wool)/255;
figure('Name','Ex-7: Lena-y & Wool Image');
subplot(1,2,1);imshow(lena);
title('Lena');
subplot(1,2,2);imshow(wool);
title('Wool');

%Filters
floyd_steinberg = (1/16) * [0 0 0; 0 0 7; 3 5 1];
stucki = (1/42) * [0 0 0 0 0; 0 0 0 0 0; 0 0 0 8 4; 2 4 8 4 2; 1 2 4 2 1];

%Threshold
res_lena = ex_7_sevilmis(lena, floyd_steinberg);
res_wool = ex_7_sevilmis(wool, floyd_steinberg);

res_lena2 = ex_7_sevilmis(lena, stucki);
res_wool2 = ex_7_sevilmis(wool, stucki);

figure('Name','Ex-7: Lena-y & Wool Error Diffusion Method Floyd-Steinberg');
subplot(1,2,1);imshow(res_lena);
title('Error Diffused Lena - Floyd-Steinberg');
subplot(1,2,2);imshow(res_wool);
title('Error Diffused Wool - Floyd-Steinberg');

figure('Name','Ex-7: Lena-y & Wool Error Diffusion Method Stucki');
subplot(1,2,1);imshow(res_lena2);
title('Error Diffused Lena - Stucki');
subplot(1,2,2);imshow(res_wool2);
title('Error Diffused Wool - Stucki');

%MSE
N1 = size(lena_old,1) * size(lena_old,2);
N2 = size(wool_old,1) * size(wool_old,2);

mse_lena = sum(((res_lena.*255.0) - double(lena_old)).^2, 'all')/N1;
mse_wool= sum(((res_wool.*255.0) - double(wool_old)).^2, 'all')/N2;

mse2_lena = sum(((res_lena2.*255.0) - double(lena_old)).^2, 'all')/N1;
mse2_wool= sum(((res_wool2.*255.0) - double(wool_old)).^2, 'all')/N2;