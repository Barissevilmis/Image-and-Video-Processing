%Lab 3 - Baris Sevilmis - 07/11/2019

%% Exercise 1 - Template Method

%Read and show original images
clear all; close all; clc;
lena_old = imread('lena.png');
rice_old = imread('rice.png');
road_old = imread('road.png');
lena = double(lena_old);
rice = double(rice_old);
road = double(road_old);
figure('Name','Ex-1: Lena, Rice & Road Images');
subplot(1,3,1);imshow(uint8(lena));
title('Lena');
subplot(1,3,2);imshow(uint8(rice));
title('Rice');
subplot(1,3,3);imshow(uint8(road));
title('Road');

%Noise and filter creation
noise = [0, 5, 11, 25];
tmp = 0;

sobel_g1 = (1/4) .* [1 0 -1; 2 0 -2; 1 0 -1];
sobel_g2 = (1/4) .* [-1 -2 -1; 0 0 0; 1 2 1];

prewitt_g1 = (1/3) .* [1 0 -1; 1 0 -1; 1 0 -1];
prewitt_g2 = (1/3) .* [-1 -1 -1; 0 0 0; 1 1 1];

roberts_g1 = [0 0 -1; 0 1 0; 0 0 0];
roberts_g2 = [-1 0 0; 0 1 0; 0 0 0];

g1 = {sobel_g1, prewitt_g1, roberts_g1};
g2 = {sobel_g2, prewitt_g2, roberts_g2};

% Different noise levels: Try each template method with provided filters
for j = 1:4
    lena_new = lena + randn(size(lena,1),size(lena,2))*noise(j);
    lena_new = mat2gray(lena_new) * 255;
    rice_new = rice + randn(size(rice,1),size(rice,2))*noise(j);
    rice_new = mat2gray(rice_new) * 255;
    road_new = road + randn(size(road,1),size(road,2))*noise(j);
    road_new = mat2gray(road_new) * 255;
    for i = 1 : 3
        [lena_res_l2, lena_res_l1] = ex_1_sevilmis(lena_new, g1{i}, g2{i}, 48, noise(j),tmp);
        tmp=tmp+1;
        [rice_res_l2, rice_res_l1] = ex_1_sevilmis(rice_new, g1{i}, g2{i}, 48, noise(j),tmp);
        tmp=tmp+1;
        [road_res_l2, road_res_l1] = ex_1_sevilmis(road_new, g1{i}, g2{i}, 48, noise(j),tmp);
        tmp=tmp+1;
    end

end
%% Exercise 2 - Compass Operator

%Read and show original images
clear all; close all; clc;
lena_old = imread('lena.png');
rice_old = imread('rice.png');
road_old = imread('road.png');
lena = double(lena_old);
rice = double(rice_old);
road = double(road_old);
figure('Name','Ex-2: Lena, Rice & Road Images');
subplot(1,3,1);imshow(uint8(lena));
title('Lena');
subplot(1,3,2);imshow(uint8(rice));
title('Rice');
subplot(1,3,3);imshow(uint8(road));
title('Road');


% Noise creation
noise = [0, 5, 11, 25];
tmp = 0;

% Different noise levels: Try each template method with provided filters
for j = 1:4
    lena_new = lena + randn(size(lena,1),size(lena,2))*noise(j);
    lena_new = mat2gray(lena_new) * 255;
    rice_new = rice + randn(size(rice,1),size(rice,2))*noise(j);
    rice_new = mat2gray(rice_new) * 255;
    road_new = road + randn(size(road,1),size(road,2))*noise(j);
    road_new = mat2gray(road_new) * 255;
    
    ex_2_sevilmis(lena_new, 180, tmp);
    tmp = tmp + 1;
    ex_2_sevilmis(rice_new, 180, tmp);
    tmp = tmp + 1;
    ex_2_sevilmis(road_new, 180, tmp);
    tmp = tmp + 1;
end

%% Exercise 3 - Laplace Operator

%Read and show original images
clear all; close all; clc;
lena_old = imread('lena.png');
rice_old = imread('rice.png');
road_old = imread('road.png');
lena = double(lena_old);
rice = double(rice_old);
road = double(road_old);
figure('Name','Ex-3: Lena, Rice & Road Images');
subplot(1,3,1);imshow(uint8(lena));
title('Lena');
subplot(1,3,2);imshow(uint8(rice));
title('Rice');
subplot(1,3,3);imshow(uint8(road));
title('Road');

% Noise
noise = [0, 5, 11, 25];
tmp = 0;

% Different noise levels: Try each template method with provided filters
for j = 1:4
    lena_new = lena + randn(size(lena,1),size(lena,2))*noise(j);
    lena_new = mat2gray(lena_new) * 255;
    rice_new = rice + randn(size(rice,1),size(rice,2))*noise(j);
    rice_new = mat2gray(rice_new) * 255;
    road_new = road + randn(size(road,1),size(road,2))*noise(j);
    road_new = mat2gray(road_new) * 255;

    ex_3_sevilmis(lena_new, 2.5, 0.7,tmp);
    tmp = tmp + 1;
    ex_3_sevilmis(rice_new, 2.5, 0.7,tmp);
    tmp = tmp + 1;
    ex_3_sevilmis(road_new, 2.5, 0.7,tmp); 
    tmp = tmp + 1;
end
%% Exercise 4 - Frei-Chen Method

%Read and show original images
clear all; close all; clc;
lena_old = imread('lena.png');
rice_old = imread('rice.png');
road_old = imread('road.png');
lena = double(lena_old);
rice = double(rice_old);
road = double(road_old);
figure('Name','Ex-3: Lena, Rice & Road Images');
subplot(1,3,1);imshow(uint8(lena));
title('Lena');
subplot(1,3,2);imshow(uint8(rice));
title('Rice');
subplot(1,3,3);imshow(uint8(road));
title('Road');

% Noise
noise = [0, 5, 11, 25];
tmp = 0;

% Different noise levels: Try each template method with provided filters
for j = 1:4
    lena_new = lena + randn(size(lena,1),size(lena,2))*noise(j);
    lena_new = mat2gray(lena_new) * 255;
    rice_new = rice + randn(size(rice,1),size(rice,2))*noise(j);
    rice_new = mat2gray(rice_new) * 255;
    road_new = road + randn(size(road,1),size(road,2))*noise(j);
    road_new = mat2gray(road_new) * 255;

    ex_4_sevilmis(lena_new, 0.12, tmp);
    tmp = tmp + 1;
    ex_4_sevilmis(rice_new, 0.12, tmp);
    tmp = tmp + 1;
    ex_4_sevilmis(road_new, 0.12, tmp);
    tmp = tmp + 1;
end