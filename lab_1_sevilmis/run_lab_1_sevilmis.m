% Lab 1 - Baris Sevilmis - 25.09.2019
 %% Exercise 2.1.1 - Read Images and Show Images
 % Display Grayscale and Original images
[trees, map_t] = imread('trees.tif');
lena = imread('lena.tif');
figure('Name','Ex-2.1.1: Read and Show-Trees');
subplot(1,2,1);imshow(trees, map_t);
title('Original Trees');
subplot(1,2,2);imshow(trees, rgb2gray(map_t));
title('Grayscale Trees');
figure('Name','Ex-2.1.1: Read and Show-Lena');
subplot(1,2,1);imshow(lena);
title('Original Lena');

[size_t, dim_t] = size(map_t);
[row_l, col_l, dim_l] = size(lena);

%RGB Lena to Grayscale lena
if dim_l > 1
    lena_gr = rgb2gray(lena);
end

%RGB Lena to Grayscale trees
if dim_t > 1
    map_gr = rgb2gray(map_t);
end

subplot(1,2,2);imshow(lena_gr);
title('Grayscale Lena');
lena_gr = double(lena_gr);

    
%% Exercise 2.1.2 - Show Images in Gray level- Invert them
[trees, map_t] = imread('trees.tif');
lena = imread('lena.tif');

[size_t, dim_t] = size(map_t);
[row_l, col_l, dim_l] = size(lena);

%RGB Lena to Grayscale lena
if dim_l > 1
    lena_gr = rgb2gray(lena);
end

%RGB Lena to Grayscale trees
if dim_t > 1
    map_gr = rgb2gray(map_t);
end
lena_gr = double(lena_gr);

%Invert Lena
lena_inv = double(255) - lena_gr;
figure('Name','Ex-2.1.2: Lena - Grayscale and Invert');
subplot(1,2,1);imshow(uint8(lena_gr));
title('Grayscale Lena');
subplot(1,2,2);imshow(uint8(lena_inv));
title('Inverted Lena');

%Invert Trees
figure('Name','Ex-2.1.2: Trees - Grayscale and Invert');
map_t_inv = double(1) - map_gr;
subplot(1,2,1);imshow(trees, map_gr);
title('Grayscale Trees');
subplot(1,2,2);imshow(trees, map_t_inv);
title('Inverted Trees');

%% Exercise 2.1.3 - Function for Modifying Color Tables
[trees, map_t] = imread('trees.tif');
lena = imread('lena.tif');

gamma_val = [0.5 1 1.5 2.0];
ex2_3_sevilmis(lena, trees, gamma_val, map_t);

%% Exercise 2.1.4 - Image of a Chess Board

% True Colors
pattern_yellow = [255 0; 0 255];
pattern_blue = [0 255; 255 0];

%CHess using RGB
chess_true = zeros(8,8,3);
chess_true(:,:,1) = repmat(pattern_yellow, [4 4]);
chess_true(:,:,2) = repmat(pattern_yellow, [4 4]);
chess_true(:,:,3) = repmat(pattern_blue, [4 4]);

% Indexed Version of Chess using Color map: Give colors to specific(2)
% intensity values according to the patterns
chess_indexed = zeros(8,8);
map_indexed = zeros(256,3);
chess_indexed(:,:) = repmat(pattern_yellow, [4 4]);
map_indexed(255,:) = [1.0 1.0 0.0];
map_indexed(1,:) = [0.0 0.0 1.0];


figure('Name','Ex-2.1.4: Chess');
subplot(2,2,1);imshow(chess_true);
title('True Color Chess');
subplot(2,2,2);imshow(chess_indexed, map_indexed);
title('Indexed Chess');

% Save images
imwrite(chess_true, "chess_truecolor.tif");
imwrite(chess_indexed,map_indexed, "chess_indexed.tif");

%% Exercise 2.2


lena_y = imread('lena-y.png');
[row_y, col_y, dim_y] = size(lena_y);

%Greyscale Lena
if  dim_y > 1
    lena_y_gr = rgb2gray(lena_y);
else
    lena_y_gr = lena_y;
end

figure('Name','Ex-2.2');
subplot(1,2,1);imshow(lena_y);
title('Original Lena');
subplot(1,2,2);imshow(lena_y_gr);
title('Grayscale Lena');

%Step size and result array
lena_y_gr = double(lena_y_gr);
step_size = [128.0 64.0 32.0 16.0 8.0 4.0 2.0];
lena_cell = cell(length(step_size),1);

%Quantization
%Array of images being used
figure('Name','Ex-2.2: Quantization');
subplot(2,4,1);imshow(uint8(lena_y_gr));
title('Non-quantized Lena');
for i = 1 : length(step_size)
    % Quantization formula as given in the sheet
    lena_cell{i} = double(floor(lena_y_gr ./ step_size(i)));
    lena_cell{i} = uint8(lena_cell{i} .* step_size(i) + step_size(i)/2.0);
    subplot(2,4,i+1);imshow(lena_cell{i});
    title(['Quantization with gray level of ',num2str(step_size(length(step_size) + 1 - i))]);
end

%% Exercise 2.3

% Filter Def.:
% Create 2D filter by 1D * 1D.T: Low Pass Filter
filter = [0.0357; 0.2411; 0.4464; 0.2411; 0.0357];
filter_2D = filter * transpose(filter);

figure('Name','Ex-2.3: Frequency Response of Filter Part - I');
filter_freq = freqz2(filter_2D);
mesh(filter_freq);

% Filter Def.:
% New filter: Edge Sharpening Filter
filter_sec2D = (1/6) .* [-1 -4 -1; -4 26 -4; -1 -4 -1];

figure('Name','Ex-2.3: Frequency Response of Filter Part - II');
filter_secfreq = freqz2(filter_sec2D);
mesh(filter_secfreq);

% Read the Image
gold = imread('gold-text.png');
[row_g, col_g, dim_g] = size(gold);

if dim_g > 1
    gold = rgb2gray(gold);
end
gold = double(gold);

% Convolve images using filters: Image -> Filter 1 -> Result 1 -> Filter 2
% -> Result 2
figure('Name','Ex-2.3: Convolution Results');
gold_conv1 = conv2(gold, filter_2D);
gold_conv2 = conv2(gold_conv1, filter_sec2D);

subplot(1,3,1);imshow(uint8(gold));
title('Original Gold-Text');
subplot(1,3,2);imshow(uint8(gold_conv1));
title('Convolved Gold-Text with Filter - I');
subplot(1,3,3);imshow(uint8(gold_conv2));
title('Convolved Gold-Text with Filter - II after Filter - I');

% Frequency Plots
figure('Name','Ex-2.3: Frequency Responses');
subplot(1,3,1);mesh(uint8(gold));
title('Plot of the Gold-Text');
subplot(1,3,2);mesh(uint8(gold_conv1));
title('Plot of the convolved Gold-Text with Filter - I');
subplot(1,3,3);mesh(uint8(gold_conv2));
title('Plot of the convolved Gold-Text with Filter - II');

%% Exercise 2.4

% Read the Image
gold = imread('gold-text.png');
[row_g, col_g, dim_g] = size(gold);

if dim_g > 1
    gold = rgb2gray(gold);
end
gold = double(gold);

% Read the Letter
gl = imread('g-letter.png');
[row_gl, col_gl, dim_gl] = size(gl);

if dim_gl > 1
    gl = rgb2gray(gl);
end
gl = double(gl);

figure('Name','Ex-2.4: Spatial Domain results(red) and Frequency Domain results(green)');
imshow(uint8(gold));
hold on;
noise_std = [0.0 5.0 10.0 25.0 40.0 50.0];

% Start with 0 noise
for j = 1 : size(noise_std,2)
    
    % Add noise randomly with requested noise std. deviation
    noise = noise_std(j) * randn(size(gold,1),size(gold,2));
    gold_new = gold + noise;
    
    % Zero nominal averages
    gold_zero = gold_new - 128;
    gl_zero = gl - 128;
    
    % Max correlation using correlation operator
    gold_corrGL = xcorr2(gold_zero, gl_zero);
    gold_corr = max(max(gold_corrGL)); 

    [gold_row, gold_col] = find(gold_corrGL == gold_corr);
    
    %Spatial Domain Result
    plot(gold_col, gold_row, 'r*','markersize',20,'LineWidth',10);

    % Frequency domain
    gold_fft = fft2(gold_zero);
    gl_fft = fft2(gl_zero, row_g, col_g);

    % Max correlation
    gold_corrfftGL = ifft2(conj(gl_fft) .* gold_fft);
    gold_corrfft = max(max(gold_corrfftGL)); 

    % Apply translation factor
    [gold_rowfft, gold_colfft] = find(gold_corrfftGL == gold_corrfft);
    gold_rowfft = gold_rowfft + size(gl_zero, 1) - 1;
    gold_colfft = gold_colfft + size(gl_zero, 2) - 1;
    
    % Frequency Domain Result
    % hold on;
    plot(gold_colfft, gold_rowfft,'g+','markersize',20,'LineWidth',10);
end


%% Exercise 2.5
% Resampling
figure('Name','Ex-2.5: Resampling)');
% Read the Image
sub4 = imread('sub4.tif');
[row_4, col_4, dim_4] = size(sub4);

if dim_4 > 1
    sub4_gr = rgb2gray(sub4);
else
    sub4_gr = sub4;
end
sub4_gr = double(sub4_gr);

% Downsample by 2 -> Image 1 -> Result 1 -> Downsample by  4 -> Result 2
sub4_grdown = sub4_gr(1:2:end,1:2:end);
sub4_grdown2 = sub4_grdown(1:4:end,1:4:end);
subplot(1,3,1);imshow(uint8(sub4_gr));
title('Greyscale Sub4');
subplot(1,3,2);imshow(uint8(sub4_grdown));
title('Downsampled Greyscale Sub4 - Factor of 2');
subplot(1,3,3);imshow(uint8(sub4_grdown2));
title('Downsampled Resulting Image - Factor of 4');

%% Exercise 2.6

% Read the Image
lena_y = imread('lena-y.png');
[row_ly, col_ly, dim_ly] = size(lena_y);
imag_disp = 1j;

if dim_ly > 1
    lena_y_gr = rgb2gray(lena_y);
else
    lena_y_gr = lena_y;
end
lena_y_gr = double(lena_y_gr);

% Part-I:
lena_y_fft = fft2(lena_y_gr);
% Real Part
lena_y_real = ifft2(real(lena_y_fft));
% Imaginary Part
lena_y_imag = ifft2(imag(lena_y_fft)*imag_disp);

% Display Lena and Fourier Domain
figure('Name','Ex-2.6: Original Image and Frequency Domian - Part I)');
subplot(1,2,1);imshow(uint8(lena_y_gr));
title('Greyscale Lena-y');
subplot(1,2,2);imagesc(log2(abs(fftshift(lena_y_fft))));
title('Lena-y in Fourier Domain with log2(abs(fftshift()))');

% Real Lena and Imaginary Lena
figure('Name','Ex-2.6: Real and Imaginary  Part of the 2DFT - Part I - Wrong range)');
subplot(1,2,1);imshow(lena_y_real,[]);
title('Lena-y after fft and ifft with only Real part');
subplot(1,2,2);imshow(uint8(ifft2(imag(lena_y_fft)*imag_disp)));
title('Lena-y after fft and ifft with only Imaginary part');

figure('Name','Ex-2.6: Real and Imaginary  Part of the 2DFT - Part I)');
subplot(1,2,1);imshow(lena_y_real,[]);
title('Lena-y after fft and ifft with only Real part');
subplot(1,2,2);imshow(lena_y_imag,[]);
title('Lena-y after fft and ifft with only Imaginary part');

%Part-II: Phase set to 0
figure('Name','Ex-2.6: Phase of the 2DFT - Part II)');

% Real part taken and absolute value used in terms of setting to phase to 0 
lena_y_phase = ifft2(abs(lena_y_fft));
lena_y_phase = real(lena_y_phase);

% Scaling with log, subtracting min and combined: Since log2 = 0, add 1 in
% any case before log
lena_y_log = log2(lena_y_phase + 1);
lena_y_min = lena_y_phase - min(lena_y_phase(:));
lena_y_logmin = log2(lena_y_phase - min(lena_y_phase(:)) + 1);

% Divide by its absolute value: Normalize/Scale image such that Magnitude = 1 
lena_y_mag = ifft2(lena_y_fft./abs(lena_y_fft));
lena_y_mag= real(lena_y_mag);

subplot(1,3,1);imshow(lena_y_min,[]);
title('Phase set to 0 Lena-y - Min');
subplot(1,3,2);imshow(lena_y_log,[]);
title('Phase set to 0 Lena-y - Log');
subplot(1,3,3);imshow(lena_y_logmin,[]);
title('Phase set to 0 Lena-y - Log2(x - min(x))');


figure('Name','Ex-2.6: Magnitude of the 2DFT - Part II)');
subplot(1,2,1);imshow(uint8(lena_y_gr));
title('Original Lena-y');
subplot(1,2,2);imshow(lena_y_mag,[]);
title('Magnitude set to 1 Lena-y');

%% Exercise 2.7

% VARIOUS SUBPLOTS TO DETECT ALPHA: READ COMMENTS FOR REQUIRED PLOTS(NOT ALL OF THEM ARE NECESSARY)
% UNNECESSARY ONES ARE COMMENTED, PLEASE UNCOMMENT TO SEE OTHER RESULTS
Lb = [10 200];
% figure('Name','Ex-2.7: Weber Law');
% tmp = 1;
% for i = 70 : 10 : 120
%     web = weber(10,i,Lb(1)); 
%     subplot(1,6,tmp);imshow(uint8(web));
%     title(['Lb = ' ,num2str(10), ', alpha = ', num2str((i - 10)/10)]);
%     tmp = tmp + 1;
% end
% figure;
% tmp2 = 1;
% for i = 70 : 10 : 120
%     web = weber(10,i,Lb(2)); 
%     subplot(1,6,tmp2);imshow(uint8(web));
%     title(['Lb = ' ,num2str(200), ', alpha = ', num2str((i - 10)/10)]);
%     tmp2 = tmp2 + 1;
% end

figure('Name','Ex-2.7: Weber Law- I');
tmp3 = 1;
% STEP SIZE CHANGED TO 4 INSTEAD OF 10: DETECTED CHANGE IN HERE(BLACK BACKGROUND) 
for i = 10 : 4 : 30
    web = weber(10,i,Lb(1)); 
    subplot(1,6,tmp3);imshow(uint8(web));
    title(['Lb = ' ,num2str(10), ', alpha = ', num2str((i - 10)/10)]);
    tmp3 = tmp3 + 1;
end
figure('Name','Ex-2.7: Weber Law- II');
tmp4 = 1;
% STEP SIZE CHANGED TO 4 INSTEAD OF 10: DETECTED CHANGE IN HERE(WHITE BACKGROUND) 
for i = 10 : 4 : 30
    web = weber(10,i,Lb(2)); 
    subplot(1,6,tmp4);imshow(uint8(web));
    title(['Lb = ' ,num2str(200), ', alpha = ', num2str((i - 10)/10)]);
    tmp4 = tmp4 + 1;
end
%figure;
%tmp5 = 1;
%for i = 130 : 10 : 180
%    web = weber(10,i,Lb(1)); 
%    subplot(1,6,tmp5);imshow(uint8(web));
%    title(['Lb = ' ,num2str(10), ', alpha = ', num2str((i - 10)/10)]);
%    tmp5 = tmp5 + 1;
%end
%figure;
%tmp6 = 1;
%for i = 130 : 10 : 180
%    web = weber(10,i,Lb(2)); 
%    subplot(1,6,tmp6);imshow(uint8(web));
%    title(['Lb = ' ,num2str(200), ', alpha = ', num2str((i - 10)/10)]);
%    tmp6 = tmp6 + 1;
%end