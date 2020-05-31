%Lab 4 - Baris Sevilmis - 13/11/2019

%% Exercise 0 - INITIALIZE REQUIRED DATA

% Data Import
clear all; close all; clc;
load('data.mat');
load('real_subj_scores.mat');

%Video Data
folder_ref = 'references/';
folder1='mp4/';
folder2='mp4_2/';

filetype_ref = fullfile(folder_ref, '*.mp4');
filetype1= fullfile(folder1, '*.mp4');
filetype2 = fullfile(folder2, '*.mp4');

videos_ref = dir(filetype_ref);
videos1 = dir(filetype1);
videos2 = dir(filetype2);
videos = [videos1; videos2];

videos_arr_ref = cell(length(videos_ref),1);
videos_arr = cell(length(videos),1);

videos_orr_ref = cell(length(videos_ref),1);
videos_orr = cell(length(videos),1);

%SI and TI initilization
SI_ref = zeros(size(videos_ref,1),1);
TI_ref = zeros(size(videos_ref,1),1);

% Sampling of videos
% Sample size in terms of video duration: First frame of every second
sample_size = 10;
sample_frame = zeros(sample_size,1);

%Used to store specific frames for each video
frames_ref = cell(sample_size,1);
orr_ref = cell(sample_size,1);

frames = cell(sample_size,1);
orr = cell(sample_size,1);

% Helpers for later
cdc = zeros(4,20);
cnts = zeros(5,16);
br = zeros(4,20);

for m = 1:5
    if m < 5
        cdc(m,:) = find(codec_lut == m);
        br(m,:) = find(bitrate_lut == m);
    end
    cnts(m,:) = find(content_lut == m);
end
%% Exercise 1 - SI & TI OF REFERENCE CONTENTS + READ AND SAVE REFERENCE VIDEOS
    
%For references
for v = 1 :size(videos_ref,1)
    
    %Process Videos 1 by 1
    filenames_ref = videos_ref(v).name;
    disp(filenames_ref);
    folder_filenames_ref = fullfile(folder_ref, filenames_ref);
    video = VideoReader(folder_filenames_ref);

    % RGB TO YUV after VideoReader & BT709_I
    % Use Y channel to compute SI & TI
    fr = 1;
    for i = 1:sample_size
        sample_frame(i) = fr;
        fr = fr + video.FrameRate;
    end
    
    curr_frame = 1;
    k = 1;
    while hasFrame(video)
        curr = readFrame(video);
        if ismember(curr_frame, sample_frame)           
            R = curr(:,:,1);
            G = curr(:,:,2);
            B = curr(:,:,3);
            [Y,~,~] = rgb2yuv(R,G,B);
            frames_ref{k} = Y;
            orr_ref{k} = curr;
            k = k + 1;
        end
        curr_frame = curr_frame + 1;
    end    
    %SI & TI + SAVE VIDEO FRAMES
    SI_ref(v) = ex_SI_sevilmis(video);
    TI_ref(v) = ex_TI_sevilmis(video);
    videos_arr_ref{v} = frames_ref;
    videos_orr_ref{v} = orr_ref;
end

%% Exercise 1.5 - READ AND SAVE PROCESSED VIDEOS
% IF PURE VIDEOS NEEDED OR NEED TO RECREATE VIDEO BASED DATA

% For processed videos
for v = 1 :size(videos,1)
    
    %Process Videos one by one
    filenames_f = filenames{v};
    disp(filenames_f);
    if exist(fullfile(folder1, filenames_f))
        folder_filenames_f = fullfile(folder1, filenames_f);
    else
        folder_filenames_f = fullfile(folder2, filenames_f);
    end
    
    fr = 1;
    for i = 1:sample_size
        sample_frame(i) = fr;
        fr = fr + video.FrameRate;
    end
    
    curr_frame = 1;
    k = 1;
    %SAVE BOTH 
    while hasFrame(video)        
        curr = readFrame(video);
        if ismember(curr_frame, sample_frame)
            
            R = curr(:,:,1);
            G = curr(:,:,2);
            B = curr(:,:,3);
            [Y,~,~] = rgb2yuv(R,G,B);
            frames{k} = Y;
            orr{k} = curr;
            k = k + 1;
        end
        curr_frame = curr_frame + 1;
    end
    %Video Data
    videos_arr{v} = frames;
    videos_orr{v} = orr;
end


%% REMOVE UNNECESSARY VARIABLES

clear frames;
clear frames_ref;
clear orr;
clear orr_ref;
clear videos1;
clear videos2;
clear videos;
clear videos_ref;

%% Exercise 2 - PLOT SI AND TI

ex_PlotSITI_sevilmis(SI_ref, TI_ref,1);

%% Exercise 3 - DETECT AND REMOVE OUTLIER

out = ex_Outlier_sevilmis(real_subj_scores);


%% Exercise 4 - COMPUTE MOS AND CI

[mos, ci] = ex_MOS_sevilmis(out);

%% Exercise 5 - PLOT MOS AND CI

for n = 1:5
    ex_PlotMOS_sevilmis(mos, ci, bitrate_lut, cnts(n,:), cdc, contents,codecs,n);
end

%% Exercise 6 - COMPUTE FR: PSNR, SSIM, MS-SSIM

peaksnr = zeros(size(content_lut,1), sample_size);
ssimval = zeros(size(content_lut,1), sample_size);
mssimval = zeros(size(content_lut,1), sample_size);

peaksnrRGB = zeros(size(content_lut,1), sample_size);
ssimvalRGB = zeros(size(content_lut,1), sample_size);
mssimvalRGB = zeros(size(content_lut,1), sample_size);

k = 1;
tmp = 1;

%MS-SSIM PARAMETERS
K = [0.01 0.03];
window = fspecial('gaussian', 11, 1.5);
level = 5;
weight = [0.0448 0.2856 0.3001 0.2363 0.1333];
method = 'product';
for i = 1:size(content_lut,1)
           
    for j = 1:sample_size
        % FOR Y
        peaksnr(i,j) = psnr(videos_arr{i}{j}, videos_arr_ref{k}{j});

        ssimval(i,j) = ssim(videos_arr{i}{j}, videos_arr_ref{k}{j});

        mssimval(i,j) = abs(ssim_mscale_new(videos_arr{i}{j}, videos_arr_ref{k}{j}, K, window, level, weight, method));

        % FOR RGB 
        peaksnrRGB(i,j) = psnr(videos_orr{i}{j}, videos_orr_ref{k}{j});

        ssimvalRGB(i,j) = ssim(videos_orr{i}{j}, videos_orr_ref{k}{j});

        mssimvalR = abs(ssim_mscale_new(videos_orr{i}{j}(:,:,1), videos_orr_ref{k}{j}(:,:,1), K, window, level, weight, method));
        mssimvalG= abs(ssim_mscale_new(videos_orr{i}{j}(:,:,2), videos_orr_ref{k}{j}(:,:,2), K, window, level, weight, method));
        mssimvalB = abs(ssim_mscale_new(videos_orr{i}{j}(:,:,3), videos_orr_ref{k}{j}(:,:,3), K, window, level, weight, method));
        mssimvalRGB(i,j) = (mssimvalR +mssimvalG + mssimvalB)/3;
    end
    
    % SWITCH REFERENCE VIDEO
    if mod(i, 16) == 0
        k = k + 1;
    end
end

peaksnr_score = mean(peaksnr,2);
peaksnrRGB_score = mean(peaksnrRGB,2);
ssim_score = mean(ssimval,2);
ssimRGB_score = mean(ssimvalRGB,2);
mssim_score = mean(abs(mssimval),2);
mssimRGB_score = mean(abs(mssimvalRGB),2);

%% Exercise 8 - PLOT FR
ttmp = 1;
for n = 1 : 5
    ex_PlotFR_sevilmis(peaksnr_score, bitrate_lut, cnts(n,:), cdc, contents,codecs, n, ttmp ,'PSNR');
    ex_PlotFR_sevilmis(peaksnrRGB_score, bitrate_lut, cnts(n,:), cdc, contents,codecs, n, ttmp + 1,'PSNRRGB');
    
    ex_PlotFR_sevilmis(ssim_score, bitrate_lut, cnts(n,:), cdc, contents,codecs,n,  ttmp + 2,'SSIM');
    ex_PlotFR_sevilmis(ssimRGB_score, bitrate_lut, cnts(n,:), cdc, contents,codecs,n, ttmp + 3,'SSIMRGB');
    
    ex_PlotFR_sevilmis(mssim_score, bitrate_lut, cnts(n,:), cdc,contents,codecs,n, ttmp + 4,'MSSIM');
    ex_PlotFR_sevilmis(mssimRGB_score, bitrate_lut, cnts(n,:), cdc,contents,codecs,n, ttmp + 5,'MSSIMRGB');
    ttmp = ttmp + 6;
end

%% Exercise 9 - COMPUTE NR: BRISQUE, NIQE, PIQE

brisque_scores = zeros(size(content_lut,1),sample_size);
niqe_scores = zeros(size(content_lut,1),sample_size);
piqe_scores = zeros(size(content_lut,1),sample_size);
tmp = 1;
for i = 1 : size(content_lut,1)
    for j = 1:sample_size
        brisque_scores(i,j) =  brisque(videos_orr{i}{j});
        niqe_scores(i,j) =  niqe(videos_orr{i}{j});
        piqe_scores(i,j) =  piqe(videos_orr{i}{j});
    end
end

brisque_score = mean(brisque_scores,2);
piqe_score = mean(piqe_scores,2);
niqe_score = mean(niqe_scores,2);

%% Exercise 10 - PLOT NR
ttmp = 1;
for n = 1 : 5
    ex_PlotNR_sevilmis(brisque_score, bitrate_lut, cnts(n,:), cdc, contents,codecs,n, ttmp,'Brisque');
    ex_PlotNR_sevilmis(niqe_score, bitrate_lut, cnts(n,:), cdc, contents,codecs,n, ttmp + 1,'Niqe');
    ex_PlotNR_sevilmis(piqe_score, bitrate_lut, cnts(n,:), cdc, contents,codecs,n, ttmp + 2,'Piqe');
    ttmp = ttmp + 3;
end

%% Exercise 11 - MOS vs OBJECTIVE SCORES FOR NO FITTING

ttmp = 1;
ex_PlotMOSvsObj_sevilmis(mos, ci, peaksnr_score,cnts,contents,ttmp,'Peaksnr');
ex_PlotMOSvsObj_sevilmis(mos, ci, peaksnrRGB_score,cnts,contents,ttmp+5,'PeaksnrRGB');
ex_PlotMOSvsObj_sevilmis(mos, ci, ssim_score, cnts,contents,ttmp+10,'SSIM');
ex_PlotMOSvsObj_sevilmis(mos, ci, ssimRGB_score, cnts,contents,ttmp+15,'SSIMRGB');
ex_PlotMOSvsObj_sevilmis(mos, ci, mssim_score, cnts,contents,ttmp+20,'MSSIM');
ex_PlotMOSvsObj_sevilmis(mos, ci, mssimRGB_score, cnts,contents,ttmp+25,'MSSIMRGB');
ex_PlotMOSvsObj_sevilmis(mos, ci, brisque_score, cnts,contents,ttmp+30,'Brisque');
ex_PlotMOSvsObj_sevilmis(mos, ci, niqe_score, cnts,contents,ttmp+35,'Niqe');
ex_PlotMOSvsObj_sevilmis(mos, ci, piqe_score, cnts,contents,ttmp+40,'Piqe');

%% Exercise 12 - SPEARMAN - PEARSON - RMSE - NO FITTING

% COMPUTE OVERALL AND CONTENT SPECIFIC PEARSON, SPEARMAN, RMSE
table_raw = zeros(5,9);
table_raw2 = zeros(5,9);
table_raw3 = zeros(5,9);
table_raw_oa = zeros(1,9);
table_raw_oa2 = zeros(1,9);
table_raw_oa3 = zeros(1,9);

for n = 1 : 5
    [table_raw(n,1), table_raw2(n,1), table_raw3(n,1)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnr_score,cnts(n,:));
    [table_raw(n,2), table_raw2(n,2), table_raw3(n,2)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnrRGB_score,cnts(n,:));
    
    [table_raw(n,3), table_raw2(n,3), table_raw3(n,3)] = ex_SpearmanPearsonRMSE_sevilmis(mos,ssim_score,cnts(n,:));
    [table_raw(n,4), table_raw2(n,4), table_raw3(n,4)] = ex_SpearmanPearsonRMSE_sevilmis(mos,ssimRGB_score,cnts(n,:));
    
    [table_raw(n,5), table_raw2(n,5), table_raw3(n,5)] = ex_SpearmanPearsonRMSE_sevilmis(mos,mssim_score,cnts(n,:));
    [table_raw(n,6), table_raw2(n,6), table_raw3(n,6)] = ex_SpearmanPearsonRMSE_sevilmis(mos,mssimRGB_score,cnts(n,:));
    
    [table_raw(n,7), table_raw2(n,7), table_raw3(n,7)] = ex_SpearmanPearsonRMSE_sevilmis(mos,brisque_score,cnts(n,:));
    
    [table_raw(n,8), table_raw2(n,8), table_raw3(n,8)] = ex_SpearmanPearsonRMSE_sevilmis(mos,niqe_score,cnts(n,:));
    
    [table_raw(n,9), table_raw2(n,9), table_raw3(n,9)] = ex_SpearmanPearsonRMSE_sevilmis(mos,piqe_score,cnts(n,:));
end

[table_raw_oa(1,1), table_raw_oa2(1,1), table_raw_oa3(1,1)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnr_score,content_lut(:));
[table_raw_oa(1,2), table_raw_oa2(1,2), table_raw_oa3(1,2)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnrRGB_score,content_lut(:));


[table_raw_oa(1,3), table_raw_oa2(1,3), table_raw_oa3(1,3)] = ex_SpearmanPearsonRMSE_sevilmis(mos,ssim_score,content_lut(:));
[table_raw_oa(1,4), table_raw_oa2(1,4), table_raw_oa3(1,4)] = ex_SpearmanPearsonRMSE_sevilmis(mos,ssimRGB_score,content_lut(:));

[table_raw_oa(1,5), table_raw_oa2(1,5), table_raw_oa3(1,5)] = ex_SpearmanPearsonRMSE_sevilmis(mos,mssim_score,content_lut(:));
[table_raw_oa(1,6), table_raw_oa2(1,6), table_raw_oa3(1,6)] = ex_SpearmanPearsonRMSE_sevilmis(mos,mssimRGB_score,content_lut(:));

[table_raw_oa(1,7), table_raw_oa2(1,7), table_raw_oa3(1,7)] = ex_SpearmanPearsonRMSE_sevilmis(mos,brisque_score,content_lut(:));

[table_raw_oa(1,8), table_raw_oa2(1,8), table_raw_oa3(1,8)] = ex_SpearmanPearsonRMSE_sevilmis(mos,niqe_score,content_lut(:));

[table_raw_oa(1,9), table_raw_oa2(1,9), table_raw_oa3(1,9)] = ex_SpearmanPearsonRMSE_sevilmis(mos,piqe_score,content_lut(:));


%% Exercise 13 - LINEAR AND CUBIC FITTING

% PSNR
% Linear Fitting
A = polyfit(peaksnr_score,mos,1); 
peaksnr_lin = polyval(A,peaksnr_score);

% Cubic Fitting
p = polyfit(peaksnr_score,mos,3);
peaksnr_cubic = polyval(p, peaksnr_score);

% PSNR RGB
% Linear Fitting
p = polyfit(peaksnrRGB_score,mos,1); 
peaksnrRGB_lin = polyval(p,peaksnrRGB_score);

% Cubic Fitting
p = polyfit(peaksnrRGB_score,mos,3);
peaksnrRGB_cubic = polyval(p, peaksnrRGB_score);

% SSIM
% Linear Fitting
p = polyfit(ssim_score,mos,1); 
ssim_lin = polyval(p,ssim_score);

%Cubic Fitting
p = polyfit(ssim_score,mos,3);
ssim_cubic = polyval(p, ssim_score);

% SSIM RGB
% Linear Fitting
p = polyfit(ssimRGB_score,mos,3);
ssimRGB_lin = polyval(p, ssimRGB_score);

%Cubic Fitting
p = polyfit(ssimRGB_score,mos,3);
ssimRGB_cubic = polyval(p, ssimRGB_score);

% MSSIM
% Linear Fitting
p = polyfit(mssim_score,mos,1); 
mssim_lin = polyval(p,mssim_score);

% Cubic Fitting
p = polyfit(mssim_score,mos,3);
mssim_cubic = polyval(p, mssim_score);

% MSSIM RGB
% Linear Fitting
p = polyfit(mssimRGB_score,mos,3);
mssimRGB_lin = polyval(p, mssimRGB_score);

% Cubic Fitting
p = polyfit(mssimRGB_score,mos,3);
mssimRGB_cubic = polyval(p, mssimRGB_score);

% BRISQUE
% Linear Fitting
p = polyfit(brisque_score,mos,1); 
brisque_lin = polyval(p,brisque_score);

% Cubic Fitting
p = polyfit(brisque_score,mos,3);
brisque_cubic = polyval(p, brisque_score);

% PIQE
% Linear Fitting
p = polyfit(piqe_score,mos,1); 
piqe_lin = polyval(p,piqe_score);

% Cubic Fitting
p = polyfit(piqe_score,mos,3);
piqe_cubic = polyval(p, piqe_score);

% NIQE
% Linear Fitting
p = polyfit(niqe_score,mos,1); 
niqe_lin = polyval(p,niqe_score);

% Cubic Fitting
p = polyfit(niqe_score,mos,3);
niqe_cubic = polyval(p, niqe_score);

%% Exercise 14 - PLOT MOS VS OBJECTIVE FOR LINEAR FITTING

% Linear Plots
ttmp = 200;

ex_PlotMOSvsObj_sevilmis(mos, ci, peaksnr_lin, cnts,contents,ttmp,'Peaksnr(Lin)');
ex_PlotMOSvsObj_sevilmis(mos, ci, peaksnrRGB_lin, cnts,contents,ttmp+1,'PeaksnrRGB(Lin)');
ex_PlotMOSvsObj_sevilmis(mos, ci, ssim_lin, cnts,contents,ttmp+2,'SSIM(Lin)');
ex_PlotMOSvsObj_sevilmis(mos, ci, ssimRGB_lin, cnts,contents,ttmp+3,'SSIMRGB(Lin)');
ex_PlotMOSvsObj_sevilmis(mos, ci, mssim_lin, cnts,contents,ttmp+4,'MSSIM(Lin)');
ex_PlotMOSvsObj_sevilmis(mos, ci, mssimRGB_lin, cnts,contents,ttmp+5,'MSSIMRGB(Lin)');
ex_PlotMOSvsObj_sevilmis(mos, ci, brisque_lin, cnts,contents,ttmp+6,'Brisque(Lin)');
ex_PlotMOSvsObj_sevilmis(mos, ci, niqe_lin, cnts,contents,ttmp+7,'Niqe(Lin)');
ex_PlotMOSvsObj_sevilmis(mos, ci, piqe_lin, cnts,contents,ttmp+8,'Piqe(Lin)');

% Cubic Plots
ttmp = 900;

ex_PlotMOSvsObj_sevilmis(mos, ci, peaksnr_cubic, cnts,contents,ttmp,'Peaksnr(Cubic)');
ex_PlotMOSvsObj_sevilmis(mos, ci, peaksnrRGB_cubic, cnts,contents,ttmp+1,'PeaksnrRGB(Cubic)');
ex_PlotMOSvsObj_sevilmis(mos, ci, ssim_cubic, cnts,contents,ttmp+2,'SSIM(Cubic)');
ex_PlotMOSvsObj_sevilmis(mos, ci, ssimRGB_cubic, cnts,contents,ttmp+3,'SSIMRGB(Cubic)');
ex_PlotMOSvsObj_sevilmis(mos, ci, mssim_cubic, cnts,contents,ttmp+4,'MSSIM(Cubic)');
ex_PlotMOSvsObj_sevilmis(mos, ci, mssimRGB_cubic, cnts,contents,ttmp+5,'MSSIMRGB(Cubic)');
ex_PlotMOSvsObj_sevilmis(mos, ci, brisque_cubic, cnts,contents,ttmp+6,'Brisque(Cubic)');
ex_PlotMOSvsObj_sevilmis(mos, ci, niqe_cubic, cnts,contents,ttmp+7,'Niqe(Cubic)');
ex_PlotMOSvsObj_sevilmis(mos, ci, piqe_cubic, cnts,contents,ttmp+8,'Piqe(Cubic)');


%% Exercise 15 - PEARSON - SPEARMAN- RMSE FOR LINEAR AND CUBIC FITTING

%Linear and Cubic together for each content
table = zeros(5,9);
table_2 = zeros(5,9);

table_3 = zeros(5,9);
table_4 = zeros(5,9);

table_5 = zeros(5,9);
table_6 = zeros(5,9);

table_oa = zeros(1,9);
table_oa3 = zeros(1,9);
table_oa5 = zeros(1,9);

table_oa2 = zeros(1,9);
table_oa4 = zeros(1,9);
table_oa6 = zeros(1,9);
for n = 1 : 5
    [table(n,1), table_3(n,1), table_5(n,1)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnr_lin,cnts(n,:));
    [table_2(n,1), table_4(n,1), table_6(n,1)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnr_cubic,cnts(n,:));
    [table(n,2), table_3(n,2), table_5(n,2)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnrRGB_lin,cnts(n,:));
    [table_2(n,2), table_4(n,2), table_6(n,2)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnrRGB_cubic,cnts(n,:));
    
    [table(n,3), table_3(n,3), table_5(n,3)] = ex_SpearmanPearsonRMSE_sevilmis(mos,ssim_lin,cnts(n,:));
    [table_2(n,3), table_4(n,3), table_6(n,3)] = ex_SpearmanPearsonRMSE_sevilmis(mos,ssim_cubic,cnts(n,:));
    [table(n,4), table_3(n,4), table_5(n,4)] = ex_SpearmanPearsonRMSE_sevilmis(mos,ssimRGB_lin,cnts(n,:));
    [table_2(n,4), table_4(n,4), table_6(n,4)] = ex_SpearmanPearsonRMSE_sevilmis(mos,ssimRGB_cubic,cnts(n,:));
    
    [table(n,5), table_3(n,5), table_5(n,5)] =ex_SpearmanPearsonRMSE_sevilmis(mos,mssim_lin,cnts(n,:));
    [table_2(n,5), table_4(n,5), table_6(n,5)] = ex_SpearmanPearsonRMSE_sevilmis(mos,mssim_cubic,cnts(n,:));
    [table(n,6), table_3(n,6), table_5(n,6)] = ex_SpearmanPearsonRMSE_sevilmis(mos,mssimRGB_lin,cnts(n,:));
    [table_2(n,6), table_4(n,6), table_6(n,6)] = ex_SpearmanPearsonRMSE_sevilmis(mos,mssimRGB_cubic,cnts(n,:));
    
    [table(n,7), table_3(n,7), table_5(n,7)] = ex_SpearmanPearsonRMSE_sevilmis(mos,brisque_lin,cnts(n,:));
    [table_2(n,7), table_4(n,7), table_6(n,7)] = ex_SpearmanPearsonRMSE_sevilmis(mos,brisque_cubic,cnts(n,:));
    
    [table(n,8), table_3(n,8), table_5(n,8)] = ex_SpearmanPearsonRMSE_sevilmis(mos,niqe_lin,cnts(n,:));
    [table_2(n,8), table_4(n,8), table_6(n,8)] = ex_SpearmanPearsonRMSE_sevilmis(mos,niqe_cubic,cnts(n,:));
    
    [table(n,9), table_3(n,9), table_5(n,9)] = ex_SpearmanPearsonRMSE_sevilmis(mos,piqe_lin,cnts(n,:));
    [table_2(n,9), table_4(n,9), table_6(n,9)] = ex_SpearmanPearsonRMSE_sevilmis(mos,piqe_cubic,cnts(n,:));
end

%OVERALL FOR LINEAR FITTING
[table_oa(1,1), table_oa3(1,1), table_oa5(1,1)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnr_lin,content_lut(:));
[table_oa(1,2), table_oa3(1,2), table_oa5(1,2)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnrRGB_lin,content_lut(:));

[table_oa(1,3), table_oa3(1,3), table_oa5(1,3)] = ex_SpearmanPearsonRMSE_sevilmis(mos,ssim_lin,content_lut(:));
[table_oa(1,4), table_oa3(1,4), table_oa5(1,4)] = ex_SpearmanPearsonRMSE_sevilmis(mos,ssimRGB_lin,content_lut(:));

[table_oa(1,5), table_oa3(1,5), table_oa5(1,5)] = ex_SpearmanPearsonRMSE_sevilmis(mos,mssim_lin,content_lut(:));
[table_oa(1,6), table_oa3(1,6), table_oa5(1,6)] = ex_SpearmanPearsonRMSE_sevilmis(mos,mssimRGB_lin,content_lut(:));

[table_oa(1,7), table_oa3(1,7), table_oa5(1,7)] = ex_SpearmanPearsonRMSE_sevilmis(mos,brisque_lin,content_lut(:));

[table_oa(1,8), table_oa3(1,8), table_oa5(1,8)] = ex_SpearmanPearsonRMSE_sevilmis(mos,niqe_lin,content_lut(:));

[table_oa(1,9), table_oa3(1,9), table_oa5(1,9)] = ex_SpearmanPearsonRMSE_sevilmis(mos,piqe_lin,content_lut(:));

%OVERALL FOR CUBIC FITTING
[table_oa2(1,1), table_oa4(1,1), table_oa6(1,1)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnr_cubic,content_lut(:));
[table_oa2(1,2), table_oa4(1,2), table_oa6(1,2)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnrRGB_cubic,content_lut(:));

[table_oa2(1,3), table_oa4(1,3), table_oa6(1,3)] = ex_SpearmanPearsonRMSE_sevilmis(mos,ssim_cubic,content_lut(:));
[table_oa2(1,4), table_oa4(1,4), table_oa6(1,4)] = ex_SpearmanPearsonRMSE_sevilmis(mos,ssimRGB_cubic,content_lut(:));

[table_oa2(1,5), table_oa4(1,5), table_oa6(1,5)] = ex_SpearmanPearsonRMSE_sevilmis(mos,mssim_cubic,content_lut(:));
[table_oa2(1,6), table_oa4(1,6), table_oa6(1,6)] = ex_SpearmanPearsonRMSE_sevilmis(mos,mssimRGB_cubic,content_lut(:));

[table_oa2(1,7), table_oa4(1,7), table_oa6(1,7)] = ex_SpearmanPearsonRMSE_sevilmis(mos,brisque_cubic,content_lut(:));

[table_oa2(1,8), table_oa4(1,8), table_oa6(1,8)] = ex_SpearmanPearsonRMSE_sevilmis(mos,niqe_cubic,content_lut(:));

[table_oa2(1,9), table_oa4(1,9), table_oa6(1,9)] = ex_SpearmanPearsonRMSE_sevilmis(mos,piqe_cubic,content_lut(:));
