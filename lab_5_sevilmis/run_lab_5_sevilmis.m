% Lab 5 - Baris Sevilmis - 27/11/2019

%% EXERCISE 1 - LOAD DATA
%clear all; close all; clc;
load('data.mat');
load('real_subj_scores.mat');

folder = 'projectedViews/';
folder_models = 'models/';

filetype= fullfile(folder, '*.png');
filetype_models= fullfile(folder_models, '*.pcd');

images = dir(filetype);
models = dir(filetype_models);

% INITIALIZE 
model_len = length(models) - 5;
img_len = length(images) - 24;
img = cell(img_len,1);
img_ref = cell(24,1);
imgY = cell(img_len,1);
imgY_ref = cell(24,1);
model = cell(model_len,1);
model_ref = cell(5,1);

% READ REFERENCE IMAGES AND DISTORTED IMAGES(PROJECTED VIEWS)
tmp = 1;
for i = 1 : length(images)
    filenames = images(i).name;
    folder_filenames= fullfile(folder, filenames);
    curr = imread(folder_filenames);
    R = curr(:,:,1);
    G = curr(:,:,2);
    B = curr(:,:,3);
    [Y,~,~] = rgb2yuv(R,G,B);
    if mod(i,102) == 1 || mod(i,102) == 2 || mod(i,102) == 3 || mod(i,102) == 4 || mod(i,102) == 5 || mod(i,102) == 6
        %disp(folder_filenames);
        img_ref{tmp} = curr;
        imgY_ref{tmp} = Y;
        tmp = tmp + 1;
    else
        imgY{i - tmp + 1} = Y;
        img{i - tmp + 1} = curr;
    end
end

% READ REFERENCE POINTCLOUD MODELS & DISTORTED POINTCLOUD MODELS
tmp = 1;
for i = 1 : length(models)
    filenames_models = models(i).name;
    folder_filenames_models = fullfile(folder_models, filenames_models);
    if mod(i,17) == 1
        model_ref{tmp} = pcread(folder_filenames_models);
        tmp = tmp + 1;
    else
        model{i - tmp + 1} = pcread(folder_filenames_models);
        %pcshow(model{j}) to display image
    end
end
model(65:end) = [];
model_ref(5:end) = [];

% HELPERS INITILAZITION
cdc = zeros(4,16);
mdls = zeros(4,16);
dgrd = zeros(4,16);

% CREATE LUT TABLES - WILL HELP LATER

cdc(1,:) = find(strcmp(data.codec, 'octree-predlift'));
cdc(2,:) = find(strcmp(data.codec, 'octree-raht'));
cdc(3,:) = find(strcmp(data.codec, 'trisoup-predlift'));
cdc(4,:) = find(strcmp(data.codec, 'vpcc'));

mdls(1,:) = find(strcmp(data.model, 'amphoriskos'));
mdls(2,:) = find(strcmp(data.model, 'biplane'));
mdls(3,:) = find(strcmp(data.model, 'longdress'));
mdls(4,:) = find(strcmp(data.model, 'loot'));

data_mdls = zeros(size(data.model));
data_mdls(strcmp(data.model, 'amphoriskos')) = mdls(1,:); 
data_mdls(strcmp(data.model, 'biplane')) = mdls(2,:); 
data_mdls(strcmp(data.model, 'longdress')) = mdls(3,:); 
data_mdls(strcmp(data.model, 'loot')) = mdls(4,:);

dgrd(1,:) = find(strcmp(data.degradation, 'D1'));
dgrd(2,:) = find(strcmp(data.degradation, 'D2'));
dgrd(3,:) = find(strcmp(data.degradation, 'D3'));
dgrd(4,:) = find(strcmp(data.degradation, 'D4'));

codecs = unique(data.codec);
contents = unique(data.model);
pnts = data.points;
brs = data.bitrate;
%% EXERCISE 1.5 - DISPLAY A SPECIFIC MODEL

rand_num = randi([1 64], 1);
ex_PlotMODELS_sevilmis(model{rand_num});

%% EXERCISE 2 - OUTLIER REMOVAL

[final_subj_scores, r1, r2] = ex_Outlier_sevilmis(real_subj_scores);

%% EXERCISE 3 - MOS & CI COMPUTATION

[mos, ci] = ex_MOS_sevilmis(final_subj_scores);

%% EXERCISE 4 - PLOT MOS & CI
for i = 1 : size(contents,1)
    ex_PlotMOS_sevilmis(mos,ci,brs,mdls(i,:),cdc,contents,codecs,i);
end

%% EXERCISE 4.5 - FIND CLOSEST POINT OF MODEL AND ITS CORRESPONDING REFERENCES

knn_ref = cell(size(model,1),1);
knn_ref2 = cell(size(model,1),1);
color_ref = cell(size(model,1),1);
color_ref2 = cell(size(model,1),1);
tmp = 1;

% knnsearch(X,Y): Find closest points in X for each point in Y
% FIND KNN FOR BOTH DISTORTED TO REFERENCE AND REFERENCE TO DISTORTED:
% SYMMETRY
% BOTH LOCATIONS AND COLOR CHANNELS
for i = 1 : size(model,1)
    IDX = knnsearch(model_ref{tmp}.Location, model{i}.Location);
    IDX2 = knnsearch(model{i}.Location, model_ref{tmp}.Location);
    
    knn_ref{i} = model_ref{tmp}.Location(IDX,:);
    color_ref{i} = model_ref{tmp}.Color(IDX,:);
    knn_ref2{i} = model{i}.Location(IDX2,:); 
    color_ref2{i} = model{i}.Color(IDX2,:);
    
    if mod(i,16) == 0
        tmp = tmp + 1;
    end
end
%% EXERCISE 5 - SYMMETRIC POINT TO POINT METRIC: MSE & HAUSDORFF
MSEPo2Po = zeros(size(model,1),1);
HAUPo2Po = zeros(size(model,1),1);
MSEPo2Pl = obj_scores.po2plane_MSE;
HAUPo2Pl = obj_scores.po2plane_HAU;

% SYMMETRIC MSE AND HAU SCORES: MAX OF THEM
% POINT TO POINT AND POINT TO PLANE
tmp = 1;
for i = 1 : size(model,1)
    [p2p_1m, p2p_1h] = ex_PointToPoint_sevilmis(model{i}.Location,knn_ref{i});
    [p2p_2m, p2p_2h] = ex_PointToPoint_sevilmis(knn_ref2{i},model_ref{tmp}.Location);
    MSEPo2Po(i)= max(p2p_1m, p2p_2m);
    HAUPo2Po(i)= max(p2p_1h, p2p_2h);
    if mod(i,16) == 0
        tmp = tmp + 1;
    end
end


%% EXERCISE 6 - NORMAL CALCULATION OF THE MODELS 

MSEPl2Pl = zeros(size(model,1),1);

% K specified as 6
% FIND NORMALS
for i = 1 : size(model,1)
    if i <= size(model_ref,1)
        model_ref{i}.Normal = pcnormals(model_ref{i}, 6);
    end
    model{i}.Normal = pcnormals(model{i}, 6);
end

% FIND PLANE TO PLANE DISTANCES: THIRD PARAM OF ANGULARSYMMETRY
% MSE USED
tmp = 1;
for i = 1 : size(model,1)
    [~, ~, MSEPl2Pl(i)] = angularSimilarity(model{i}, model_ref{tmp}, 'MSE');
    if mod(i,16) == 0
        tmp = tmp + 1;
    end
end

%% EXERCISE 7 - PSNR FOR RGB & YUV

PSNR_RGB_log = zeros(size(model,1),1);
PSNR_RGB_avg = zeros(size(model,1),1);
PSNR_YUV = zeros(size(model,1),1);

% UTILIZE FORMULA 12, 13 & 14: FIND COLOR BASED SYMMETRIC PSNR SCORES
% RGB AND YUV
tmp = 1;
for i = 1 : size(model,1)
    psnr_log_1 = ex_PSNR_sevilmis(model{i}.Color, color_ref{i});
    psnr_log_2 = ex_PSNR_sevilmis(color_ref2{i}, model_ref{tmp}.Color);
    PSNR_RGB_log(i) = min(psnr_log_1, psnr_log_2); 
    
    psnr_avg_1 = ex_PSNRAVG_sevilmis(model{i}.Color, color_ref{i});
    psnr_avg_2 = ex_PSNRAVG_sevilmis(color_ref2{i}, model_ref{tmp}.Color);
    PSNR_RGB_avg(i) = min(psnr_avg_1, psnr_avg_2);
    
    psnr_yuv_1 = ex_PSNRYUV_sevilmis(model{i}.Color, color_ref{i}); 
    psnr_yuv_2 = ex_PSNRYUV_sevilmis(color_ref2{i}, model_ref{tmp}.Color);
    PSNR_YUV(i) = min(psnr_yuv_1, psnr_yuv_2);
    
    if mod(i,16) == 0
        tmp = tmp + 1;
    end
end

%% EXERCISE 8 - PSNR, SSIM, MS-SSIM

peaksnr = zeros(64,6);
peaksnrY = zeros(64,6);
ssimY = zeros(64,6);
mssimY = zeros(64,6);

% FR METRICS AS IN LAB 4: LUMINANCE PLANES UTILIZED IN EACH OF THEM
% RGB & Y DISTINCTLY ONLY FOR PSNR 
K = [0.01 0.03];
window = fspecial('gaussian', 11, 1.5);
level = 5;
weight = [0.0448 0.2856 0.3001 0.2363 0.1333];
method = 'product';
tmp = 1;
indFR = 1;
for i = 1 : 6 : size(img,1)
    for j = 1 : 6
        % DO BOTH
        peaksnr(indFR,j) = psnr(img{i + j - 1}, img_ref{j + tmp - 1});
        peaksnrY(indFR,j) = psnr(imgY{i + j - 1}, imgY_ref{j + tmp - 1});

        % Y CHANNEL ONLY
        ssimY(indFR,j) = ssim(imgY{i + j - 1}, imgY_ref{j + tmp - 1});

        mssimY(indFR,j) = ssim_mscale_new(double(imgY{i + j - 1}), double(imgY_ref{j + tmp - 1}),...
                                     K, window, level, weight, method);
    end
    if (mod(i - 1, 96) == 0) && ((i - 1) ~= 0)
        tmp = tmp + 6;
    end
    indFR = indFR + 1;
end

%% EXERCISE 8.5 - WEIGHTING PROJECTIONS - COMPUTE NON-GRAY AMOUNT AND USE AS WEIGHT
weights = zeros(4,6);
weightsSum = zeros(4,1);
peaksnrf = zeros(64,1);
peaksnrYf = zeros(64,1);
ssimYf = zeros(64,1);
mssimYf = zeros(64,1);

% WEIGHTING OF MODELS ACCORDING TO 6 DIFFERENT REFERENCE MODELS
% COUNT AMOUNT OF NON-GRAY PIXELS(!= 192)
% MULTIPLY CORRESPONDING SCORES AND SUM: DIVIDE BY SUM OF WEIGHTS AT END
k = 1;
indW = 1;
for i = 1 : size(weights,1)
    for j = 1 : size(weights , 2)
        weights(i,j) = sum(sum(sum(img_ref{k} ~= 192)));
        k = k + 1;
    end
    weightsSum(i) =  sum(weights(i,:));
    peaksnrf(indW : indW + 15) = (peaksnr(indW : indW + 15, :) * weights(i,:)') / weightsSum(i);
    peaksnrYf(indW : indW + 15) = (peaksnrY(indW : indW + 15, :) * weights(i,:)') / weightsSum(i);
    ssimYf(indW : indW + 15) = (ssimY(indW : indW + 15, :) * weights(i,:)') / weightsSum(i);
    mssimYf(indW : indW + 15) = (mssimY(indW : indW + 15, :) * weights(i,:)') / weightsSum(i);
    indW = indW + 16;
end



%% EXERCISE 9 - PLOT OBJECTIVE SCORE VS BITRATES

ttmp = 1;
for i = 1 : size(contents,1)
    ex_PlotObj_sevilmis(MSEPo2Po, brs,mdls(i,:),cdc,contents,codecs,i, ttmp ,'Point-To-Point-MSE');
    ex_PlotObj_sevilmis(MSEPo2Pl, brs,mdls(i,:),cdc,contents,codecs,i, ttmp + 1,'Point-To-Plane-MSE');
    ex_PlotObj_sevilmis(MSEPl2Pl, brs,mdls(i,:),cdc,contents,codecs,i, ttmp + 2,'Plane-To-Plane-MSE');
    
    ex_PlotObj_sevilmis(HAUPo2Po, brs,mdls(i,:),cdc,contents,codecs,i, ttmp + 3 ,'Point-To-Point-HAU');
    ex_PlotObj_sevilmis(HAUPo2Pl, brs,mdls(i,:),cdc,contents,codecs,i, ttmp + 4,'Point-To-Plane-HAU');
    
    ex_PlotObj_sevilmis(PSNR_RGB_log, brs,mdls(i,:),cdc,contents,codecs,i, ttmp + 5,'PSNR-RGB-LOG');
    ex_PlotObj_sevilmis(PSNR_RGB_avg, brs,mdls(i,:),cdc,contents,codecs,i, ttmp + 6,'PSNR-RGB-AVG');
    ex_PlotObj_sevilmis(PSNR_YUV, brs,mdls(i,:),cdc,contents,codecs,i, ttmp + 7,'PSNR-YUV');
    
    ex_PlotObj_sevilmis(peaksnrf, brs,mdls(i,:),cdc,contents,codecs,i, ttmp + 8,'PSNR-RGB');
    ex_PlotObj_sevilmis(peaksnrYf, brs,mdls(i,:),cdc,contents,codecs,i, ttmp + 9,'PSNR-Y');
    ex_PlotObj_sevilmis(ssimYf, brs,mdls(i,:),cdc,contents,codecs,i,  ttmp + 10,'SSIM-Y');
    ex_PlotObj_sevilmis(mssimYf, brs,mdls(i,:),cdc,contents,codecs,i, ttmp + 11,'MS-SSIM-Y');

    ttmp = ttmp + 12;
end

%% EXERCISE 10 - SPEARMAN - PEARSON - RMSE - NO FITTING

% COMPUTE OVERALL AND CONTENT SPECIFIC PEARSON, SPEARMAN, RMSE
table_raw = zeros(size(contents,1),12);
table_raw2 = zeros(size(contents,1),12);
table_raw3 = zeros(size(contents,1),12);

table_raw_oa = zeros(1,12);
table_raw_oa2 = zeros(1,12);
table_raw_oa3 = zeros(1,12);

for n = 1 : size(contents,1)
    
    [table_raw(n,1), table_raw2(n,1), table_raw3(n,1)] = ex_SpearmanPearsonRMSE_sevilmis(mos,MSEPo2Po,mdls(n,:));
    [table_raw(n,2), table_raw2(n,2), table_raw3(n,2)] = ex_SpearmanPearsonRMSE_sevilmis(mos,MSEPo2Pl,mdls(n,:));  
    [table_raw(n,3), table_raw2(n,3), table_raw3(n,3)] = ex_SpearmanPearsonRMSE_sevilmis(mos,MSEPl2Pl,mdls(n,:));
    
    [table_raw(n,4), table_raw2(n,4), table_raw3(n,4)] = ex_SpearmanPearsonRMSE_sevilmis(mos,HAUPo2Po,mdls(n,:));
    [table_raw(n,5), table_raw2(n,5), table_raw3(n,5)] = ex_SpearmanPearsonRMSE_sevilmis(mos,HAUPo2Pl,mdls(n,:));  

    [table_raw(n,6), table_raw2(n,6), table_raw3(n,6)] = ex_SpearmanPearsonRMSE_sevilmis(mos,PSNR_RGB_log,mdls(n,:));
    [table_raw(n,7), table_raw2(n,7), table_raw3(n,7)] = ex_SpearmanPearsonRMSE_sevilmis(mos,PSNR_RGB_avg,mdls(n,:));  
    [table_raw(n,8), table_raw2(n,8), table_raw3(n,8)] = ex_SpearmanPearsonRMSE_sevilmis(mos,PSNR_YUV,mdls(n,:));

    [table_raw(n,9), table_raw2(n,9), table_raw3(n,9)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnrf,mdls(n,:));
    [table_raw(n,10), table_raw2(n,10), table_raw3(n,10)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnrYf,mdls(n,:));   
    [table_raw(n,11), table_raw2(n,11), table_raw3(n,11)] = ex_SpearmanPearsonRMSE_sevilmis(mos,ssimYf,mdls(n,:));
    [table_raw(n,12), table_raw2(n,12), table_raw3(n,12)] = ex_SpearmanPearsonRMSE_sevilmis(mos,mssimYf,mdls(n,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[table_raw_oa(1,1), table_raw_oa2(1,1), table_raw_oa3(1,1)] = ex_SpearmanPearsonRMSE_sevilmis(mos,MSEPo2Po,data_mdls);
[table_raw_oa(1,2), table_raw_oa2(1,2), table_raw_oa3(1,2)] = ex_SpearmanPearsonRMSE_sevilmis(mos,MSEPo2Pl,data_mdls);  
[table_raw_oa(1,3), table_raw_oa2(1,3), table_raw_oa3(1,3)] = ex_SpearmanPearsonRMSE_sevilmis(mos,MSEPl2Pl,data_mdls);

[table_raw_oa(1,4), table_raw_oa2(1,4), table_raw_oa3(1,4)] = ex_SpearmanPearsonRMSE_sevilmis(mos,HAUPo2Po,data_mdls);
[table_raw_oa(1,5), table_raw_oa2(1,5), table_raw_oa3(1,5)] = ex_SpearmanPearsonRMSE_sevilmis(mos,HAUPo2Pl,data_mdls);  

[table_raw_oa(1,6), table_raw_oa2(1,6), table_raw_oa3(1,6)] = ex_SpearmanPearsonRMSE_sevilmis(mos,PSNR_RGB_log,data_mdls);
[table_raw_oa(1,7), table_raw_oa2(1,7), table_raw_oa3(1,7)] = ex_SpearmanPearsonRMSE_sevilmis(mos,PSNR_RGB_avg,data_mdls);  
[table_raw_oa(1,8), table_raw_oa2(1,8), table_raw_oa3(1,8)] = ex_SpearmanPearsonRMSE_sevilmis(mos,PSNR_YUV,data_mdls);

[table_raw_oa(1,9), table_raw_oa2(1,9), table_raw_oa3(1,9)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnrf,data_mdls);
[table_raw_oa(1,10), table_raw_oa2(1,10), table_raw_oa3(1,10)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnrYf,data_mdls);   
[table_raw_oa(1,11), table_raw_oa2(1,11), table_raw_oa3(1,11)] = ex_SpearmanPearsonRMSE_sevilmis(mos,ssimYf,data_mdls);
[table_raw_oa(1,12), table_raw_oa2(1,12), table_raw_oa3(1,12)] = ex_SpearmanPearsonRMSE_sevilmis(mos,mssimYf,data_mdls);

%% EXERCISE 11 - MOS vs Objective Metrics

ttmp = 1;
ex_PlotMOSvsObj_sevilmis(mos, ci, MSEPo2Po,mdls,contents,ttmp,'Point-To-Point-MSE');
ex_PlotMOSvsObj_sevilmis(mos, ci, MSEPo2Pl,mdls,contents,ttmp+11,'Point-To-Plane-MSE');
ex_PlotMOSvsObj_sevilmis(mos, ci, MSEPl2Pl, mdls,contents,ttmp+22,'Plane-To-Plane-MSE');

ex_PlotMOSvsObj_sevilmis(mos, ci, HAUPo2Po,mdls,contents,ttmp+33,'Point-To-Point-HAU');
ex_PlotMOSvsObj_sevilmis(mos, ci, HAUPo2Pl,mdls,contents,ttmp+44,'Point-To-Plane-HAU');

ex_PlotMOSvsObj_sevilmis(mos, ci, PSNR_RGB_log,mdls,contents,ttmp+55,'PSNR-RGB-LOG');
ex_PlotMOSvsObj_sevilmis(mos, ci, PSNR_RGB_avg,mdls,contents,ttmp+66,'PSNR-RGB-AVG');
ex_PlotMOSvsObj_sevilmis(mos, ci, PSNR_YUV, mdls,contents,ttmp+77,'PSNR-YUV');

ex_PlotMOSvsObj_sevilmis(mos, ci, peaksnrf, mdls,contents,ttmp+88,'PSNR-RGB');
ex_PlotMOSvsObj_sevilmis(mos, ci, peaksnrYf, mdls,contents,ttmp+99,'PSNR-Y');
ex_PlotMOSvsObj_sevilmis(mos, ci, ssimYf, mdls,contents,ttmp+110,'SSIM-Y');
ex_PlotMOSvsObj_sevilmis(mos, ci, mssimYf, mdls,contents,ttmp+121,'MS-SSIM-Y');


%% EXERCISE 12 - LINEAR AND CUBIC FITTING

% MSEPo2Po
% LINEAR FITTING
p = polyfit(MSEPo2Po,mos,1); 
MSEPo2Po_lin = polyval(p,MSEPo2Po);

% CUBIC FITTING
p = polyfit(MSEPo2Po,mos,3);
MSEPo2Po_cub = polyval(p, MSEPo2Po);

% MSEPo2Pl
% LINEAR FITTING
p = polyfit(MSEPo2Pl,mos,1); 
MSEPo2Pl_lin = polyval(p,MSEPo2Pl);

% CUBIC FITTING
p = polyfit(MSEPo2Pl,mos,3);
MSEPo2Pl_cub = polyval(p, MSEPo2Pl);

% MSEPl2Pl
% LINEAR FITTING
p = polyfit(MSEPl2Pl,mos,1); 
MSEPl2Pl_lin = polyval(p,MSEPl2Pl);

% CUBIC FITTING
p = polyfit(MSEPl2Pl,mos,3);
MSEPl2Pl_cub = polyval(p, MSEPl2Pl);

% HAUPo2Po
% LINEAR FITTING
p = polyfit(HAUPo2Po,mos,1); 
HAUPo2Po_lin = polyval(p,HAUPo2Po);

% CUBIC FITTING
p = polyfit(HAUPo2Po,mos,3);
HAUPo2Po_cub = polyval(p, HAUPo2Po);

% HAUEPo2Pl
% LINEAR FITTING
p = polyfit(HAUPo2Pl,mos,1); 
HAUPo2Pl_lin = polyval(p,HAUPo2Pl);

% CUBIC FITTING
p = polyfit(HAUPo2Pl,mos,3);
HAUPo2Pl_cub = polyval(p, HAUPo2Pl);

% PSNR_RGB_log
% LINEAR FITTING
p = polyfit(PSNR_RGB_log,mos,1); 
PSNR_RGB_log_lin = polyval(p,PSNR_RGB_log);

% CUBIC FITTING
p = polyfit(PSNR_RGB_log,mos,3);
PSNR_RGB_log_cub = polyval(p, PSNR_RGB_log);

% PSNR_RGB_avg
% LINEAR FITTING
p = polyfit(PSNR_RGB_avg,mos,1); 
PSNR_RGB_avg_lin = polyval(p,PSNR_RGB_avg);

% CUBIC FITTING
p = polyfit(PSNR_RGB_avg,mos,3);
PSNR_RGB_avg_cub = polyval(p, PSNR_RGB_avg);

% PSNR_YUV
% LINEAR FITTING
p = polyfit(PSNR_YUV,mos,1); 
PSNR_YUV_lin = polyval(p,PSNR_YUV);

% CUBIC FITTING
p = polyfit(PSNR_YUV,mos,3);
PSNR_YUV_cub = polyval(p, PSNR_YUV);

% PSNR-RGB
% LINEAR FITTING
p = polyfit(peaksnrf,mos,1); 
peaksnr_lin = polyval(p,peaksnrf);

% CUBIC FITTING
p = polyfit(peaksnrf,mos,3);
peaksnr_cub = polyval(p, peaksnrf);

% PSNR-Y
% LINEAR FITTING
p = polyfit(peaksnrYf,mos,1); 
peaksnrY_lin = polyval(p,peaksnrYf);

% CUBIC FITTING
p = polyfit(peaksnrYf,mos,3);
peaksnrY_cub = polyval(p, peaksnrYf);

% ssimY Y
% LINEAR FITTING
p = polyfit(ssimYf,mos,3);
ssimY_lin = polyval(p, ssimYf);

% CUBIC FITTING
p = polyfit(ssimYf,mos,3);
ssimY_cub = polyval(p, ssimYf);

% MssimY Y
% LINEAR FITTING
p = polyfit(mssimYf,mos,1); 
mssimY_lin = polyval(p,mssimYf);

% CUBIC FITTING
p = polyfit(mssimYf,mos,3);
mssimY_cub = polyval(p, mssimYf);

%% EXERCISE 13 - PLOT MOS VS OBJECTIVE FOR LINEAR FITTING AND CUBIC FITTING

% LINEAR PLOTS
ttmp = 300;
ex_PlotMOSvsObj_sevilmis(mos, ci, MSEPo2Po_lin, mdls,contents,ttmp,'Point-To-Point(Lin)-MSE');
ex_PlotMOSvsObj_sevilmis(mos, ci, MSEPo2Pl_lin, mdls,contents,ttmp+1,'Point-To-Plane(Lin)-MSE');
ex_PlotMOSvsObj_sevilmis(mos, ci, MSEPl2Pl_lin, mdls,contents,ttmp+2,'Plane-To-Plane(Lin)-MSE');

ex_PlotMOSvsObj_sevilmis(mos, ci, HAUPo2Po_lin, mdls,contents,ttmp+3,'Point-To-Point(Lin)-HAU');
ex_PlotMOSvsObj_sevilmis(mos, ci, HAUPo2Pl_lin, mdls,contents,ttmp+4,'Point-To-Plane(Lin)-HAU');

ex_PlotMOSvsObj_sevilmis(mos, ci, PSNR_RGB_log_lin, mdls,contents,ttmp+5,'PSNR-RGB-LOG(Lin)');
ex_PlotMOSvsObj_sevilmis(mos, ci, PSNR_RGB_avg_lin, mdls,contents,ttmp+6,'PSNR-RGB-AVG(Lin)');
ex_PlotMOSvsObj_sevilmis(mos, ci, PSNR_YUV_lin, mdls,contents,ttmp+7,'PSNR-YUV(Lin)');

ex_PlotMOSvsObj_sevilmis(mos, ci, peaksnr_lin, mdls,contents,ttmp+8,'PSNR-RGB(Lin)');
ex_PlotMOSvsObj_sevilmis(mos, ci, peaksnrY_lin, mdls,contents,ttmp+9,'PSNR-Y(Lin)');
ex_PlotMOSvsObj_sevilmis(mos, ci, ssimY_lin, mdls,contents,ttmp+10,'SSIM-Y(Lin)');
ex_PlotMOSvsObj_sevilmis(mos, ci, mssimY_lin, mdls,contents,ttmp+11,'MS-SSIM-Y(Lin)');

% CUBIC PLOTS
ttmp = 500;

ex_PlotMOSvsObj_sevilmis(mos, ci, MSEPo2Po_cub, mdls,contents,ttmp,'Point-To-Point(Cubic)-MSE');
ex_PlotMOSvsObj_sevilmis(mos, ci, MSEPo2Pl_cub, mdls,contents,ttmp+1,'Point-To-Plane(Cubic)-MSE');
ex_PlotMOSvsObj_sevilmis(mos, ci, MSEPl2Pl_cub, mdls,contents,ttmp+2,'Plane-To-Plane(Cubic)-MSE');

ex_PlotMOSvsObj_sevilmis(mos, ci, HAUPo2Po_cub, mdls,contents,ttmp+3,'Point-To-Point(Cubic)-HAU');
ex_PlotMOSvsObj_sevilmis(mos, ci, HAUPo2Pl_cub, mdls,contents,ttmp+4,'Point-To-Plane(Cubic)-HAU');

ex_PlotMOSvsObj_sevilmis(mos, ci, PSNR_RGB_log_cub, mdls,contents,ttmp+5,'PSNR-RGB-LOG(Cubic)');
ex_PlotMOSvsObj_sevilmis(mos, ci, PSNR_RGB_avg_cub, mdls,contents,ttmp+6,'PSNR-RGB-AVG(Cubic)');
ex_PlotMOSvsObj_sevilmis(mos, ci, PSNR_YUV_cub, mdls,contents,ttmp+7,'PSNR-YUV(Cubic)');

ex_PlotMOSvsObj_sevilmis(mos, ci, peaksnr_cub, mdls,contents,ttmp+8,'PSNR-RGB(Cubic)');
ex_PlotMOSvsObj_sevilmis(mos, ci, peaksnrY_cub, mdls,contents,ttmp+9,'PSNR-RGB(Cubic)');
ex_PlotMOSvsObj_sevilmis(mos, ci, ssimY_cub, mdls,contents,ttmp+10,'SSIM-Y(Cubic)');
ex_PlotMOSvsObj_sevilmis(mos, ci, mssimY_cub, mdls,contents,ttmp+11,'MS-SSIM-Y(Cubic)');


%% EXERCISE 14 - SPEARMAN - PEARSON - RMSE - LINEAR AND CUBIC FITTING

% COMPUTE OVERALL AND CONTENT SPECIFIC PEARSON, SPEARMAN, RMSE
table = zeros(size(contents,1),12);
table2 = zeros(size(contents,1),12);
table3 = zeros(size(contents,1),12);

table4 = zeros(size(contents,1),12);
table5 = zeros(size(contents,1),12);
table6 = zeros(size(contents,1),12);

table_oa = zeros(1,12);
table_oa2 = zeros(1,12);
table_oa3 = zeros(1,12);

table_oa4 = zeros(1,12);
table_oa5 = zeros(1,12);
table_oa6 = zeros(1,12);

for n = 1 : size(contents,1)
    [table(n,1), table2(n,1), table3(n,1)] = ex_SpearmanPearsonRMSE_sevilmis(mos,MSEPo2Po_lin,mdls(n,:));
    [table(n,2), table2(n,2), table3(n,2)] = ex_SpearmanPearsonRMSE_sevilmis(mos,MSEPo2Pl_lin,mdls(n,:));  
    [table(n,3), table2(n,3), table3(n,3)] = ex_SpearmanPearsonRMSE_sevilmis(mos,MSEPl2Pl_lin,mdls(n,:));
    
    [table(n,4), table2(n,4), table3(n,4)] = ex_SpearmanPearsonRMSE_sevilmis(mos,HAUPo2Po_lin,mdls(n,:));
    [table(n,5), table2(n,5), table3(n,5)] = ex_SpearmanPearsonRMSE_sevilmis(mos,HAUPo2Pl_lin,mdls(n,:));  
    
    [table(n,6), table2(n,6), table3(n,6)] = ex_SpearmanPearsonRMSE_sevilmis(mos,PSNR_RGB_log_lin,mdls(n,:));
    [table(n,7), table2(n,7), table3(n,7)] = ex_SpearmanPearsonRMSE_sevilmis(mos,PSNR_RGB_avg_lin,mdls(n,:));  
    [table(n,8), table2(n,8), table3(n,8)] = ex_SpearmanPearsonRMSE_sevilmis(mos,PSNR_YUV_lin,mdls(n,:));

    [table(n,9), table2(n,9), table3(n,9)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnr_lin,mdls(n,:));
    [table(n,10), table2(n,10), table3(n,10)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnrY_lin,mdls(n,:));
    [table(n,11), table2(n,11), table3(n,11)] = ex_SpearmanPearsonRMSE_sevilmis(mos,ssimY_lin,mdls(n,:));
    [table(n,12), table2(n,12), table3(n,12)] = ex_SpearmanPearsonRMSE_sevilmis(mos,mssimY_lin,mdls(n,:));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [table4(n,1), table5(n,1), table6(n,1)] = ex_SpearmanPearsonRMSE_sevilmis(mos,MSEPo2Po_cub,mdls(n,:));
    [table4(n,2), table5(n,2), table6(n,2)] = ex_SpearmanPearsonRMSE_sevilmis(mos,MSEPo2Pl_cub,mdls(n,:));  
    [table4(n,3), table5(n,3), table6(n,3)] = ex_SpearmanPearsonRMSE_sevilmis(mos,MSEPl2Pl_cub,mdls(n,:));
    
    [table4(n,4), table5(n,4), table6(n,4)] = ex_SpearmanPearsonRMSE_sevilmis(mos,HAUPo2Po_cub,mdls(n,:));
    [table4(n,5), table5(n,5), table6(n,5)] = ex_SpearmanPearsonRMSE_sevilmis(mos,HAUPo2Pl_cub,mdls(n,:));  

    [table4(n,6), table5(n,6), table6(n,6)] = ex_SpearmanPearsonRMSE_sevilmis(mos,PSNR_RGB_log_cub,mdls(n,:));
    [table4(n,7), table5(n,7), table6(n,7)] = ex_SpearmanPearsonRMSE_sevilmis(mos,PSNR_RGB_avg_cub,mdls(n,:));  
    [table4(n,8), table5(n,8), table6(n,8)] = ex_SpearmanPearsonRMSE_sevilmis(mos,PSNR_YUV_cub,mdls(n,:));

    [table4(n,9), table5(n,9), table6(n,9)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnr_cub,mdls(n,:));
    [table4(n,10), table5(n,10), table6(n,10)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnrY_cub,mdls(n,:));
    [table4(n,11), table5(n,11), table6(n,11)] = ex_SpearmanPearsonRMSE_sevilmis(mos,ssimY_cub,mdls(n,:));
    [table4(n,12), table5(n,12), table6(n,12)] = ex_SpearmanPearsonRMSE_sevilmis(mos,mssimY_cub,mdls(n,:));
end

[table_oa(1,1), table_oa2(1,1), table_oa3(1,1)] = ex_SpearmanPearsonRMSE_sevilmis(mos,MSEPo2Po_lin,data_mdls);
[table_oa(1,2), table_oa2(1,2), table_oa3(1,2)] = ex_SpearmanPearsonRMSE_sevilmis(mos,MSEPo2Pl_lin,data_mdls);  
[table_oa(1,3), table_oa2(1,3), table_oa3(1,3)] = ex_SpearmanPearsonRMSE_sevilmis(mos,MSEPl2Pl_lin,data_mdls);

[table_oa(1,4), table_oa2(1,4), table_oa3(1,4)] = ex_SpearmanPearsonRMSE_sevilmis(mos,HAUPo2Po_lin,data_mdls);
[table_oa(1,5), table_oa2(1,5), table_oa3(1,5)] = ex_SpearmanPearsonRMSE_sevilmis(mos,HAUPo2Pl_lin,data_mdls);  

[table_oa(1,6), table_oa2(1,6), table_oa3(1,6)] = ex_SpearmanPearsonRMSE_sevilmis(mos,PSNR_RGB_log_lin,data_mdls);
[table_oa(1,7), table_oa2(1,7), table_oa3(1,7)] = ex_SpearmanPearsonRMSE_sevilmis(mos,PSNR_RGB_avg_lin,data_mdls);  
[table_oa(1,8), table_oa2(1,8), table_oa3(1,8)] = ex_SpearmanPearsonRMSE_sevilmis(mos,PSNR_YUV_lin,data_mdls);

[table_oa(1,9), table_oa2(1,9), table_oa3(1,9)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnr_lin,data_mdls);
[table_oa(1,10), table_oa2(1,10), table_oa3(1,10)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnrY_lin,data_mdls);   
[table_oa(1,11), table_oa2(1,11), table_oa3(1,11)] = ex_SpearmanPearsonRMSE_sevilmis(mos,ssimY_lin,data_mdls);
[table_oa(1,12), table_oa2(1,12), table_oa3(1,12)] = ex_SpearmanPearsonRMSE_sevilmis(mos,mssimY_lin,data_mdls);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[table_oa4(1,1), table_oa5(1,1), table_oa6(1,1)] = ex_SpearmanPearsonRMSE_sevilmis(mos,MSEPo2Po_cub,data_mdls);
[table_oa4(1,2), table_oa5(1,2), table_oa6(1,2)] = ex_SpearmanPearsonRMSE_sevilmis(mos,MSEPo2Pl_cub,data_mdls);  
[table_oa4(1,3), table_oa5(1,3), table_oa6(1,3)] = ex_SpearmanPearsonRMSE_sevilmis(mos,MSEPl2Pl_cub,data_mdls);

[table_oa4(1,4), table_oa5(1,4), table_oa6(1,4)] = ex_SpearmanPearsonRMSE_sevilmis(mos,HAUPo2Po_cub,data_mdls);
[table_oa4(1,5), table_oa5(1,5), table_oa6(1,5)] = ex_SpearmanPearsonRMSE_sevilmis(mos,HAUPo2Pl_cub,data_mdls);  

[table_oa4(1,6), table_oa5(1,6), table_oa6(1,6)] = ex_SpearmanPearsonRMSE_sevilmis(mos,PSNR_RGB_log_cub,data_mdls);
[table_oa4(1,7), table_oa5(1,7), table_oa6(1,7)] = ex_SpearmanPearsonRMSE_sevilmis(mos,PSNR_RGB_avg_cub,data_mdls);  
[table_oa4(1,8), table_oa5(1,8), table_oa6(1,8)] = ex_SpearmanPearsonRMSE_sevilmis(mos,PSNR_YUV_cub,data_mdls);

[table_oa4(1,9), table_oa5(1,9), table_oa6(1,9)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnr_cub,data_mdls);
[table_oa4(1,10), table_oa5(1,10), table_oa6(1,10)] = ex_SpearmanPearsonRMSE_sevilmis(mos,peaksnrY_cub,data_mdls);   
[table_oa4(1,11), table_oa5(1,11), table_oa6(1,11)] = ex_SpearmanPearsonRMSE_sevilmis(mos,ssimY_cub,data_mdls);
[table_oa4(1,12), table_oa5(1,12), table_oa6(1,12)] = ex_SpearmanPearsonRMSE_sevilmis(mos,mssimY_cub,data_mdls);