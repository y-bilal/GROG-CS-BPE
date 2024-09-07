clc; clear all; close all;

setpath_ESP;

%%  Parameters

% Radial
N = 256;
proj =400;
Nc = 8;
os =1;
af=8;
%BPE
NBP = 81;
k_max = 0.4;
rand_sys = 1; % make it pseudo random
p_fac = 2;
% Gridding
s = 2; %oversampling Factor
step_size = 0.1;  %step size

load ('1_5T_Human_head_data');
im_ref=abs(sos(temp)); colormap(gray(256));
%% %% %%%%%%%%%%%%% Changing the dimensions(N x N x Nc )to (Nc x N x N) %%%%%%%%%%%%%%%%%%%%%%%%%%
im_low = permute(temp,[3,1,2]);
%% % Make Radial trajectoy
[kx,ky,dcf]=make_radial_traj(proj,N,os);  

%% %  First we have to Degrid Cartesain data to non-Cartesian data using Fessler toolbox
G=make_G_object(kx,ky,N);
sig_rad=reshape(degrid_Fessler(im_low,G),Nc,os*N,proj);
% Make density compensation function for radial data
dcf=repmat(abs(kx(:,1)),1,proj);

%%  Undersampled radial data 
sig_rad_u=sig_rad(:,:,1:af:end);
kxu=kx(:,1:af:end);
kyu=ky(:,1:af:end);
dcfu=dcf(:,1:af:end);
%% %%%%%%%%%%%%%%%%%% GROG-BPE Generation%%%%%%%%%%%%%%%%%%%%%%
% Use GROG_BPE to gererated BPE signal with trajectory
tic;
[kspave1,kx_blip,ky_blip,dcf_blip]=GROG_BPE(sig_rad_u,kxu,kyu,dcfu,NBP,k_max, rand_sys, p_fac);
t_gen = toc
%% %%%%%%%%%%%%%%%%%% Gridding of GROG-BPE Data%%%%%%%%%%%%%%%%%%%%%%
% Use GROG to grid the BPE sinal on cartesian grid
%Use original radial signal to calculate the weight sets first
[wsy,wsx]=SelfCalibratingWeightsRadial(kxu,kyu,sig_rad_u);
tic
[imnew,kspave,wsx,wsy]=sc_grog_radial(kspave1,kx_blip,ky_blip,step_size,s,wsx,wsy);
t_grid=toc
%% %%%%%%%%%%%%% Changing the dimensions (8x515x515) to (515x515x8)%%%%%%%%%%%%%%%%%%%%%%%%%%
k_sps=permute(kspave,[2,3,1]);

%% %%%%%%%%%%%%% Trimming of Gridded Kspace%%%%%%%%%%%%%%%%%%%%%%%%%%
% For s = 2, after gridding, the size of k_sp should be 512x512x8 but it is
% 515x515x8 or even larger in case of higher k_max. To make it an exact sxN, we have to trim the extra points
diff_size = size(k_sps,1) - s*N;
trim = floor(diff_size/2);
if(s==1)
    k_sp_grid_trim = k_sps(trim+1:end-trim,trim+1:end-trim,:);
end% for s = 1
 if (s==2)
k_sp_grid_trim = k_sps(trim+2:end-trim,trim+1:end-trim-1,:);
 end% for s = 2
% load  BPE_81_af_8_gridded;
data_mc = k_sp_grid_trim;%BPE Kspace after gridding and trimming from 515x515x8 to 512x512x8
mask = sum(data_mc,3)~=0;
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATAc = data_mc;
clearvars -except  DATAc mask   af  s N im_ref proj im_us NBP k_max t_grid t_gen
%%

%%ESPIRIT parameters
ncalib = 30; 
ksize = [6,6]; % ESPIRiT kernel-window-size
eigThresh_k = 0.02; % threshold of eigenvectors in k-space
eigThresh_im = 0.9; % threshold of eigenvectors in image space
mask = repmat(mask,[1,1,8]);
%% parameters for L1-reconstruction with splitting
nIterCG = 5;       % number of CG iterations for the PI part
nIterSplit = 30;    % number of splitting iterations for CS part
splitWeight = 0.4;  % reasonable value
lambda = 0.000025;    % L1-Wavelet threshold

%%Extracting calibration data
[sx,sy,Nc] = size(DATAc);
calib = crop_4d(DATAc,[ncalib,ncalib,Nc]);

%%
% Display coil images: 
im = ifft2c(DATAc);


figure, imshow3(abs(im),[],[1,Nc]); 
title('magnitude of physical coil images');
colormap((gray(256))); colorbar;

figure, imshow3(angle(im),[],[1,Nc]); 
title('phase of physical coil images');
colormap('default'); colorbar;

%% Compute Eigen-Value Maps
% Maps are computed in two steps. 


% compute Calibration matrix, perform 1st SVD and convert singular vectors
% into k-space kernels

[k,S] = dat2Kernel(calib,ksize);

idx = max(find(S >= S(1)*eigThresh_k));



%%
% crop kernels and compute eigen-value decomposition in image space to get
% maps
[M,W] = kernelEig(k(:,:,:,1:idx),[sx,sy]);

%% Compute Soft-SENSE ESPIRiT Maps 
% crop sensitivity maps according to eigenvalues==1. Note that we have to
% use 2 sets of maps. Here we weight the 2 maps with the eigen-values

maps = M(:,:,:,end-1:end);

% Weight the eigenvectors with soft-senses eigen-values
weights = W(:,:,end-1:end) ;
weights = (weights - eigThresh_im)./(1-eigThresh_im).* (W(:,:,end-1:end) > eigThresh_im);
weights = -cos(pi*weights)/2 + 1/2;

% create and ESPIRiT operator
ESP = ESPIRiT(maps,weights);

%% Reconsturctions
% ESPIRiT CG reconstruction with soft-sense and 1 sets of maps
XOP = Wavelet('Daubechies_TI',4,6);
FT = p2DFT(mask,[sx,sy,Nc]);

disp('Performing ESPIRiT reconstruction from 2 maps')
tic; [reskESPIRiT, resESPIRiT] = cgESPIRiT(DATAc,ESP, nIterCG*3, 0.01,DATAc*0); toc

disp('Performing L1-ESPIRiT reconstruction from 2 maps')
tic
[resL1ESPIRiT] = cgL1ESPIRiT(DATAc, resESPIRiT*0, FT, ESP, nIterCG,XOP,lambda,splitWeight,nIterSplit);
T_RECON = toc

im_out = crop(abs(sos(resL1ESPIRiT)),256,256);
 figure(100), imshow(im_out,[]);
x_sig = [80;120];
y_sig = [80;120];
x_noise = [25;65];
y_noise = [25;65]; 
norm_ref = mat2gray(im_ref);
norm_recon=mat2gray(im_out);

SNR_out=modified_SNR_evaluation(norm_recon, x_sig, y_sig, x_noise, y_noise)
ap_out = artifact_power(norm_recon,norm_ref)
rmse_out = RMSE(norm_ref,norm_recon)


