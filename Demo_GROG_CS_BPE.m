%% Setting Up path for Fessler toolbox
%  1) Open irt folder
%  2) Open and run 'setup' file
%  3) you are good to go

close all;clear all;clc;
setpath_CS;
%%  Parameters

% Radial Data
N = 256;
proj =400;
Nc = 8;
os =1;
af=8;

%BPE
NBP = 81;
k_max = 0.4
rand_sys = 1; % 1 for Random Blipping 0 for systemetic blipping
p_fac = 2;

% Gridding
s = 2; %oversampling Factor
step_size = 0.1;  %step size

load 1_5T_Human_head_data;

im_ref=abs(sos(temp)); colormap(gray(256));
figure(1);
imshow(im_ref,[]);
title('Orignal Image in cartesian trajectory');
%%  Degridding Cartesain data to Radial data%% 
im_low = permute(temp,[3,1,2]); %Changing the dimensions(N x N x Nc )to (Nc x N x N)
[kx,ky,dcf]=make_radial_traj(proj,N,os);  % Make Radial trajectoy
G=make_G_object(kx,ky,N);
sig_rad=reshape(degrid_Fessler(im_low,G),Nc,os*N,proj);

%% Display Radial data (fully sampled)
% Grids radial data onto Cartesain grid using Fessler Gridding
im_rad1=grid_Fessler(sig_rad,dcf(:),G,N);
img_rad=permute(im_rad1,[2,3,1]);%Changing the dimensions (Nc x N x N)to(N x N x Nc )
im_rad= abs(sos(img_rad));
figure(), imshow(im_rad ,[]);
title(['Fully sampled Radial image']);

%%  Undersampling of radial data 
sig_rad_u=sig_rad(:,:,1:af:end);
kxu=kx(:,1:af:end);
kyu=ky(:,1:af:end);
dcfu=dcf(:,1:af:end);

%% Display the undersampled data

 % Grids radial data onto Cartesain grid for display
Gu=make_G_object(kxu,kyu,N);
im_rad1=grid_Fessler(sig_rad_u,dcfu(:),Gu,N);

im_rad_u1=permute(im_rad1,[2,3,1]);%Changing the dimensions (Nc x N x N)to(N x N x Nc )
im_rad_u = abs(sos(im_rad_u1));
figure(), imshow(im_rad_u ,[]);
title(['Undersampled Radial image after Fesslar gridding. AF = ',  num2str(af)]);

%% %%%%%%%%%%%%%%%%%% GROG-BPE Generation%%%%%%%%%%%%%%%%%%%%%%
% Use GROG_BPE to gererated BPE signal from the undersampled radial signal
[kspave1,kx_blip,ky_blip,dcf_blip]=GROG_BPE(sig_rad_u,kxu,kyu,dcfu,NBP,k_max, rand_sys, p_fac);

%% %%%%%%%%%%%%%%%%%% Gridding of GROG-BPE Data%%%%%%%%%%%%%%%%%%%%%%
% Use GROG to grid the BPE sinal on cartesian grid
[wsy,wsx]=SelfCalibratingWeightsRadial(kxu,kyu,sig_rad_u); 
[imnew,kspave,wsx,wsy]=sc_grog_radial(kspave1,kx_blip,ky_blip,step_size,s,wsx,wsy);

%% %%%%%%%%%%%%% Changing the dimensions (8x515x515) to (515x515x8)%%%%%%%%%%%%%%%%%%%%%%%%%%
k_sps=permute(kspave,[2,3,1]);

%% %%%%%%%%%%%%% Trimming of Gridded Kspace%%%%%%%%%%%%%%%%%%%%%%%%%%
% For s = 2, after gridding, the size of k_sp should be 512x512x8 but it is
% 515x515x8 or even larger in case of higher k_max. To make it an exact sxN, we have to trim the extra points
diff_size = size(k_sps,1) - s*N;
trim = floor(diff_size/2);
k_sp_grid_trim = k_sps(trim+2:end-trim,trim+1:end-trim-1,:);

data_mc = k_sp_grid_trim;   %BPE Kspace after gridding and trimming from 515x515x8 to 512x512x8
mask = sum(data_mc,3)~=0;

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except  data_mc mask af NBP k_max s N im_ref 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%
% CS Recon Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Nc = size(data_mc,3);
    TVWeight = 0.001
    xfmWeight =0.001;	

for(coil_index = 1:Nc)
    data = data_mc(:,:,coil_index);
    N1 = size(data); 	% image Size
    Itnlim = 15;		% Number of iterations
    pdf=1; 
    %generate Fourier sampling operator
    FT = p2DFT(mask,N1,1,2);

    % scale data
    im_dc = FT'*(data.*mask./pdf);
    data = data/max(abs(im_dc(:)));
    im_dc = im_dc/max(abs(im_dc(:)));

    %generate transform operator
    XFM = Wavelet('Daubechies',6,4);	% Wavelet
    param = init;
    param.FT = FT;
    param.XFM = XFM;
    param.TV = TVOP;
    param.data = data;
    param.TVWeight =TVWeight;     % TV penalty 
    param.xfmWeight = xfmWeight;  % wavelet penalty
    param.Itnlim = Itnlim;
    res = XFM*im_dc;

    % do iterations
    for n=1:5
        res = fnlCg(res,param);
        im_res = XFM'*res;
    end
    im_res_mc(:,:,coil_index) = im_res ;
end
toc
%% %%%%%%%%%%%%% Find Compsosite Images %%%%%%%%%%%%%%%%%%%%%%%%%%%%
im_res_comp = abs(sos(im_res_mc));
%% %%%%%%%%%%%%% Final Cropped images in case of Oversampled data%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im_out = crop(im_res_comp,N,N);

%% %%%%%%%%%%%%%%Normalized Images%%%%%%%%%%%%
norm_ref = mat2gray(im_ref);
norm_recon=mat2gray(im_out);

%% %%%%%%%%%%%%% Quantifying parameters of reconstructed image %%%%%%%%%%%%%%%%%%%%%%%%%%%%
AP_out = artifact_power(norm_recon,norm_ref)
RMSE_out = RMSE(norm_ref,norm_recon)
x_sig = [80;120];
y_sig = [80;120];
x_noise = [25;65];
y_noise = [25;65]; 
SNR_out=modified_SNR_evaluation(im_out, x_sig, y_sig, x_noise, y_noise)
figure();imshow(norm_recon,[]); title('Reconstructed Image');
