%Analysis of the Reconstructed Image
function SNR=SNR_evaluation(image)
%{
    This function is used to find out the Signal to Noise ration of the
reconstructed image.
Input: recon_image  This is the reconstructed image whose SNR is to be estimated
Output: SNR  The evaluated value of SNR will be returned

Code Developed by: Hammad Omer, Department of Bioengineering, Imperial College London

Dated: 09/07/2009
%}

recon_image=image;
figure,colormap(gray(128));imagesc(abs(recon_image)); %Display the image so that the areas representing noise and signal can be selected
title('Signal to Noise Ratio Calculation');

%------------------Signal Calculation-----------------------
disp('Please select the region of interest to calculate the mean signal intensity');
[x,y]=ginput(2); %'x' is an array containing the x coordinates of both the points selected and 'y' is an array containing the y coordinates of both points selected
recon_image(round(y(1):y(2)),round(x(1)))=max(recon_image(:)); %The next five lines of code are just to highlight the area of selection in the image
recon_image(round(y(1):y(2)),round(x(2)))=max(recon_image(:));
recon_image(round(y(1)),round(x(1):x(2)))=max(recon_image(:));
recon_image(round(y(2)),round(x(1):x(2)))=max(recon_image(:));
colormap(gray(128));imagesc(abs(recon_image));

signal_region=image(round(y(1):y(2)),round(x(1)):x(2)); %Selection of the signal region to calculate the 'mean signal value'

%------------------Noise Calculation-------------------------

disp('Please select the background area representing Noise');
[x,y]=ginput(2); %'x' is an array containing the x coordinates of both the points selected and 'y' is an array containing the y coordinates of both points selected
recon_image(round(y(1):y(2)),round(x(1)))=max(recon_image(:));%The next five lines of code are just to highlight the area of selection in the image
recon_image(round(y(1):y(2)),round(x(2)))=max(recon_image(:));
recon_image(round(y(1)),round(x(1):x(2)))=max(recon_image(:));
recon_image(round(y(2)),round(x(1):x(2)))=max(recon_image(:));
colormap(gray(128));imagesc(abs(recon_image));

noise_region=image(round(y(1):y(2)),round(x(1)):x(2)); %Selection of the noise region to calculate the 'standard deviation of noise'

signal=mean2(abs(signal_region));
noise=std2(abs(noise_region));

SNR=abs(10*log10((signal/noise).^2)); %Calculation of Signal to Noise value
return