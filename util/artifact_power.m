function artifact_pwr=artifact_power(recon_image,ref_image)
%{
    This function is used to find out the Artifact Power of the
reconstructed image.
Input: recon_image  This is the reconstructed image whose SNR is to be estimated, ref_image:  This is the reference image 
Output: artifact_pwr    The artifact power e.g. '0' artifact power means there is no artifact present and so on

Code Developed by: Hammad Omer, Department of Bioengineering, Imperial College London

Dated: 09/07/2009
%}

recon_image=abs(recon_image);
ref_image=abs(ref_image);
rows=size(ref_image,1);
columns=size(ref_image,2);

if rows~=size(recon_image,1)||columns~=size(recon_image,2) %Check if the dimensions of recon_image and ref_image are not equal then fill zeros to equate the sizes
    recon_image=zerofilling(recon_image,rows,columns);
end
image=ref_image;
max_image=max(image(:)); %Find the maximum value of the reference image
image=image./max_image; %Normalize the reference image with respect to its maximum value

max_recon=max(recon_image(:)); %Find the maximum value of the reconstructed image
recon_image=recon_image./max_recon; %Normalize the reference image with respect to its maximum value

diff_image=image-recon_image; %Find the difference image between the two normalized images
sq_diff_image=diff_image.^2; %Compute the square of the difference image
sum_diff_image=sum(sq_diff_image(:)); %Take Sum of the values in the difference image

sq_ref_image=image.^2; %Compute the square of the normalized reference image
sum_ref_image=sum(sq_ref_image(:));%Take Sum of the values in the square reference image 

artifact_pwr=sum_diff_image./sum_ref_image; %Compute the artifact power