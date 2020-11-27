

function [MSE, NMSE, PSNR,CNR, SSIM,L1error] = metrics1(Reconst_image,Ground_truth, K, window,L)

% for the window , use window function--it gives a gaussian window
% according to filter size and sigma.
%  K=[0.01 0.03] is used generally. L=255 (2^bitrate-1).

Reconst_image=im2double(Reconst_image);%%reconstructed image 
Ground_truth=im2double(Ground_truth);%%%% reference image
row=size(Ground_truth,1);
col=size(Ground_truth,2);
MSE=sum(sum(abs(Reconst_image-Ground_truth).^2))/(col*row);

%%%Normalized MSE
NMSE=MSE/sqrt(sum(sum(abs(Ground_truth).^2)));

%%%%% PSNR
Ratio_SNR=255/MSE;
PSNR=10*log10(Ratio_SNR);

%%%%%CNR---- It is not used as a metric in FastMRI
mean_reconst=mean(mean(Reconst_image));
mean_ground=mean(mean(Ground_truth));
CNR_ratio=(mean_reconst-mean_ground)/std(Ground_truth(:));
CNR = 20*log10(CNR_ratio);

%%%%% Structural Similarity Index

SSIM=ssim(Reconst_image,Ground_truth,'radius',1);
% SSIM= ssim4(Reconst_image,Ground_truth, K, window, L);
%%%%L1 error

L1error=sum(sum(abs(Reconst_image-Ground_truth)));

%%%%

end
