h5info('file1000052.h5'); 
Struct1=h5read('file1000052.h5','/kspace'); 
Kspace_image=Struct1.r+i*Struct1.i; 
LAMBDA=0.075;
ITER=100;
 
K=[0.01 0.03];
L=255 ;
sigma=1;
filtersize=7;
 
% % imshow(fftshift(ifft2(ifftshift(Kspace_image(:,:,30))),[])) 
% W0=fftshift(ifft2(ifftshift(Kspace_image(:,:,30))));
W0=Kspace_image(:,:,30);
W1=Kspace_image(:,:,30:31);
[row_kspace,column_kspace,size1]=size(Kspace_image); 
W0_x=fftshift(ifft2(ifftshift(Kspace_image(:,:,30))));
W1_x=fftshift(ifft2(ifftshift(Kspace_image(:,:,30:31))));
% 
cartesian_mask=imread('Cartesian Mask'.png');
cartesian_mask=logical(cartesian_mask);
% cartesian_mask=imrotate(cartesian_mask,90);
cartesian_mask=double(imresize(cartesian_mask,[row_kspace,column_kspace]));
% 
cartesian_mask_kspace=(Kspace_image(:,:,30)).*(cartesian_mask);
size2=2;

for i=1:size2
  s(:,:,i)= cartesian_mask;    
end
W1_zf=W1.*s;
% sqrt sum-of-squares of k-space
und2x2 = W1_zf;
ksp_rss = bart('rss 012', und2x2);

% zero-filled reconstruction sqrt-sum-of-squares
zf_coils = bart('fft -i 3', W1_zf);

zf_rss = bart('rss 012', zf_coils);

ksp_rss = squeeze(ksp_rss);
zf_coils = squeeze(zf_coils);
zf_rss = squeeze(zf_rss);

figure,subplot(1,2,1), imshow(abs(ksp_rss).^0.125, []); title('k-space')
subplot(1,2,2), imshow(abs(zf_rss), []); title('zero-filled recon')

sens= bart('ecalib -m1', und2x2(:,:,1));
MASK = str2num(evalc("bart('bitmask 0 1')"));
reco = bart(sprintf('pics -R W:%i:0:%g -i %i',MASK,LAMBDA,ITER),und2x2(:,:,1),sens);

sense_recon = squeeze(reco);
figure, imshow(abs(sense_recon), []); title('ESPIRiT Reconstruction')

% nrmse_l1_wav=bart('nrmse',W0_x,sense_recon)

[window] = window1(sigma,filtersize);


[MSE, NMSE, PSNR,CNR, SSIM,L1error] = metrics1(Reconst_image,Ground_truth, K, window,L)
