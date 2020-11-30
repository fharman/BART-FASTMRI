h5info('file1000052.h5'); 
Struct1=h5read('file1000052.h5','/kspace'); 
Kspace_image=Struct1.r+i*Struct1.i; 
LAMBDA=0.075;
ITER=100;
LAGRANGIAN=10;
  
K=[0.01 0.03];
L=255 ;
sigma=1;
filtersize=7;
 
  
INPLANE_ACC=2;
AUTO_CAL=32;
 
% % imshow(fftshift(ifft2(ifftshift(Kspace_image(:,:,30))),[])) 
% W0=fftshift(ifft2(ifftshift(Kspace_image(:,:,30))));
W0=Kspace_image(:,:,30);
W1=Kspace_image(:,:,30:31);
[row_kspace,column_kspace,size1]=size(Kspace_image); 
W0_x=fftshift(ifft2(ifftshift(Kspace_image(:,:,30))));
W1_x=fftshift(ifft2(ifftshift(Kspace_image(:,:,30:31))));
SIZEx=row_kspace;
SIZEy=column_kspace;
  
tmp_poisson=bart(sprintf('poisson -Y %i  -y %i  -Z %i -z %i  -C %i -e',SIZEx,INPLANE_ACC,SIZEy ,INPLANE_ACC,AUTO_CAL ));
  
   
tmp_poisson1=squeeze(tmp_poisson) ;
  
W1_zf= bart('fmac' ,W1, tmp_poisson1); 
  
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

sens= bart('ecalib -m 1 ', und2x2(:,:,1));

reco = bart('pics -r 0.01',und2x2(:,:,1),sens);


sense_recon = squeeze(reco);
figure, imshow(abs(sense_recon), []); title('SENSE Reconstruction')


Reconst_image=abs(sense_recon);%%reconstructed image 

Ground_truth= abs(W0_x);%%%% reference image


[window] = window1(sigma,filtersize);


[MSE, NMSE, PSNR,CNR, SSIM,L1error] = metrics1(Reconst_image,Ground_truth, K, window,L)
