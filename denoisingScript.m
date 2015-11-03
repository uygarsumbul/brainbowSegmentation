info    = imfinfo('bbVol.tif');
thisVol = zeros(info(1).Height, info(1).Width, numel(info));
for zz = 1:numel(info);
  thisVol(:,:,zz) = imread('traced.tif','Index',zz);
end
CHANNELCOUNT = 3;
bbVol = zeros(size(thisVol,1),size(thisVol,2),size(thisVol,3)/CHANNELCOUNT, CHANNELCOUNT); %%%% 4 CHANNELS
for kk=1:size(bbVol,3);
  for mm = 1:CHANNELCOUNT
    bbVol(:,:,kk,mm)=thisVol(:,:,CHANNELCOUNT*(kk-1)+mm);
  end
end;
clear thisVol;
cd ~/bb/BM4D_v3p2/; sigma = 2000; for kk=1:4; [tmp, ~] = bm4d(squeeze(bbVol(:,:,:,kk)), 'Gauss', sigma); bbVol(:,:,:,kk) = tmp; end;
cd ~/bb/data;
save -v7.3 dcai_brainbow3_bm4d_sigma2000.mat bbVol
