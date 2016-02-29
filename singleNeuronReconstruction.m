function singleNeuronReconstruction(fn, th)
%fn = 'denoised_boyden_Brainbow_2x4_1_w1Conf488_bm4d_sigma2000';
%th = 0.1;
load(['/home/uygar/bb/data/' fn '.mat']);
bbVol = bbVol - min(bbVol(:));
bbVol = bbVol / max(bbVol(:));
binarybbVol=bbVol>th;
CC = bwconncomp(binarybbVol, 26);
compsizes=cellfun(@numel,CC.PixelIdxList);
[sorted,order]=sort(compsizes,'descend');
tmp=false(size(bbVol)); tmp(CC.PixelIdxList{order(1)})=true;
while true
  neighbors = setdiff(find(imdilate(tmp, true(5,5,5))), find(tmp));
  toRemove = [];
  for kk = 2:numel(order)
    if ~isempty(intersect(neighbors, CC.PixelIdxList{order(kk)}))
      tmp(CC.PixelIdxList{order(kk)}) =true;
      toRemove(end+1) = order(kk);
    end
  end
  if isempty(toRemove)
    break;
  else
    CC.PixelIdxList(toRemove) = [];
    compsizes=cellfun(@numel,CC.PixelIdxList);
    [sorted,order]=sort(compsizes,'descend');
  end
end

fn = [fn num2str(th)];
cd ~/bb/results/
imwrite(max(tmp,[],3), ['cc_' fn '.jpg']);
imwrite(double(tmp(:,:,1)), ['cc_' fn '.tif']); for kk=2:size(tmp,3); imwrite(double(tmp(:,:,kk)), ['cc_' fn '.tif'], 'WriteMode', 'append'); end
CCtmp = bwconncomp(tmp, 26);  
cd ~/bb/distinguishable_colors/; myColors = distinguishable_colors(CCtmp.NumObjects, [0 0 0]); cd ~/bb; tmpColor = zeros(size(tmp, 1), size(tmp, 2), size(tmp, 3), 3); voxelCount = prod(size(tmp));
for kk=1:CCtmp.NumObjects; tmpColor(CCtmp.PixelIdxList{kk})=myColors(kk,1); tmpColor(CCtmp.PixelIdxList{kk}+voxelCount)=myColors(kk,2); tmpColor(CCtmp.PixelIdxList{kk}+2*voxelCount)=myColors(kk,3); end;
cd ~/bb/results/tmp; for kk=1:size(tmpColor,3); imwrite(squeeze(tmpColor(:,:,kk,:)),['colored_cc_' fn '_' num2str(kk) '.jpg']); end;

