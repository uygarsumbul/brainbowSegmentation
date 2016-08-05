function [denoised_img] = denoise_BM4D(TIFFstack4D,sigma)

% size(TIFFstack4D) = [numYpixels,numXpixels,numSlices,CHANNELCOUNT]
% requires bm4d

denoised_img = zeros(size(TIFFstack4D));

parfor i = 1:size(TIFFstack4D,4)
    denoised_img(:,:,:,i) = bm4d(TIFFstack4D(:,:,:,i), 'Gauss', sigma);
end

end
