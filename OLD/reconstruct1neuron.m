% script for reconstructing a single neuron
% dependencies: denoise_bm4d, singleNeuronReconstruction
% must be passed imgFilePath
% must be passed threshold (0-1)
% must be passed sigma

%% parameters
if ~exist('imgFilePath','var')
    if usejava('desktop')
        imgFilePath = input('\nPlease enter the full path to the TIFF file:\n');
    else
        disp('No imgFilePath provided. Quitting MATLAB...');
        exit
    end
end

if ~exist('threshold','var')
    if usejava('desktop')
        threshold = input('\nPlease enter a threshold value (0-1):\n');
    else
        disp('No threshold provided. Arbitrarily choosing setting threshold = 0.1');
        threshold = .1;
    end
end

if ~exist('sigma','var')
    if usejava('desktop')
        threshold = input('\nPlease enter a value for sigma:\n');
    else
        disp('No sigma provided. Arbitrarily setting sigma = 500');
        sigma = 500;
    end
end

CHANNELCOUNT               = 3;


%% housekeeping
[imgDirpath,imgFileName,~] = fileparts(imgFilePath);
imgDirpath                 = [imgDirpath,'/',imgFileName,'_single_neuron_segmentation/'];
mkdir(imgDirpath);
paramDir = sprintf('%s/%s_sigma%d_threshold%0.4f/',imgDirpath,imgFileName,sigma,threshold);mkdir(paramDir);

info        = imfinfo(imgFilePath);
numXpixels  = info(1).Width;
numYpixels  = info(1).Height;
stackHeight = numel(info);

%% read in TIFF stack from imgFilePath
disp('Reading in TIFF stack')
try
    load(sprintf('%s%s_TIFFstack.mat',imgDirpath,imgFileName));
catch
    TIFFstack = zeros(numYpixels,numXpixels,stackHeight,CHANNELCOUNT);
    if ndims(imread(imgFilePath,'Index',1)) > 2 %#ok<ISMAT>
        parfor i = 1:stackHeight;
            TIFFstack(:,:,i,:) = imread(imgFilePath,'Index',i);
        end
    else
        for i = i:CHANNELCOUNT:stackHeight;
            TIFFstack(:,:,i,1) = imread(imgFilePath,'Index',i);
            TIFFstack(:,:,i,2) = imread(imgFilePath,'Index',i+1);
            TIFFstack(:,:,i,3) = imread(imgFilePath,'Index',i+2);
        end
    end
    save(sprintf('%s%s_TIFFstack.mat',imgDirpath,imgFileName),'TIFFstack','-v7.3');
end
disp('...done.')

%% denoise TIFF stack
disp('Denoising TIFF stack')
try 
    load(sprintf('%s/%s_denoised_sigma%d.mat',imgDirpath,imgFileName,sigma))
catch
    [bbVol] = denoise_BM4D(TIFFstack,sigma);
    save(sprintf('%s%s_denoised_sigma%d.mat',paramDir,imgFileName,sigma),'bbVol','-v7.3');
end
disp('...done.')

%% reconstruct single neuron
disp('Reconstructing single neuron')
[tmp, tmpColor] = singleNeuronReconstruction(bbVol, threshold);
disp('...done.')

%% write out max projection Z
disp('Writing out max Z projection')
imwrite(max(tmp,[],3), sprintf('%s%s_sigma%d_threshold%0.4f_max_Z_project.jpg',paramDir,imgFileName,sigma,threshold));
disp('...done.')

%% write out TIFF stack (cc)
disp('Writing out TIFF stack')
imwrite(double(tmp(:,:,1)), sprintf('%s/%s_sigma%d_threshold%0.4f.tiff',paramDir,imgFileName,sigma,threshold)); 
for kk = 2:size(tmp,3)
    imwrite(double(tmp(:,:,kk)), sprintf('%s%s_sigma%d_threshold%0.4f.tiff',paramDir,imgFileName,sigma,threshold), 'WriteMode', 'append'); 
end
disp('...done.')

%% write out JPGs (colored cc)
disp('Writing out colored JPGs')
mkdir(sprintf('%s/slices_sigma%d_threshold%0.4f/',paramDir,sigma,threshold));
for kk=1:size(tmpColor,3)
    imwrite(squeeze(tmpColor(:,:,kk,:)),sprintf('%s/slices_sigma%d_threshold%0.4f/%s_colored_sigma%d_threshold%0.4f_slice%d.jpg',paramDir,sigma,threshold,imgFileName,sigma,threshold,kk));
end
disp('...done.')
