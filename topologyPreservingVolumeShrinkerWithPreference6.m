function labels = topologyPreservingVolumeShrinkerWithPreference6(template, skeleton, erosionRadius, brightness)
% Software developed by: Uygar Sümbül <uygar@stat.columbia.edu, uygar@mit.edu> unless otherwise indicated in the function headers
% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE.
% IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DAMAGES WHATSOEVER.
%
% Grow the neuron on the trace skeleton by preserving the topology of the skeleton and the geometry of the neuron

% skeleton: 3d binary image stack where 1 represents skeleton voxels and 0 represents background voxels
% template: the template is supposed to shrink to the skeleton, as much as possible without changing topology
% dilationRadius: number of iterations the routine will run for.

% use 6-neighborhood to inflate the volume conservatively
sixNeigh = [5 11 13 14 15 17 23];
blur_filter=false(3,3,3); blur_filter(sixNeigh) = true; %blur_filter([14 sixNeigh(6)])=true;
%blur_filter(2,2,2)=1;blur_filter(2,2,1)=1;blur_filter(2,2,3)=1;
%blur_filter(2,1,2)=1;blur_filter(2,3,2)=1;blur_filter(1,2,2)=1;blur_filter(3,2,2)=1;

labels = template; oldLabels = labels; erosion = imerode(template,blur_filter); region = labels & ~erosion;
for kk = 1:erosionRadius
  % use 6-neighborhood in simple point testing to accept new voxels conservatively
  labels = simple_point_warp_6(labels, skeleton, region, brightness);
%  blur_filter=false(3,3,3); blur_filter([14 sixNeigh(rem(kk-1,6)+1)])=true;
  erosion = imerode(labels,blur_filter);
  region = labels & ~erosion;
  tmp = xor(labels, oldLabels);
  if ~any(tmp(:))
    break;
  end
  oldLabels = labels;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function source=simple_point_warp_6(source, target, mask, brightness)
% uses simple point relaxation to warp 3d source into 3d target. source is only modified at nonzero locations in the mask
% original authors: Viren Jain, Mark Richardson, M.I.T.
% modified by Uygar Sumbul

% make sure source and target are binary images
source=source>0.5; target=target>0.5;

diffCount = -Inf;
mask(1,:,:)=0; mask(:,:,1)=0; mask(:,1,:)=0; mask(end,:,:)=0; mask(:,:,end)=0; mask(:,end,:)=0;

while 1
    missclass_points_image = (mask>0) .* ( source ~= target );
    diff_before = diffCount;
    diffCount = sum(sum(sum(missclass_points_image)));
    if diff_before == diffCount
        break;
    end
    lin_ind_sp = find(missclass_points_image ~= 0);
    thisBrightness = brightness(lin_ind_sp);
    [~, idx] = sort(thisBrightness); lin_ind_sp = lin_ind_sp(idx);
    for ii = 1:length(lin_ind_sp)
        [x y z] = ind2sub(size(mask),lin_ind_sp(ii));
        patch = source(x-1:x+1, y-1:y+1, z-1:z+1);
         if simple3d(patch,6)
             source(x, y, z) = ~source(x, y, z);
         end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s=simple3d(im,n)
% Decides whether the central point of IM is a simple point
% 6-simple means that T_6(IM)=1 and T_26(~IM)=1 -- 26-simple is the opposite
% authors: Viren Jain, Mark Richardson, M.I.T.                                     

if ndims(im) ~=3
    error('image patch must be 3d')
end
if any(size(im)~=[3 3 3])
    error('must be a 3x3x3 image patch')
end

switch n
    case 6
        if topo(im,6)==1 & topo(1-im,26)==1
            s=1;
        else
            s=0;
        end
    case 26
        if topo(im,26)==1 & topo(1-im,6)==1
            s=1;
        else
            s=0;
        end
    otherwise
        error('n must be 4 or 8')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t=topo(im,n)
% Computes topological numbers for the central point of an image patch.
% These numbers can be used as the basis of a topological classification.
% T_4 and T_8 are used when IM is a 2d image patch of size 3x3
% T_6 and T_26 are used when IM is a 3d image patch of size 3x3x3,
% defined on p. 172 of Bertrand & Malandain, Patt. Recog. Lett. 15, 169-75 (1994).
% authors: Viren Jain, Mark Richardson, M.I.T.

switch n
    case 4
        % number of 4-connected components in the 8-neighborhood of the
        % center that are 4-adjacent to the center
        if ndims(im) ~= 2
            error('n=4 is valid for a 2d image')
        end
        if any(size(im)~=[3 3])
            error('must be 3x3 image patch')
        end
        neighbor4=[0 1 0; 1 0 1; 0 1 0];
        im(2,2)=0;    % ignore the central point
        components=bwlabel(im,4).*neighbor4;  % zero out locations that are not in the four-neighborhood
    case 8
        % number of 8-connected components in the 8-neighborhood of the
        % center (adjacency is automatic)
        if ndims(im) ~= 2
            error('n=8 is valid for a 2d image')
        end
        if any(size(im)~=[3 3])
            error('must be 3x3 image patch')
        end
        im(2,2)=0;  % ignore the central point
        components=bwlabel(im,8);
    case 6
        % number of 6-connected components in the 18-neighborhood of the center
        % that are 6-adjacent to the center
        if ndims(im) ~= 3
            error('n=6 is valid for a 3d image')
        end
        if any(size(im)~=[3 3 3])
            error('must be 3x3x3 image patch')
        end
        neighbor6=conndef(3,'min');  % the nearest six neighbors
        neighbor18=ones(3,3,3); neighbor18(1:2:3,1:2:3,1:2:3)=0; neighbor18(2,2,2)=0;   % the nearest 18 neighbors
        components=bwlabeln(neighbor18.*im,6);  % 6-connected components in the 18 neighborhood of the center
        components=components.*neighbor6;  % keep components that are 6-adjacent to the center
    case 26
        % number of 26-components in the 26-neighborhood of the center
        % (adjacency is automatic)
        if ndims(im) ~= 3
            error('n=26 is valid for a 3d image')
        end
        if any(size(im)~=[3 3 3])
            error('must be 3x3x3 image patch')
        end
        im(2,2,2)=0;
        components=bwlabeln(im,26);
    otherwise
        error('n must be 4, 8, 6, or 26')
end
t=length(unique(nonzeros(components)));
