function [label, centroid, dis] = colorInitializedWeightedKmeans(X, k, initialLabels, options)

% ADAPTED BY UYGAR SUMBUL FROM ...

% FKMEANS Fast K-means with optional weighting and careful initialization.
% [L, C, D] = FKMEANS(X, k) partitions the vectors in the n-by-p matrix X
% into k (or, rarely, fewer) clusters by applying the well known batch
% K-means algorithm. Rows of X correspond to points, columns correspond to
% variables. The output k-by-p matrix C contains the cluster centroids. The
% n-element output column vector L contains the cluster label of each
% point. The k-element output column vector D contains the residual cluster
% distortions as measured by total squared distance of cluster members from
% the centroid.
%
% FKMEANS(X, C0) where C0 is a k-by-p matrix uses the rows of C0 as the
% initial centroids instead of choosing them randomly from X.
%
% FKMEANS(X, k, options) allows optional parameter name/value pairs to 
% be specified. Parameters are:
%
%   'weight' - n-by-1 weight vector used to adjust centroid and distortion
%              calculations. Weights should be positive.
%   'careful' - binary option that determines whether "careful seeding"
%               as recommended by Arthur and Vassilvitskii is used when
%               choosing initial centroids. This option should be used
%               with care because numerical experiments suggest it may
%               be counter-productive when the data is noisy.
%
% Notes
% (1) The careful seeding procedure chooses the first centroid at random
% from X, and each successive centroid from the remaining points according
% to the categorical distribution with selection probabilities proportional
% to the point's minimum squared Euclidean distance from the already chosen
% centroids. This tends to spread the points out more evenly, and, if the
% data is made of k well separated clusters, is likely to choose an initial
% centroid from each cluster. This can speed convergence and reduce the
% likelihood of getting a bad solution [1]. However, in experiments where
% 5% uniformly distributed noise data was added to such naturally clustered
% data the results were frequently worse then when centroids were chosen at
% random.
% (2) If, as is possible, a cluster is empty at the end of an iteration,
% then there may be fewer than k clusters returned. In practice this seems
% to happen very rarely.
% (3) Unlike the Mathworks KMEANS this implementation does not perform a
% final, slow, phase of incremental K-means ('onlinephase') that guarantees
% convergence to a local minimum. 
%
% References
% [1] "k-means++: The Advantages of Careful Seeding", by David Arthur and
% Sergei Vassilvitskii, SODA 2007.

n = size(X,1);

% option defaults
weight = 0; % uniform unit weighting
maxIter = 200;

if nargin == 4
    if isfield(options, 'weight')
        weight = options.weight;
    end
    if isfield(options,'maxIter')
        maxIter = options.maxIter;
    end
end

% generate initial labeling of points
label = initialLabels;

last = zeros(size(label));
iter = 0;
if ~weight
    % code defactoring for speed
    while any(label ~= last) & iter<maxIter
        % remove empty clusters
        [~,~,label] = unique(label);
        % transform label into indicator matrix
        ind = sparse(label,1:n,1,k,n,n);
        % compute centroid of each cluster
        centroid = (spdiags(1./sum(ind,2),0,k,k)*ind)*X;
        % compute distance of every point to each centroid
        distances = bsxfun(@minus,centroid*X',0.5*sum(centroid.^2,2));
        % assign points to their nearest centroid
        last = label;
        [~,label] = max(distances);
        label = label';
        iter = iter + 1;
    end
    dis = ind*(sum(X.^2,2) - 2*max(distances)');
else


    while any(label ~= last) & iter<maxIter
        % remove empty clusters
        [~,~,label] = unique(label);
        % transform label into indicator matrix
        ind = sparse(label,1:n,weight,k,n,n);
        % compute centroid of each cluster
        centroid = (spdiags(1./sum(ind,2),0,k,k)*ind)*X;
        % compute distance of every point to each centroid
        distances = bsxfun(@minus,centroid*X',0.5*sum(centroid.^2,2));
        % assign points to their nearest centroid
        last = label;
        [~,label] = max(distances);
        label = label';
        iter = iter + 1;
end
dis = ind*(sum(X.^2,2) - 2*max(distances)');
end
label = label';

