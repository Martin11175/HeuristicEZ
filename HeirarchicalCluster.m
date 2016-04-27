function [ Z ] = HeirarchicalCluster( X, rank, compare_type, ...
    threshold_similarity, min_overlap, num_max)
%HEIRARCHICALCLUSTER Clustering as specified by EZ for APs and locations
%   X [in] - Normalised data matrix, clusters by columns
%   rank [in] - Ranking of all columns indicating preference when selecting 
%       for cluster representative (higher rank = higher preference)
%   compare_type [optional] - Data sets to compare for clustering:
%       'all' (default) - Compare all values per column, including 0
%       'overlap' - Compare only values > 0 in either column
%       'max' - Compare only similar maximum strength values
%   threshold_similarity [optional] - Threshold at which to cluster (0 - 1)
%   min_overlap [optional] - Number of measurements necessary to compare columns
%   num_max [optional] - Number of maximum measurements to compare
%   Z [out] - Array of columns elected to represent clusters

% Author: Martin Higgs (martinhiggs85@googlemail.com)

Z = 1:size(X,2); % Mapping of all columns to clusters
max_similarity = 1;
if exist('theshold_similarity', 'var') == 0
    threshold_similarity = 0.9;
end
if exist('min_overlap', 'var') == 0
    min_overlap = 3;
end

if exist('compare_type', 'var') == 0 || strcmp(compare_type, 'all')
    compare_type = 1;
elseif strcmp(compare_type, 'overlap')
    compare_type = 2;
elseif strcmp(compare_type, 'max')
    compare_type = 3;
    if exist('num_max', 'var') == 0
        num_max = min_overlap + 1;
    elseif num_max < min_overlap
        error('Minimum overlap requested is less than the number of strongest' ...
            + ' measurements specified');
    end
else
    error('Unrecognised comparison type for clustering');
end

% Prepare sets of maximum strength values for comparison
if compare_type == 3
    max_X = zeros(num_max, size(X,2));
    for i = 1:size(X,2)
        [~, order] = sort(X(:,i), 'descend');
        max_X(:,i) = order(1:num_max);
    end
end

% Immediately cluster all empty columns
if compare_type == 1
    Z(mean(X,1) < (max_similarity - threshold_similarity)) = 0;
else
    Z(sum(X ~= 0,1) < min_overlap) = 0;
end

% Iterate until no more similar clusters
while max_similarity > threshold_similarity
    max_similarity = threshold_similarity; % Reset maximum observed similarity between clusters
    C = sort(unique(Z(Z ~= 0))); % Clusters to be evaluated
    
    if size(C) == 1 % Clustered to a single point
        break;
    end
    
    % Calculate similarity and mark maximum on each pass
    for i = nchoosek(C, 2)'
        % Similarity between clusters is average of all inter-cluster pairs
        pair_sim = zeros(sum(Z == i(1)), sum(Z == i(2)));
        count_x = 1;
        for x = find(Z == i(1))
            count_y = 1;
            for y = find(Z == i(2))
                if compare_type == 1
                    overlap = true(size(X,1),1);
                elseif compare_type == 2
                    overlap = (X(:,x) > 0) & (X(:,y) > 0);
                    if sum(overlap) < min_overlap
                        continue; % If not enough overlap, pair_sim = 0
                    end
                elseif compare_type == 3
                    overlap = ismember(max_X(:,x), max_X(:,y));
                    if sum(overlap) < min_overlap
                        continue; % If not enough overlap, pair_sim = 0
                    end
                    overlap = max_X(overlap,x);
                end
                pair_sim(count_x, count_y) = 1 - mean(abs(X(overlap,x) - X(overlap,y)));
                count_y = count_y + 1;
            end
            count_x = count_x + 1;
        end
        similarity = mean2(pair_sim);
        
        if similarity > max_similarity
            max_similarity = similarity;
            cluster = i; % Pair of clusters to merge
        end
            
        % Nothing can be more similar than exact (speed boost)
        if max_similarity == 1
            break;
        end
    end
    
    % Merge most similar cluster
    if max_similarity > threshold_similarity
        Z(Z == cluster(2)) = cluster(1);
        %fprintf('Merging %d / %d. %d clusters remaining. Similarity: %f\n', cluster(1), cluster(2), size(unique(Z), 2), max_similarity);
    end
end

% Select representative column for each cluster
for c = sort(unique(Z(Z ~= 0)))
    % First select based on the max rank for the cluster
    rep = find(rank == max(rank(Z == c)) & (Z == c)');
    
    % Check for similarity to cluster if multiple with equivalent rank
    if size(rep,1) > 1
        best_avg_similarity = 0;
        
        % Calculate average similarity to the rest of the cluster
        for i = rep'
            avg_similarity = 0;
            for j = find(Z == c)
                if j ~= i
                    if compare_type == 1
                        overlap = true(size(X,1),1);
                    elseif compare_type == 2
                        overlap = (X(:,i) > 0) & (X(:,j) > 0);
                        if sum(overlap) < min_overlap
                            continue; % If not enough overlap, pair_sim = 0
                        end
                    elseif compare_type == 3
                        overlap = ismember(max_X(:,i), max_X(:,j));
                        if sum(overlap) < min_overlap
                            continue; % If not enough overlap, pair_sim = 0
                        end
                        overlap = max_X(overlap,i);
                    end
                    avg_similarity = avg_similarity + (1 - mean(abs(X(overlap,i) - X(overlap,j))));
                end
            end
            avg_similarity = avg_similarity / (sum(Z == c) - 1);
            
            if avg_similarity > best_avg_similarity
                best_avg_similarity = avg_similarity;
                rep = i;
            end
        end
    end
    
    % Replace cluster representative
    Z(Z == c) = rep;
end

Z = sort(unique(Z(Z ~= 0)));

end

