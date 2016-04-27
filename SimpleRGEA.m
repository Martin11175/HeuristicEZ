function [ G ] = SimpleRGEA( J, K, ...
    min_AP_strength, num_strongest_APs, min_strong_overlap )
%SimpleRGEA Simple Relative Gain Estimation Algorithm using strongest RSSIs
%   J [in] - Matrix of RSSI measurements (APs as columns)
%   K [in] - Vector of device IDs relating rows to device
%   G [out] - Pairs of device IDs (1) to estimated gain (2)

% Author: Martin Higgs (martinhiggs85@googlemail.com)

if exist('num_strongest_APs', 'var') == 0
    num_strongest_APs = 4;
end
if exist('min_strong_overlap', 'var') == 0
    min_strong_overlap = 3;
end
if exist('min_AP_strength', 'var') == 0
    min_AP_strength = -90;
end

D = sort(unique(K)); % List of device IDs
J(J == 100) = -100; % Replace positive invisibility markers
J(J < min_AP_strength) = -100;

% Isolate strongest RSSI measurements per location
max_J = zeros(size(J,1), num_strongest_APs);
for i = 1:size(J,1)
    [~, order] = sort(J(i,:), 'descend');
    max_J(i,:) = order(1:num_strongest_APs);
    %J(i, order(num_strongest_APs:size(J,2))) = 0;
end
%J = sparse(J);

% Calculate relative gain between pairs of devices
deltaG = zeros(size(D,1));
sigma_deltaG = zeros(size(D,1)); % Uncertainty (estimated standard deviation)

for i = nchoosek(1:size(D,1), 2)'
    % Test pairs of measurements for proximity by similarity
    max_1 = max_J(K == D(i(1)), :);
    k_1 = J(K == D(i(1)), :);
    max_2 = max_J(K == D(i(2)), :);
    k_2 = J(K == D(i(2)), :);
    
    avg_diff = zeros(size(max_1, 1), size(max_2, 1));
    prox = false(size(max_1, 1), size(max_2, 1));

    for m = 1:size(max_1,1)
        if sum(k_1(m, max_1(m,:))) > (num_strongest_APs * min_AP_strength)
        for n = 1:size(max_2,1)
            if sum(ismember(max_1(m,:), max_2(n,:))) > min_strong_overlap
                prox(m, n) = true;
            end
        end
        end
    end
    [m, n] = find(prox);
    
    if size(m, 1) > 0
        for p = [m n]'
            % Only compare visible APs
            vis = (k_1(p(1),:) > -100) & (k_2(p(2),:) > -100);
            if sum(vis) > 0
                avg_diff(p(1), p(2)) = mean(k_1(p(1), vis) - k_2(p(2), vis));
            end
        end
        
        deltaG(i(1), i(2)) = mean(avg_diff(prox));
        sigma_deltaG(i(1), i(2)) = (1 / size(m, 1)) * sqrt(sum((avg_diff(prox) - deltaG(i(1), i(2))).^2));
    end
    fprintf('%f | %d:%d / %d\n', deltaG(i(1), i(2)), i(1), i(2), size(D,1))
end

% Solve least mean squares set of simultaneous equations
[i, j] = find(deltaG);
C = zeros(size(i,1), size(D,1));
d = zeros(size(i,1),1);

% Prepare system of simultaneous equations weighted by estimated standard deviation
for k = 1:size(i,1)
    % TODO: Scale weight, sigma^-1 increases exponentially as reaching peak
    % of bell curve (relatively little confidence improvement)!
    weight = 1 / sigma_deltaG(i(k), j(k));
    C(k, [i(k) j(k)]) = 1 * weight;
    d(k) = deltaG(i(k),j(k)) * weight;
end

warning('off','MATLAB:rankDeficientMatrix');
G = [D (C \ d)];
warning('on','MATLAB:rankDeficientMatrix');

end
