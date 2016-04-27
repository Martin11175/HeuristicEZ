function [ G ] = GroundRGEA( J, K, X, prox_dst )
%GroundRGEA Ground Truth Relative Gain Estimation Algorithm using location data
%   J [in] - Matrix of RSSI measurements (APs as columns)
%   K [in] - Vector of device IDs relating rows to device
%   X [in] - Matrix of locations for RSSI measurements
%       [latitude, longitude, floor]
%   G [out] - Pairs of device IDs (1) to estimated gain (2)

% Author: Martin Higgs (martinhiggs85@googlemail.com)

D = sort(unique(K)); % List of device IDs
J(J == 100) = -100; % Replace positive invisibility markers to prevent skew

% Calculate relative gain between pairs of devices
deltaG = zeros(size(D,1));
sigma_deltaG = zeros(size(D,1)); % Uncertainty (estimated standard deviation)
tot_prox = 0;
if exist('prox_dst', 'var') == 0
    prox_dst = 1;
end

for i = nchoosek(1:size(D,1), 2)'
    k_1 = J(K == D(i(1)), :);
    x_1 = X(K == D(i(1)), :);
    k_2 = J(K == D(i(2)), :);
    x_2 = X(K == D(i(2)), :);
    
    % Test pairs of measurements for proximity by similarity
    avg_diff = zeros(size(k_1, 1), size(k_2, 1));
    
    % Proximate on the same floor within 1m
    prox = (abs(bsxfun(@minus, x_1(:,1), x_2(:,1)')) < prox_dst ...
        & abs(bsxfun(@minus, x_1(:,2), x_2(:,2)')) < prox_dst ...
        & bsxfun(@eq, x_1(:,3), x_2(:,3)'));
    [m, n] = find(prox);
    
    if size(m, 1) > 0
        tot_prox = tot_prox + size(m,1);
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
    %fprintf('%f | %d:%d / %d\n', deltaG(i(1), i(2)), i(1), i(2), size(D,1))
end

% Print out device observation connectedness
%{
fprintf('\n');
deltaG
fprintf('\n');
for i = 1:size(deltaG,1);
    fprintf('%d: ', i);
    for m = sort(unique(([find(deltaG(i,:)) find(deltaG(:,i))'])))
        fprintf('%d ',m);
    end
    fprintf('\n');
end
%}

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
