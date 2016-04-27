% Main Matlab script file
% Author: Martin Higgs (martinhiggs85@googlemail.com)

%----------------------------Parameter Setup------------------------------%

% Load raw data from CSV
% TODO: Enable dynamic size scaling
% Algorithms generally set for 520 APs, should accept 100 (positive
% invisibility marker) for unused access points as a workaround for smaller
% data sets.
if exist('srcData', 'var') == 0
    srcData = csvread('UJIndoorLoc/trainingData.csv', 1);
end
if exist('testData', 'var') == 0
    testData = csvread('UJIndoorLoc/validationData.csv', 1);
end
if exist('chance', 'var') == 0
    chance = 100;
end

% Algorithm types
if exist('RGEA_type', 'var') == 0
    RGEA_type = 'ground';
end
if exist('apselect_type', 'var') == 0
    apselect_type = 'overlap';
end
if exist('locselect_type', 'var') == 0
    locselect_type = 'none';
end
if exist('ips_type', 'var') == 0
    ips_type = 'all';
end

% Scope definition
if exist('bounds', 'var') == 0
    bounds = 50;
end
if exist('floors', 'var') == 0
    floors = 0;
end
if exist('buildings', 'var') == 0
    buildings = 0;
end
if exist('thresholds', 'var') == 0
    thresholds = -100;
end

% Scale known location data to a local space
min_lat = min(srcData(:,521));
min_long = min(srcData(:,522));
srcData(:,521) = srcData(:,521) - (min_lat - bounds);
srcData(:,522) = srcData(:,522) - (min_long - bounds);
testData(:,521) = testData(:,521) - (min_lat - bounds);
testData(:,522) = testData(:,522) - (min_long - bounds);

%-----------------------------EZ-ALGORITHM--------------------------------%

% Iterate over 2D subspaces
for floor = floors
for building = buildings
for threshold = thresholds

% Isolate a 2D sub-space (single floor of single building)
rows = ((srcData(:,523) == floor) & (srcData(:,524) == building));
dataSet = srcData(rows, :);
normalisedData = 1 - (abs(dataSet(:, 1:520)) / 100);
%normalisedData = normalisedData - mean(normalisedData);
original_locations = dataSet(:, 521:522);

% Relative Gain Estimation Algorithm
if strcmp(RGEA_type,'none') == 0
    tic
    if strcmp(RGEA_type,'rgea') == 1
        rgea_num = 1;
        G = RGEA(dataSet(:, 1:520), dataSet(:, 528));
    elseif strcmp(RGEA_type,'simple') == 1
        rgea_num = 2;
        G = SimpleRGEA(dataSet(:, 1:520), dataSet(:, 528));
    else
        rgea_num = 3;
        G = GroundRGEA(dataSet(:, 1:520), dataSet(:, 528), dataSet(:,521:523));
    end
    for d = G'
        visMeas = (dataSet(:,1:520) ~= 100) & (dataSet(:,1:520) ~= -100);
        deviceMeas = dataSet(:,528) == d(1);
        visMeas(~deviceMeas, :) = false;
        totMeas = [visMeas, false(size(visMeas,1), 9)];
        dataSet(totMeas) = dataSet(totMeas) - d(2);
    end
    toc
else
    rgea_num = 0;
end

% Artificial GPS restriction
for j = 1:size(dataSet,1)
    if randi(100) > chance
        dataSet(j, 521:522) = [0 0];
    end
end

% APSelect Algorithm
if strcmp(apselect_type,'none') == 0
    tic
    ap_rank = zeros(520, 1);
    for i = 1:520 % Rank by number of visible known locations
        ap_rank(i) = sum(dataSet(normalisedData(:,i) > 0, 521) > 0);
    end
    if strcmp(apselect_type,'all') == 1
        ap_num = 1;
        APSelect = HeirarchicalCluster(normalisedData, ap_rank, 'all');
    elseif strcmp(apselect_type,'max') == 1
        ap_num = 2;
        APSelect = HeirarchicalCluster(normalisedData, ap_rank, 'max');
    else
        ap_num = 3;
        APSelect = HeirarchicalCluster(normalisedData, ap_rank, 'overlap');
    end
    toc
else
    ap_num = 0;
    APSelect = 1:520;
end

% LocSelect Algorithm
if strcmp(locselect_type,'none') == 0
    tic
    % Prefer GPS localised measurements (binary rank)
    loc_rank = (dataSet(:,521) ~= 0 | dataSet(:,522) ~= 0);
    if strcmp(locselect_type,'all') == 1
        loc_num = 1;
        LocSelect = HeirarchicalCluster(normalisedData', loc_rank, 'all');
    elseif strcmp(locselect_type,'max') == 1
        loc_num = 2;
        LocSelect = HeirarchicalCluster(normalisedData', loc_rank, 'max');
    else
        loc_num = 3;
        LocSelect = HeirarchicalCluster(normalisedData', loc_rank, 'overlap');
    end
    toc
else
    loc_num = 0;
    LocSelect = 1:size(normalisedData,1);
end
LocSelect = ismember(1:size(normalisedData,1), LocSelect)';

% Iterate until no new parameters can be determined
APparams = zeros(520, 5);
new_value_flag = true;
while new_value_flag
    new_value_flag = false;

    % Determine new AP parameters from localised measurements
    for i = APSelect(APparams(APSelect, 4) == 0)
        % Find which observations include our desired AP
        J = ((threshold < dataSet(:,i)) & (dataSet(:,i) < 0) ...
            & ((dataSet(:,521) ~= 0) | (dataSet(:,522) ~= 0)) ...
            & LocSelect);
        O = dataSet(J,[i 521 522]);
        
        if size(O,1) > 4
            new_value_flag = true;
            % Objective function looking to minimise RSSI error
            objective = @(C) mean((O(:,1) - C(3) + ((10 * C(4)) * ...
                log10(sqrt((O(:,2) - C(1)).^2 + (O(:,3) - C(2)).^2)))).^2);
            
            % Lower and upper bounds for parameter search:
            %   Lat and long (in metres): bounds outside of observation locations
            %   Transmit power: 0 to -50dBm
            %   Path loss: 1.5 to 6.0 from EZ's trials
            lower = [min(O(:,2)) - bounds, min(O(:,3)) - bounds, -50, 1.5];
            upper = [max(O(:,2)) + bounds, max(O(:,3)) + bounds, 0, 6.0];
            
            % Perform simulated annealing, starting from average values
            C0 = [mean(O(:,2)), mean(O(:,3)), -25, 3.0];
            options = optimset('display','off');
            fprintf('AP: %d, %d observations...', i, size(O,1));
            med_ap_err = Inf;
            try_count = 0;
            tic
            % Avoid obviously terrible results
            while med_ap_err > 10 && try_count < 10
                APparams(i,1:4) = simulannealbnd(objective, C0, lower, upper, options);
                % Trust metric
                med_ap_err = median(abs((sqrt((O(:,2) - APparams(i,1)).^2 + (O(:,3) - APparams(i,2)).^2) ...
                    - 10.^((APparams(i,3) - O(:,1))./(10*APparams(i,4))))));
                fprintf('Median localisation error: %fm\n', med_ap_err);
                try_count = try_count + 1;
            end
            toc
            if med_ap_err > 10
                APparams(i,1:5) = [0 0 0 0 0];
            else
                APparams(i, 5) = med_ap_err;
            end
        end
    end
    
    % If new APs are parameterised, check to localise new measurements
    if new_value_flag
        new_value_flag = false;
        
        for j = dataSet((dataSet(:,521) == 0) & (dataSet(:,522) == 0) & LocSelect)'
            % Find which APs we wish to evaluate
            I = ((threshold < j(1:520)) & (j(1:520) < 0) & (APparams(:,4) ~= 0));
    	    if (strcmp(ips_type,'max') == 1) && (sum(I) > 2)
                tmp = find(I);
                [~, order] = sort(j(I), 'descend');
                I = false(520,1);
                I(tmp(order(1:3))) = true;
            end
            O = j(I);
            if strcmp(RGEA_type,'none') == 0
                if sum(G(:,1) == j(528)) > 0
                    O = O - G(G(:,1) == j(528),2);
                end
            end
            C = APparams(I, :);
            
            if size(O,1) > 2
                new_value_flag = true;
                % Calculate distance and set equation to solve
                D = 10.^((C(:,3) - O)./(10*C(:,4)));
                objective = @(J) D - sqrt((J(1) - C(:,1)).^2 + (J(2) - C(:,2)).^2);
                
                % Bounds on area to search (outside AP locations)
                lower = [min(C(:,1)) - bounds, min(C(:,1)) - bounds];
                upper = [max(C(:,2)) + bounds, max(C(:,2)) + bounds];
                
                % Perform simulated annealing, starting from average values
                J0 = [mean(C(:,1)), mean(C(:,2))];
                options = optimset('display','off');
                dataSet(j,521:522) = lsqnonlin(objective, J0, lower, upper, options);
            end
        end
    end
end

%------------------------------EVALUATION---------------------------------%

% Isolate a 2D sub-space (single floor of single building)
rows = ((testData(:,523) == floor) & (testData(:,524) == building));
test_dataSet = testData(rows, :);

% Localise validation data
IPSresults = zeros(size(test_dataSet,1),2);
IPSerror = zeros(size(test_dataSet,1),1);
n = 1;
for j = test_dataSet'
    % Find which APs we wish to evaluate
    I = ((threshold < j(1:520)) & (j(1:520) < 0) & (APparams(:,4) ~= 0));
    if (strcmp(ips_type,'max') == 1) && (sum(I) > 2)
        tmp = find(I);
        [~, order] = sort(j(I), 'descend');
        I = false(520,1);
        I(tmp(order(1:3))) = true;
    end
    O = j(I);
    if strcmp(RGEA_type,'none') == 0
        if sum(G(:,1) == j(528)) > 0
            O = O - G(G(:,1) == j(528),2);
        end
    end
    C = APparams(I, :);
        
    if size(O,1) > 2
        % Calculate distance and set equation to solve
        D = 10.^((C(:,3) - O)./(10*C(:,4)));
        objective = @(J) D - sqrt((J(1) - C(:,1)).^2 + (J(2) - C(:,2)).^2);
        
        % Bounds on area to search (outside AP locations)
        lower = [min(C(:,1)) - bounds, min(C(:,1)) - bounds];
        upper = [max(C(:,2)) + bounds, max(C(:,2)) + bounds];
    
        % Perform simulated annealing, starting from average values
        J0 = [mean(C(:,1)), mean(C(:,2))];
        options = optimset('display','off');
        IPSresults(n,:) = lsqnonlin(objective, J0, lower, upper, options);
        
        % Calculate localisation error
        IPSerror(n) = sqrt((IPSresults(n,1) - test_dataSet(n,521))^2 + (IPSresults(n,2) - test_dataSet(n,522))^2);
    end
    n = n + 1;
end

% Summary results output
if strcmp(ips_type, 'all') == 1
    ips_num = 0;
else
    ips_num = 1;
end

fprintf('%d: %d/%d localised. %f median IPS error\n', ...
    threshold, sum(IPSerror > 0), size(IPSerror,1), median(IPSerror(IPSerror > 0)));
if exist('summary.csv', 'file') == 0
    fprintf(fopen('summary.csv', 'w'), ...
        'building, floor, rgea, ap, loc, ips, threshold, bounds, chance, num_test, num_loc, median_error\n');
end
dlmwrite('summary.csv', ...
    [building, floor, rgea_num, ap_num, loc_num, ips_num, ...
    threshold, bounds, chance, ...
    size(IPSerror,1), sum(IPSerror > 0), median(IPSerror(IPSerror > 0))], ...
    'delimiter', ',', '-append');

% Create results directory structure
if exist('results', 'dir') == 0
    mkdir('results');
end
building_dir = sprintf('results/building_%d', building);
if exist(building_dir, 'dir') == 0
    mkdir(building_dir);
end
floor_dir = sprintf('results/building_%d/floor_%d', building, floor);
if exist(floor_dir, 'dir') == 0
    mkdir(floor_dir);
end

% Output csv files of results
src_filename = sprintf('results/building_%d/floor_%d/src_rgea(%s)_ap(%s)_loc(%s)_ips(%s)_threshold(%d)_bounds(%d)_chance(%d).csv', ...
    building, floor, RGEA_type, apselect_type, locselect_type, ips_type, threshold, bounds, chance);
csvwrite(src_filename, [original_locations dataSet(:, 521:522) LocSelect]);

tst_filename = sprintf('results/building_%d/floor_%d/tst_rgea(%s)_ap(%s)_loc(%s)_ips(%s)_threshold(%d)_bounds(%d)_chance(%d).csv', ...
    building, floor, RGEA_type, apselect_type, locselect_type, ips_type, threshold, bounds, chance);
csvwrite(tst_filename, [test_dataSet(:,521:522) IPSresults]);

ap_filename = sprintf('results/building_%d/floor_%d/ap_rgea(%s)_ap(%s)_loc(%s)_ips(%s)_threshold(%d)_bounds(%d)_chance(%d).csv', ...
    building, floor, RGEA_type, apselect_type, locselect_type, ips_type, threshold, bounds, chance);
csvwrite(ap_filename, APparams);
end
end
end
