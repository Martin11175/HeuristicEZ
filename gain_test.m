% Relative Gain Estimation Test
% Author: Martin Higgs (martinhiggs85@googlemail.com)

% Load raw data from CSV
srcData = csvread('UJIndoorLoc/trainingData.csv', 1);

% Scale known location data to a local space
min_lat = min(srcData(:,521));
min_long = min(srcData(:,522));
srcData(:,521) = srcData(:,521) - min_lat;
srcData(:,522) = srcData(:,522) - min_long;

% Prepare summary files
if exist('ez_relgain_summary.csv', 'file') == 0
    fprintf(fopen('ez_relgain_summary.csv', 'w'), ...
        'floor,building,prox_threshold,min_AP_overlap,gain_error\n');
end
if exist('h_relgain_summary.csv', 'file') == 0
    fprintf(fopen('h_relgain_summary.csv', 'w'), ...
        'floor,building,num_strongest_APs,min_AP_strength,min_strong_overlap,gain_error\n');
end

% Parameters to test
% EZ
prox_threshold = [ 2 3 4 5 7 10 ];
min_AP_overlap = [ 1 2 3 4 5 ];
% Heuristic
min_AP_strength = [ -100 -90 -80 -70 ];
num_strongest = [ 2 3 4 5 7 10 ];
min_strong_overlap = [ 2 3 4 5 6 ];

% Test accuracy -----------------------------------------------------------
for building = [ 0 1 2 ]
for floor = [ 0 1 2 3 ]
    
    dataSet = ((srcData(:,523) == floor) & (srcData(:,524) == building));
    %dataSet = (srcData(:,524) == building);
    %dataSet = (srcData(:,523) == floor);
    
    % Set data up for test
    J = srcData(dataSet, 1:520);
    K = srcData(dataSet, 528);
    D = sort(unique(K));
    X = srcData(dataSet, 521:523);
    
    % Calculate base line estimates
    g = GroundRGEA(J,K,X);
    
    % Run algortihms and write summaries
    for prox = prox_threshold
        for overlap = min_AP_overlap
            fprintf('F(%d) B(%d) | EZ: %d %d\n', floor, building, prox, overlap);
            ez = RGEA(J,K,prox,overlap);
            
            dlmwrite('ez_relgain_summary.csv', ...
                [floor, building, prox, overlap, sum(sum(abs(g - ez)))], ...
                'delimiter', ',', '-append');
        end
    end

    for strength = min_AP_strength
        for num_str = num_strongest
            for str_overlap = min_strong_overlap
                fprintf('F(%d) B(%d) | H: %d %d %d\n', floor, building, num_str, strength, str_overlap);
                h = SimpleRGEA(J,K, strength, num_str, str_overlap);
                
                dlmwrite('h_relgain_summary.csv', ...
                    [floor, building, num_str, strength, str_overlap, sum(sum(abs(g - h)))], ...
                    'delimiter', ',', '-append');
            end
        end
    end

end
end
