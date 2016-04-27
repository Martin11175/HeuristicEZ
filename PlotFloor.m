function PlotFloor( I, J )
%PLOTFLOOR Plot a floor of AP and measurement locations
%   I [in] - Vector of AP parameters:
%       c longtitude, c latitude, transmit power, path loss rate
%   J [in] - Matrix of measurement locations (X_long, X_lat)

% Author: Martin Higgs (martinhiggs85@googlemail.com)

figure
hold on

% Plot measured locations as points
for j = J'
    plot(j(1), j(2), '.k')
end

% Colour and topography ring definition
colours = ['r' 'y' 'g' 'c' 'b'];
angle=0:0.01:2*pi;

% Plot Access Points as topographical rings
I = I((I(:,4) ~= 0),:); % Remove undefined APs
for i = I'
    plot(i(1), i(2), '*r')
    
    n = 1;
    for rss = [-40 -55 -70 -85 -100]
        % Distance = strength difference divided by path loss rate
        radius = (rss - i(3)) / i(4);
        xrotation = radius * cos(angle);
        yrotation = radius * sin(angle);
        plot(i(1) + xrotation, i(2) + yrotation, colours(n));
        
        n = n + 1;
    end
end

hold off

end

