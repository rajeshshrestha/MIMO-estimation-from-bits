%Define parameters
function xy =angle_generation(nPoints, lim, minDist)
% nPoints = 20;       %number of coordinates
% lim = [100,1900];   %bounds of random numbers
% minDist = 20;      %minimum distance between turbines

% Create random coordinates and continually replace coordinates
% that are too close to another point.  Stop when minimum distance
% is satisfied or after making nPoints*100 attempts.
xy = nan(nPoints,1);
c = 0; %Counter
while any(isnan(xy(:))) && c<(nPoints*10000000)
    % Fill NaN values with new random coordinates
    xy(isnan(xy)) = rand(1,sum(isnan(xy(:)))) * (lim(2)-lim(1)) + lim(1);
    % Identify rows that are too close to another point
    [~,isTooClose] = find(triu(squareform(pdist(xy)) < minDist,1));
    % Replace too-close coordinates with NaN
    xy(isTooClose,:) = NaN; 
    c = c+1;
end

% Throw error if the loop had to quit and missing values remain
if any(isnan(xy(:)))
    error('The while-loop gave up. There are %d coordinates with missing values.',sum(isnan(xy(:,1))))
end
% 
% % Display number of attempts 
% fprintf('%d number of attempts.\n', c)
% 
% Show the minimum distance 
distances = squareform(pdist(xy)); 
% fprintf('Min distance = %.2f\n', min(distances(distances~=0))/pi*180)

% % Plot results
%  figure() 
%  plot(xy(:,1),xy(:,2), 'ks', 'MarkerSize', 10)
%  grid on