function [xF,yF] = helper_plotbestfitline(x,y)
ind = ~isnan(x) & ~isnan(y);
% Find fit
coefficients = polyfit(x(ind), y(ind), 1);
% Create a new x axis with exactly 1000 points 
xF = linspace(min(x(ind)), max(x(ind)), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yF = polyval(coefficients , xF);