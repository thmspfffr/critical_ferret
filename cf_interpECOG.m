function [xi yi zi] = is_interpECOG(imat,interpMethod)
% This function interpolates ECOG data.
% Input:    - imat(64,n-Observations)
% Output:   - xi = x coordinates
%           - yi = y coordinates
%           - zi = interpolated ECOG values
% I.S. 2013

xpos = [4 6 8 10 12 ...          % first row
        3 5 7 9 11 13 15 ...     % second row
        4 6 8 10 12 14 16 ...    % third row
        3 5 7 9 11 13 15 17 ...  % fourth row
        2 4 6 8 10 12 14 16 ...  % fifth row
        3 5 7 9 11 13 15 17 ...  % sixth row
        2 4 6 8 10 12 14 16 ...  % seventh row
        1 3 5 7 9 11 13 15 ...   % eigth row
        4 6 8 10 12];            % ninth row

ypos = [1 1 1 1 1 ...          % first row
        2 2 2 2 2 2 2 ...      % second row
        3 3 3 3 3 3 3 ...      % third row
        4 4 4 4 4 4 4 4 ...    % fourth row
        5 5 5 5 5 5 5 5 ...    % fifth row
        6 6 6 6 6 6 6 6 ...    % sixth row
        7 7 7 7 7 7 7 7 ...    % seventh row
        8 8 8 8 8 8 8 8 ...    % eigth row
        9 9 9 9 9];            % ninth row    
    
    % Change coordinates to mm's
    xpos = xpos*0.75;
    ypos = ypos*1.4142;
    
    % get meshgrid
    [xi yi] = meshgrid(0:0.05:(max(xpos)+0.2),0:0.05:(max(ypos)+1.4142));
    
    % interpolate data
    zi = griddata(xpos,ypos,imat,xi,yi,interpMethod);
end