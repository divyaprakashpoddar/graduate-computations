function [isn, icn, xb, yb, nb, X] = defineEllipse(a,b,nx,ny,ci,cj,theta)
% [isn, icn, xb, yb, nb, X] = defineEllipse(a,b,nx,ny,ci,cj,theta)
%   Identifies the boundary nodes of an ellipse and stores their locations
%
% input:
%   a           = Semi-major axis
%   b           = Semi-minor axis
%   nx          = Grid size in x-direction
%   ny          = Grid size in y-direction
%   ci          = Ellipse center (x-coordinate)
%   cj          = Ellipse center (y-coordinate)
%   theta       = Angle by which the ellipse has been rotated
% output:
%   isn         = A matrix containing flags (1 for wall and cylinder and 0 for fluid)
%   icn         = A matrix containing flags (1 for cylinder and 0 for fluid)
%   xb          = A vector containing the x-coordinate of the boundary nodes
%   yb          = A vector containing the y-coordinate of the boundary nodes
%   nb          = Total number of boundary nodes
%   X           = A vector containing the x and y coordinates of the cylinder center
%
% Author: Divyaprakash (2020AMZ8461), PhD Candidate
%         Applied Mechanics Department, IIT Delhi
% e-mail: divyaprakash.poddar@gmail.com
% Date  : 14 January 2021
    
    X = [ci,cj];
    % Get cylinder and wall coordinates
    [isn, icn] = updateGrid(nx,ny,ci,cj,a,b,theta);

    % LBM Parameters (D2Q9)
    na      = 9;                                    
    ex      = [0 1 0 -1 0 1 -1 -1 1];
    ey      = [0 0 1 0 -1 1 1 -1 -1];

    % Initialize Cylinder Boundary Coordinates
    xb = 0; 
    yb = 0;
    
    % Get boundary node locations
    for i = 3:nx-3 % Limits so as to avoid the walls
        for j = 3:ny-3 % Limits so as to avoid the walls
            if (isn(i,j) == 0)
                for a = 1:na
                    ia = i + ex(a);
                    ja = j + ey(a);
                    if (isn(ia,ja) == 1) % It means a boundary is encountered
                        cx = i + 1/2*ex(a);
                        cy = j + 1/2*ey(a);
                        xb = [xb cx];
                        yb = [yb cy];
                    end
                end
            end
        end
    end

    xb = xb(2:end);
    yb = yb(2:end);

    % Delete the multiple boundary coordinate entries
    idD = 0;
    for i = 1:numel(xb) % Number of total iterations (Equal to original array size)
        if i <= numel(xb) % To avoid invalid index error as the array reduces in size
            p1 = xb(i); p2 = yb(i);
            idD = getSimilarPairsId(i,p1,p2,xb,yb);
            idD = idD(2:end);
            xb(idD) = [];
            yb(idD) = [];
        end
    end
    
    % Number of boundary nodes
    nb = numel(xb);
end
