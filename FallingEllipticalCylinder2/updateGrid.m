function [isn, icn] = updateGrid(nx,ny,ci,cj,a,b,theta)
% [isn, icn] = updateGrid(nx,ny,ci,cj,a,b,theta)
%   Generates matrices containing flags for wall, solid and fluid nodes
%
% input:
%   nx          = Grid size in x-direction
%   ny          = Grid size in y-direction
%   ci          = Ellipse center (x-coordinate)
%   cj          = Ellipse center (y-coordinate)
%   a           = Semi-major axis
%   b           = Semi-minor axis
%   theta       = Angle by which the ellipse has been rotated
% output:
%   isn         = A matrix containing flags (1 for wall and cylinder and 0 for fluid)
%   icn         = A matrix containing flags (1 for cylinder and 0 for fluid)
%
% Author: Divyaprakash (2020AMZ8461), PhD Candidate
%         Applied Mechanics Department, IIT Delhi
% e-mail: divyaprakash.poddar@gmail.com
% Date  : 14 January 2021

    icn     = zeros(nx,ny); % [0--> fluid, 1--> Fictitious Fluid]
    isn     = ones(nx,ny);  % [0--> Fluid, 1--> Solid]
    isn(:,2:end-1) = 0;
    d       = @(x,y,t) ((x-ci)*cos(t) + (y-cj)*sin(t))^2/a^2 + ((x-ci)*sin(t) - (y-cj)*cos(t))^2/b^2;
    for i = 1:nx
        for j = 1:ny
            if d(i,j,theta) <= 1
                isn(i,j) = 1; 
                icn(i,j) = 1; % Cylinder (Fictitious Fluid) Nodes
            end
        end
    end
end
