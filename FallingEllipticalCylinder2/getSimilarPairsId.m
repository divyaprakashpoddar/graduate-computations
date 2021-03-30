function id = getSimilarPairsId(i,p1,p2,CX,CY)
% id = getSimilarPairsId(i,p1,p2,CX,CY)
%   Returns the index of the duplicates of a coordinate pair (p1,p2)
%
% input:
%   i           = Index of the coordinate which is compared
%   p1          = A boundary node's x-coordinate which is compared
%   p2          = A boundary node's y-coordinate which is compared
%   CX          = Vector of boundary node's x-ccordinate
%   CY          = Vector of boundary node's y-ccordinate
% output:
%   id          = Index of the duplicates of coordinate pair (p1,p2)
%
% Author: Divyaprakash (2020AMZ8461), PhD Candidate
%         Applied Mechanics Department, IIT Delhi
% e-mail: divyaprakash.poddar@gmail.com
% Date  : 14 January 2021

    id = 0;
    for j = 1:numel(CY)
        if j ~= i
            if (abs(p1-CX(j)) < 1e-5 && abs(p2-CY(j)) < 1e-5)
            %if (p1 == CX(j) && p2 == CY(j))    
                id = [id j];
            end
        end
    end
end
