function dictionary = getDictionaryIntoV1(V1,V2,Coords,bounds)
%   dictionary = getDictionaryIntoV1(V1,Coords)
%
%   This function will help us by creating a dictionary that will connect
%   indices in Coords to indices in V1. 
%
% INPUTS:
%
%   V1: an array of size (2,M) of x coordinates for the edges in the voronoi
%   diagram.
%
%   NOTE: V2 can be used in place of V1. (If the x coordinates is found in
%   a column of V1, the corresponding y coordinate can be found in the same
%   column of V2)
%   
%   Coords: an array of size (K,2), where K is the number of vertices that
%   are to be mapped into V1. Each row contains the x and y coordinate of
%   a vertex of the voronoi diagram.
%
%   NOTE: Coords should either be intV or B(:,2:3).
%
% OUTPUTS:
%
%   dictionary: a dictionary of size (K,Z), where K is the number of
%   vertices that are mapped into V1, and Z is the number of edges the ith
%   vertex contributes to in V1.

sz = size(Coords,1);
V1 = round(V1,4);
V2 = round(V2,4);
Coords = round(Coords,4);


Coords = Coords';
dict = reshape(Coords,2,1,sz);
% compare x coordinates
temp1 = dict(1,:,:) == V1(1,:) & dict(1,:,:) > bounds(1) & dict(1,:,:) < bounds(2);
temp2 = dict(1,:,:) == V1(2,:) & dict(1,:,:) > bounds(1) & dict(1,:,:) < bounds(2);

% compare y coordinates
temp3 = dict(2,:,:) == V2(1,:) & dict(2,:,:) > bounds(1) & dict(2,:,:) < bounds(2);
temp4 = dict(2,:,:) == V2(2,:) & dict(2,:,:) > bounds(1) & dict(2,:,:) < bounds(2);




dictTop = temp1+temp3;
dictBottom = temp2+temp4;

dictionary = cell(sz,3);

for i = 1:sz
    [~,topCols] = find(dictTop(:,:,i) ~= 0);
    [~,bottomCols] = find(dictBottom(:,:,i) ~= 0);
    szt = length(topCols);
    szb = length(bottomCols);
    
    if szt ~= 0 
       for j = 1:szt
           dictionary{i,j} = [1,topCols(j)];
       end
    end
    
    if szb ~= 0 
        start = szt + 1;
        stop = szt + szb;
        indx = 1;
       for j = start:stop
           dictionary{i,j} = [2,bottomCols(indx)];
           indx = indx + 1;
       end
    end
    
    

end

end


