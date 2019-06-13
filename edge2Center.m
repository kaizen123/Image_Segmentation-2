function E2C = edge2Center(V1,vDepX,D1,BX,D2)
% E2C = edge2Center(V1,vDepX,D1,BX,D2)
%
%   This function will look through each edge [V1(i,j),V2(i,j)] and output
%   a list of indices into X. The indices correspond to points which
%   influence the coordinate for the given edge.
%
%
% INPUTS:
%
%   V1: an array of size (2,M) of x coordinates for the edges in the voronoi
%       diagram.
%
%   vDepX: an array of size (K,3), where K is the number of interior vertices in V3,
%          whose rows are the indices of centers in X that generate a
%          vertex.
%
%   D1: a dictionary, for interior vertices of size (K,Z), where K is the number of
%       vertices that are mapped into V1, and Z is the number of edges the ith
%       vertex contributes to in V1.
%
%   BX: an array of size (b,2), where b is the number of boundary
%       coordinates, and the elements in each row correspond to the two
%       indices for points in X that generate the perpindicular bisector
%       that intersects the voronoi boundary at the points in B.
%
%   D2: a dictionary, for boundary vertices of size (K,Z), where K is the number of
%       vertices that are mapped into V1, and Z is the number of edges the ith
%       vertex contributes to in V1.
%
%   NOTE: Also see help for:
%         BX = boundaryCenters(X,B)
%         [intV,vDepX,xGenV] = intVDependsOnX(V3,C,bounds)
%         dictionary = getDictionaryIntoV1(V1,Coords)
%
% OUTPUTS:
%
%   E2C: A cell array of same size as V1. In place of each coordinate is a
%       list of indices into X, where X is the matrix of centers that
%       generate the voronoi diagram. If a point on a given edge is an
%       interior vertex, the corresponding entry in E2C is a triplet. If a
%       point on a given edge is a boundary vertex, the corresponding entry
%       in E2C is a tuple.



E2C = cell(size(V1));

for i = 1:size(D1,1) % ith interior vertex
    for j = 1:size(D1,2) % jth edge where ith vertex appears
     
        edge = D1{i,j}(2);
        row = D1{i,j}(1);
        E2C{row,edge} = [E2C{row,edge},vDepX(i,:)];
    
    end
end

for i = 1:size(D2,1) % ith boundary vertex
     
        edge = D2{i,1}(2);
        row = D2{i,1}(1);
        E2C{row,edge} = [E2C{row,edge},BX(i,:)];
    
end

