function [IND] = getDependentVertices(X,V1,V2,V3,C)
%   IND = getDependentVertices(X,V1,V2,V3,C)
% 
%   This function will output a cell array IND, where IND{i,j} contains
%   the indices in V1(j,:) & V2(j,:) that are dependent on X(i,:).
% 
% INPUTS:
%
%   X: an array of coordinates having size (N,2) that will be used
%      to generate voronoi centers where N is the number of coordinates.
%
%   V1: an array of size (2,M) of x coordinates for the edges in the voronoi
%       diagram. V1(1,:) is connected to V1(2,:)
%
%   V2: an array of size (2,M) of y coordinates for the edges in the voronoi
%       diagram. V2(1,:) is connected to V2(2,:)
%
%   V3: an array of size (M,2) of coordinates of vetices in the voronoi diagram
%       used by C.
%
%   C:  a cell array that holds the vertices that are dependent on the
%       points in X.
%       i.e. C{i} = [3,5,7] => [V3(1,3),V3(1,3)],[V3(1,5),V3(1,5)],[V3(1,7),V3(1,7)]
%       are dependent on X(i,:).
%
%   NOTE: [V1,V2] = voronoi(X(:,1),X(:,2)) & [V3,C] = voronoin(X)
%
% OUTPUT:
%
%   IND: a cell array of size (N,2), where N is the number of points in X,
%   where the ith row contains the indices for V1 & V2 which are the
%   vertices that are dependent on the ith point in X.
%
% e.g.:
%
%   If X is of size (5,2):
%
%   IND{3,1} = [ind1, ind2, ...] where [V1(3,ind1),V2(3,ind1)],
%   [V1(1,ind2),V2(1,ind2)], [V1(1,ind3),V2(1,ind3)], ..., are dependent on
%   the point X(3,:).
%
%   IND{3,2} = [ind1, ind2, ...] where [V1(2,ind1),V2(2,ind1)],
%   [V1(2,ind2),V2(2,ind2)], [V1(2,ind3),V2(2,ind3)], ..., are dependent on
%   the point X(3,:).
%   
%   The reason for two different cells:
%   for a given column index of V1,V2, (say 3), the points
%   [V1(1,3),V2(1,3)] & [V1(2,3),V2(2,3)] are connected by an edge.
%

szx = size(X,1); % number of point in X


% All of the interior vertices are contained in the first rows
% of V1 & V2. Extract those to search for dependencies.

intVx1 = V1(1,:); % interior vertices: x coordinates (row)
intVy1 = V2(1,:); % interior vertices: y coordinates (row)

intVx2 = V1(2,:); % interior vertices: x coordinates (row)
intVy2 = V2(2,:); % interior vertices: y coordinates (row)
IND = {szx,2};

for i = 1:szx
    C2 =  C{i}(C{i}~=1); % remove the 1 indicating an inf vertex
    
    kvx = V3(C2,1); % get the known x coordinates for vertices depend. on X
    kvy = V3(C2,2); % get the known y coordinates for vertices depend. on X
    
    % The following commented lines are left in to make sure optimized
    % version works well.
%     [indx_x1,~] = ismembertol(intVx1,kvx,1e-6); % get logical array of intersecting elem.
%     indx_x1 = find(indx_x1 == 1); % get the indices where intersection is true.
%     
%     [indx_x2,~] = ismembertol(intVx2,kvx,1e-6); % get logical array of intersecting elem.
%     indx_x2 = find(indx_x2 == 1); % get the indices where intersection is true.
%     
    
    % FOR OPTIMIZING
    sz_intVx1 = size(intVx1,2);
    kvx_mat = repmat(kvx,1,sz_intVx1); % replicate matrix for logical comparison
    indx_x1 = sum(intVx1 == kvx_mat,1); % get logical array of intersectinf elem.
    %nnz(indx_x1 == indx_x1); for debugging
    indx_x1 = find(indx_x1 == 1); % get indices where intersection is true.
    
    sz_intVx2 = size(intVx2,2);
    kvx_mat = repmat(kvx,1,sz_intVx2);% replicate matrix for logical comparison
    indx_x2 = sum(intVx2 == kvx_mat,1); % get logical array of intersectinf elem.
    %nnz(indx_x1 == indx_x1);for debugging
    indx_x2 = find(indx_x2 == 1);% get indices where intersection is true.
    
    IND{i,1} = indx_x1;
    IND{i,2} = indx_x2;
end


end
