function XPQTuples = getXPQTuples(xGenV,vDepX)
%   XPQTuples = getXPQTuples(xGenV,vDepX)
%
%   This function helps with identifying which three centers in the voronoi
%   diagram generate a specific interior vertex. Useful in computing
%   gradients of vertices.
%
%
% INPUTS:
%   
%   xGenV: a cell array of size (N,1), where N is the number of centers
%           used to generate the voronoi diagram. Each row in xGenV is a
%           list of indices for intV
%   vDepX: an array of size (K,3), where K is the number of interior vertices in V3,
%          whose rows are the indices of centers in X that generate a
%          vertex.
%
% OUTPUTS:
%
%   XPQTuples: a cell array of size (N,Z), where N is the number of centers
%              used to generate the voronoi diagram and Z is the max number
%              of interior vertices  dependent on any one center. Each row
%              is made up of arrays of size (1,4) whose entries correspond
%              to indices into X and intV. The first three indices in this
%              array correspond to the three points in X that generate a
%              vertex and the fourth index corresponds to this vertex in
%              intV.
%   e.g.
%       if XPQTuples(1,:) = {[1,4,6,23],[1,2,9,54],[1,3,7,34],[],[],[]}
%       then:
%       X(1,:),X(4,:), and X(6,:) generate intV(23,:)
%       X(1,:),X(2,:), and X(9,:) generate intV(54,:)
%       X(1,:),X(3,:), and X(7,:) generate intV(34,:)
%   
%   NOTE: XPQTuples was constructed as a dynamic cell, therefore some of
%   the entries may be empty. These empty cells get converted to zero
%   matrices when passed through the gradient calculation.
%
%
%   NOTE: Also see help for:
%   [intV,vDepX,xGenV] = intVDependsOnX(V3,C,bounds)
%   [V3,C] = voronoin(X)


sz = size(xGenV,1);

XPQTuples = {};
for i = 1:sz
    for j = 1:size(xGenV{i},1)
        vertex = xGenV{i}(j); % get jth vertex that xi generates
        cindx = vDepX(vertex,:)~=i; % find the indices for centers not equal to i
                              % that also generate vertex.
        centers = vDepX(vertex,cindx);
        XPQTuples{i,j} = [i,centers,vertex];
    end
end





