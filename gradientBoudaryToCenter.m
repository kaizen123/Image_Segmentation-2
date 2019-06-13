function [DBDX, xGenB] = gradientBoudaryToCenter(X, B, BX)
% DBX = gradientBoudaryToCenter(X, B, BX)
%
%   This function finds the gradient of the boundary points with respect of
%   the coordinates of the two generating points (centers). 
%
% INPUTS:
%
%   X: an array of coordinates having size (N,2), where N is the number
%      of coordinates in X, that was used to generate the voronoi diagram.
%
%   B: an array of size (b,3), where b is the number of boundary
%       coordinates. B(:,1) is an array containing the types, T, of
%       the boundary vertices. B(:,2:3) is an array containg the
%       coordinates of the boundary vertices.
%
%   BX: an array of size (b,2), where b is the number of boundary
%       coordinates, and the elements in each row correspond to the two
%       indices for points in X that generate the perpindicular bisector
%       that intersects the voronoi boundary at the points in B.
%
%   NOTE: [V1,V2,B] = findBounds(V1,V2,bounds), BX = boundaryCenters(B(:,2:3))
%          and [V1 ,V2] = voronoi(X(:,1),X(:,2)).
%
% OUTPUT:
%
%   DBX<no longer return>: an array of size (b,8), where b is the number of 
%       boundary coordinates, and elements in each row correspond to the
%       gradients of two coordinates of the boundary point with respect of
%       the four coordinates of the two centers. For a given B(i,:) = (b1, 
%       b2), BX(i,:) is [j k] and then the corresponding centers are X(j,:)
%       and X(k,:). DBX(i,:) = [b1//X(j,1), b1//X(j,2), b1//X(k,1),
%       b1//X(k,2), b2//X(j,1), b2//X(j,2), b2//X(k,1), b2//X(k,2)], where
%       // is the partial derivative. In most cases, half of the row is
%       zero because either b1 or b2 is a constant. 
%
%   DBDX: N x b. 
%   DBDX{i,j}=[B(z,1)//X(i,1),B(z,1)//X(i,2);B(z,2)//X(i,1),B(z,2)//X(i,2)]
%   where xGenB{i, j} = z. 
%
%   xGenB: cell array of N x 1, ith row is a list of indices of vertices
%   related to the ith polygen center


% b is the number of boundary points for which we wish to find centers for
b = size(B,1);
% n is the number of points in X to search through
n = size(X,1);

% Types code for boundary points
TYPE_ON_X_AXIS = 1; % Constant value for the 1st coordinate
TYPE_ON_Y_AXIS = 2; % Constant value for the 2nd coordinate

% xj_diff =   [(X(j,1)-b1), (X(j,2)-b2)]
% xk_diff = - [(X(k,1)-b1), (X(k,2)-b2)]
% div_2 = (X(j,2)-X(k,2))
% div_1 = (X(j,1)-X(k,1))
% For b1 is constant
% DBX(i,5) = b2//X(j,1)=   (X(j,1)-b1) / (X(j,2)-X(k,2))
% DBX(i,6) = b2//X(j,2)=   (X(j,2)-b2) / (X(j,2)-X(k,2))
%                     =>    xj_diff    / div_2
% DBX(i,7) = b2//X(k,1)= - (X(k,1)-b1) / (X(j,2)-X(k,2))
% DBX(i,8) = b2//X(k,2)= - (X(k,2)-b2) / (X(j,2)-X(k,2))
%                     =>    xk_diff    / div_2
%
% For b2 is constant
% DBX(i,1) = b1//X(j,1)=   (X(j,1)-b1) / (X(j,1)-X(k,1))
% DBX(i,2) = b1//X(j,2)=   (X(j,2)-b2) / (X(j,1)-X(k,1))
%                     =>    xj_diff    / div_1
% DBX(i,3) = b1//X(k,1)= - (X(k,1)-b1) / (X(j,1)-X(k,1))
% DBX(i,4) = b1//X(k,2)= - (X(k,2)-b2) / (X(j,1)-X(k,1))
%                     =>    xk_diff    / div_1

div = zeros(b,1);
% For b2 is constant
ID_1 = B(:,1)==TYPE_ON_Y_AXIS; % id for boundary points with constant b2
div(ID_1) = X(BX(ID_1, 1), 1) - X(BX(ID_1, 2), 1); % div1
%     X(    j    , 1) - X(    k    , 1)
% For b1 is constant
ID_2 = B(:,1)==TYPE_ON_X_AXIS; % id for boundary points with constant b1
div(ID_2) = X(BX(ID_2, 1), 2) - X(BX(ID_2, 2), 2); % div2
%     X(    j    , 2) - X(    k    , 2)

DBX = zeros(b,8);
DBX(:, 1:2) =   X(BX(:, 1), :) - B(:, 2:3); % xj_diff
%               X(j       , 1) - b1 
%               X(j       , 2) - b2
DBX(:, 3:4) = - X(BX(:, 2), :) + B(:, 2:3); % xk_diff
%             - X(k       , 1) + b1 
%             - X(k       , 2) + b2
DBX = DBX./div;
DBX(ID_2,5:8) = DBX(ID_2,1:4);
DBX(ID_2,1:4) = 0;

%%%%%%%%%%%%%%

DBDX = cell(n,b);
xGenB = cell(n, 1);

for i = 1 : b
    row = DBX(i, :);
    for j = 1:2
        x = BX(i, j);
        xGenB{x} = [xGenB{x}, i];
        len = size(xGenB{x}, 2);
        DBDX{x, len, :} = [row(j*2-1:j*2);row(j*2+3:j*2+4)];
    end
end
end


