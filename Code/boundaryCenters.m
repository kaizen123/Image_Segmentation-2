function BX = boundaryCenters(X,B)
% BX = boundaryCenters(X,B)
%
%   This function finds the two generating points (centers) whose
%   perpindicular bisector intersects the voronoi diagram at the points
%   given in B. This function finds the distances between all combinations
%   of points in X and B. For a given B(i,:), this function searches for
%   X(j,:) and X(k,:) such that |B(i,:) - X(j,:)| == |B(i,:) - X(k,:)|.
%
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
% OUTPUTS:
%
%   BX: an array of size (b,2), where b is the number of boundary
%       coordinates, and the elements in each row correspond to the two
%       indices for points in X that generate the perpindicular bisector
%       that intersects the voronoi boundary at the points in B.
%
%   e.g. BX(i,:) = [j,k] such that the vector perpindicular to
%        [ X(j,:) - X(k,:) ] intersects the voronoi boundary at B(i,:)
%
%   NOTE: Also see help for findBounds. [V1,V2,B] = findBounds(V1,V2,bounds)
%

% n is the number of boundary points for which we wish to find centers for
n = size(B,1);

% m is the number of points in X to search through
m = size(X,1);

% The rows in B_mat(:,:,1) are copies of [B(1,1), B(2,1), ..., B(n,1)]
% i.e. x coordinates of the points in B
% The rows in B_mat(:,:,2) are copies of [B(1,2), B(2,2), ..., B(n,2)]
% i.e. y coordinates of the points in B
B_mat = zeros(m,n,2);
B_mat(:,:,1) = ones(m,1)*B(:,1)';
B_mat(:,:,2) = ones(m,1)*B(:,2)';

% The columns in x_mat(:,:,1) are copies of [X(1,1), X(2,1), ..., X(m,1)]'
% i.e. x coordinates of the points in X
% The columns in x_mat(:,:,2) are copies of [X(1,2), X(2,2), ..., X(m,2)]'
% i.e. y coordinates of the points in X
x_mat = zeros(m,n,2);
x_mat(:,:,1) = X(:,1)*ones(1,n);
x_mat(:,:,2) = X(:,2)*ones(1,n);

% Compute the elementwise (along dim = 3) 2 norm
diff_mat = B_mat - x_mat;
s_mat = diff_mat.^2;
mat = sum(s_mat,3);
mat = mat.^(1/2);

vals = round(mat,8); % accuracy is 1e-8 and can be adjusted if needed.
vals1 = vals; % vectorize
vals2 = vals;


%%%% THIS SECTION HAS A BUG. FIXED SECTION BELOW
% Since exactly two points generate a perpindicular bisector that
% intersects at the boundary, then there should be exactly two elements in
% each column of vals that are equal. Find the unique elements in vals and
% set these elements equal to 0.
% NOTE: unique(A) will return first occurences of unique elements in A.
% NOTE: unique(A,'last') will return last occurences of unique elements in A.
% [~,IA1,~]= unique(vals1);
% vals1(IA1) = 0; % leaves out last occurences of repeated elements.
% [~,IA2,~]= unique(vals2,'last'); % leaves out first occurences
% vals2(IA2) = 0;


% Since exactly two points generate a perpindicular bisector that
% intersects at the boundary, then there should be exactly two elements in
% each column of vals that are equal. Find the unique elements in vals and
% set these elements equal to 0.
% NOTE: unique(A) will return first occurences of unique elements in A.
% NOTE: unique(A,'last') will return last occurences of unique elements in A.
for i = 1:size(vals1,2)
    [~,IA1,~]= unique(vals1(:,i));
    vals1(IA1,i) = 0; % leaves out last occurences of repeated elements.
    [~,IA2,~]= unique(vals2(:,i),'last'); % leaves out first occurences
    vals2(IA2,i) = 0;
end





% combine results above resulting in both first and last occurences of
% repeated elements along with their linear indices.
vals = vals1 + vals2;
vals = reshape(vals,m,n);
indices = find(vals ~=0);
[I,~] = ind2sub([m,n],indices);
sz = size(I,1);
BX = reshape(I,2,sz/2)';
end

