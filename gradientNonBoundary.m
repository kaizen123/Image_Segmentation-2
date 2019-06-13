function DVDX = gradientNonBoundary(X,XPQTuples)
%  DVDX = gradientNonBoundary(X,XPQTuples)
%
%   Outputs a cell array where each row corresponds to a center in X, used 
%   to generate the voronoi diagram, and the entries in the row are
%   gradient matrices.
%
% INPUTS:
%
%   X: an array of coordinates having size (N,2) that will be used
%      to generate voronoi centers where N is the number of coordinates.
%
%   XPQTuples: a cell array of size (N,Z), where N is the number of centers
%              used to generate the voronoi diagram and Z is the max number
%              of interior vertices  dependent on any one center. Each row
%              is made up of arrays of size (1,4) whose entries correspond
%              to indices into X and intV. The first three indices in this
%              array correspond to the three points in X that generate a
%              vertex and the fourth index corresponds to this vertex in
%              intV.
%
%   NOTE: also see help for:
%       XPQTuples = getXPQTuples(xGenV,vDepX)
%       [intV,vDepX,xGenV] = intVDependsOnX(V3,C,bounds)
%
% OUTPUTS:
%
%   DVDX: a cell array the same size as XPQTuples. The ith row will contain
%   the gradient matrices for each vertex that is dependent on center i.
%
%   e.g.
%       if xGenV{i,:} = [5,7] and 
%
%       DVDX{i,:} = [0.3521, 0.4610; 0.2142, 0.2805],
%                   [0.4231, 0.7843; 0.1231, 0.5643],
%                   [0,0;0,0],...,[0,0;0,0] then,
%
%       [0.3521, 0.4610; 0.2142, 0.2805] is the gradient matrix for
%       intV(5,:) with respect to X(i,:) and, 
%
%       [0.4231, 0.7843; 0.1231, 0.5643] is the gradient matrix for 
%       intV(7,:) with respect to X(i,:).
%       
%       in general each matrix has the following layout:
%       [ intV(z,1)//X(i,1), intV(z,1)//X(i,2) ]
%       [ intV(z,2)//X(i,1), intV(z,2)//X(i,2) ]
%       where // denotes the partial derivative
%
%   NOTE: some of the entries in DVDX will be zero matrices becuase of the
%   way XPQTuples was constructed (see note in getXPQTuples). 


DVDX = cellfun(@(z)gradV(X,z),XPQTuples,'UniformOutput',false);

end



function gradientV = gradV(Centers,XPQCELL)
if isempty(XPQCELL)
    gradientV = zeros(2,2);
else
    ind_x = XPQCELL(1);
    ind_p = XPQCELL(2);
    ind_q = XPQCELL(3);
    
    X = Centers(ind_x,:);
    P = Centers(ind_p,:);
    Q = Centers(ind_q,:);
    
    % get  vectors P - X1 & perp
    p_x = P - X;
    p_x = p_x /norm(p_x);
    p_x_perp = [p_x(2),-p_x(1)];
    
    % get  vectors Q - X1 & perp
    q_x = Q - X;
    q_x = q_x /norm(q_x);
    q_x_perp = [q_x(2),-q_x(1)];
    
    % get the vectors Q - P & perp
    q_p = Q - P;
    q_p = q_p/norm(q_p);
    q_p_perp = [q_p(2),-q_p(1)];
    
    T = [p_x',q_x'];
    
    dv1 = 1/(2 * p_x * q_p_perp'); % change in v corresponding to basis element 1
    dv2 = 1/(2 * q_p_perp * q_x'); % change in v corresponding to basis element 2
    DV1 = dv1*q_p_perp;
    DV2 = dv2*q_p_perp;
    CBDV = [DV1',DV2'];
    gradientV = CBDV*inv(T);
end
end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The double comment lines of code below are how we calculated the change
    % in V given a dx
    
    % % form the change of basis matrix. This will be used to get the
    % % perturbation of x (DX) in terms of the directions of:
    % % (i) p-x and (ii) (q-x).
    % T = inv([p_x',q_x']);
    % CB_DX = T*dx
    % % find the change in vertex v along the direction of q_p_perp as a result
    % % of moving x along p_x and then update q_x
    % dv1 = CB_DX(1)/(2 * p_x * q_p_perp'); % change in v corresponding to basis element 1
    % q_x = x(2,:) - x(1,:) - CB_DX(1)*p_x;
    % q_x = q_x /norm(q_x);
    % q_x_perp = [q_x(2),-q_x(1)];
    % dv2 = CB_DX(2)/(2 * q_p_perp * q_x'); % change in v corresponding to basis element 2
    % % calculate the new vertex & plot with a red square
    % v3_dv = [V1(1,3), V2(1,3)]+(dv1 + dv2)*q_p_perp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

