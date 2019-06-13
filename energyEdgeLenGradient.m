function [energy, grad] = energyEdgeLenGradient(V1, V2)
% [energy, grad] = energyEdgeLenGradient(V1, V2)
%
% INPUTS:
%   V1: an array of size (2,M) of x coordinates for the edges in the voronoi
%       diagram. V1(1,:) is connected to V1(2,:)
%
%   V2: an array of size (2,M) of y coordinates for the edges in the voronoi
%       diagram. V2(1,:) is connected to V2(2,:)
%
% OUTPUTS:
%   energy: scalar, the square sum of the edges. 
%   sum_{i.j}[(xi-xj)^2 + (yi-yj)^2] where i,j indicate some joined vertices
%
%   grad:an array of size (2,2,M). Record some partial derivative. 
%   Denote L(e)=(xi-xj)^2 + (yi-yj)^2 where e is the edge index,  
%         i.e. xi = V1(1,e), xj=V1(2,e), yi = V2(1,e), yj=V2(2,e).
%   energy = sum_{e} L(e)
%   grad(d,p,e) = (energy//L(e)) * (L(e)//Vd(p,e))
%               =              1 * 2(Vd(p,e)-Vd(p',e))

diff1 = V1(1,:) - V1(2,:);
diff2 = V2(1,:) - V2(2,:);
n = size(diff1, 2);

len = sqrt(diff1.*diff1 + diff2.*diff2);
energy = sum(len);

grad = zeros(2,2,n);
grad(1,:,:) = [diff1; -diff1];
grad(2,:,:) = [diff2; -diff2];
grad = grad ./ reshape(len, 1,1,n);
end