function gradient = chainRule(xGenV, xGenB, DVDX, DBDX, grad, dictV, dictB)
% gradient = chainRule(DVDX, DBDX, xGenV, xGenB, gradV, gradB)
%
% INPUTS:
%   xGenV: N x 1. The ith row is the list of inside vertices related to center i.
%   xGenB: N x 1. The ith row is the list of boundary vertices related to center i.
%
%   DVDX: N x S. 
%   DVDX{i,j}=[V(z,1)//X(i,1),V(z,1)//X(i,2);V(z,2)//X(i,1),V(z,2)//X(i,2)]
%   where xGenV{i, j} = z. 
%   DBDX: N x T. 
%   DBDX{i,j}=[B(z,1)//X(i,1),B(z,1)//X(i,2);B(z,2)//X(i,1),B(z,2)//X(i,2)]
%   where xGenB{i, j} = z. 
%
%   grad: 2 x 2 x E. 
%   grad(d,pt,e)=(energy//edge(i))*(edge(e)//Vd(pt,e)). where Vd refers to V1 or V2

n = size(xGenV, 1);
gradientV = zeros(n, 2);

for x = 1:n
    DxV = [0 0];
    
    
    vs = xGenV{x}; % get the indices for interior vertices influenced by this x
    if isempty(vs)
       continue; 
    end
    m = size(vs, 1); % num of interior vertices
    for j = 1:m % j is the index for interior vertex
        dvdx = DVDX{x,j};  % 2 x 2
        location = dictV(vs(j),:); % row num, col num in V1,V2
        l = size(location, 2); % number of edges x influences
        for i = 1:l
            info = location{i};
            if isempty(info)
                break
            end
            pt = info(1); edge = info(2);
        
            djdv = grad(:, pt, edge);  % 2 x 1
            djdx = dvdx.*djdv;
            DxV = DxV + sum(djdx, 1);
        end
    end
    gradientV(x,:) = DxV;
end

n = size(xGenB, 1); 
gradientB = zeros(n,2);

for x = 1:n
    DxB = [0,0];

    bs = xGenB{x};
    if isempty(bs)
       continue; 
    end
    m = size(bs, 2);
    for j = 1:m
        dbdx = DBDX{x,j};  % 2 x 2
        location = dictB(bs(j),:);
        l = size(location, 2);
        % skip empty
        for i = 1:l
            info = location{i};
            if isempty(info)
                break
            end
            pt = info(1); edge = info(2);
        
            djdb = grad(:, pt, edge);  % 2 x 1
            djdx = dbdx.*djdb; 
            DxB = DxB + sum(djdx, 1);

            
        end
    end
    
    gradientB(x, :) = DxB;
    
end

gradient = gradientB + gradientV;

