function [energy, gradient] = energyDiffGradient(xGenV, xGenB, V, B, X, DVDX, DBDX, VX, BX, pic, delta)
% % energyDiffGradient(xGenV, xGenB, DVDX, DBDX, grad1, grad2, dictV, dictB, E2C)
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

%%%%%%%%%%%%%%integration function to check. 

n = size(xGenV, 1);
energy = 0;
gradient = zeros(n, 2);

% Get the grids ready for inpolygon to find the inside pixels 
[W, H] = size(pic);
grids = zeros(H*W, 2);
for j = 0: (H-1)
    for i = 0:(W-1)
        grids(j*W+i+1, :) = [i, j];
    end
end

% Decide which x_center is closest to every corner
corners = [0, 0; W-1, 0; 0, H-1; W-1, H-1];
corner2x = [];
ncorner = size(corners, 1);
for c = 1:ncorner
    dis = sum((X-corners(c, :)).^2, 2);
    [~, x] = min(dis);
    % corner2x(i) is the idx of x for corners(i)
    corner2x = [corner2x, x];
end 

% CHECK IF THE CORNERS ARE CORRECT
% for i = 1:size(corner2x,2)
%    myVoronoi(X,[0,10],'p');
%    hold on;
%    plot(corners(i,1),corners(i,2),'rx','MarkerSize',10)
%    hold on;
%    xindx = corner2x(i);
%    plot(X(xindx,1),X(xindx,2),'rx','MarkerSize',10)
%    clf
% end


cvalues = zeros(n, 1);
coordinates_cache = cell(n, 1);
deriIdx_cache = cell(n, 1);
coorIdx_cache = cell(n, 1);

for x = 1:n
    % order the vertex, anti-clockwise
    % idx is the index of every ordered vertex, which is used for DVDX,DBDX
    % E.G. idx = [2, -2, 0, 0, -1, 1] means the vertices are ordered as:
    % xGenV{x,2}, xGenB{x,2}, a corner, a corner, xGenB{x,1}, xGenV{x,1}
    [deriIdx, coorIdx, coordinates] = sortV(transpose(xGenV{x}), xGenB{x}, find(corner2x==x),V, B, corners, X(x, :));
    coordinates_cache{x} = coordinates;
    deriIdx_cache{x} = deriIdx;
    coorIdx_cache{x} = coorIdx;
    [c, loss] = polygonLoss(coordinates, pic, grids);
    cvalues(x) = c;
    energy = energy + loss;
    
        % CHECK IF VERTICES ARE IN THE CORRECT ORDER FOR EACH POLYGON
%     myVoronoi(X,[0,10],'p');
%     hold on
%     plot(X(x,1),X(x,2),'rx','MarkerSize',10)  
%     hold on
%     for v = 1:size(coordinates,1)
%        plot(coordinates(v,1),coordinates(v,2),'rx','MarkerSize',10)
%        hold on;
%        text(coordinates(v,1),coordinates(v,2),num2str(v),'FontSize',15);
%        hold on
%     end
%     clf
end



for x = 1:n
    coordinates = coordinates_cache{x};
    deriIdx = deriIdx_cache{x};
    coorIdx = coorIdx_cache{x};
    
    m = size(deriIdx, 2);
    c = cvalues(x);
    
    for i = 1:m
        % The idx for the second point of this edge
        i_ = mod(i, m) + 1;
        % Get the coordinates of the edge
        X1 = coordinates(i, 1);  Y1 = coordinates(i, 2);
        X2 = coordinates(i_, 1);  Y2 = coordinates(i_, 2);
        
        % Get the partial derivative dv/dx for the first point
        deriIdx1 = deriIdx(i);
        coorIdx1 = coorIdx(i);
        if deriIdx1 > 0
            dvdx1 = DVDX{x, deriIdx1};
            xs1 = VX(coorIdx1,:);
        elseif deriIdx1 < 0
            dvdx1 = DBDX{x, -deriIdx1};
            xs1 = BX(-coorIdx1,:);
        else
            dvdx1 = [0 0; 0 0];
            xs1 = x;
        end
        
        % Get the partial derivative dv/dx for the second point
        deriIdx2 = deriIdx(i_);
        coorIdx2 = coorIdx(i_);
        if deriIdx2 > 0
            dvdx2 = DVDX{x, deriIdx2};
            xs2 = VX(coorIdx2,:);
        elseif deriIdx2 < 0
            dvdx2 = DBDX{x, -deriIdx2};
            xs2 = BX(-coorIdx2,:);
        else
            dvdx2 = [0 0; 0 0];
            xs2 = x;
        end
        xs = intersect(xs1, xs2);
        neighbor = xs(xs~=x);
        %endpt = [X1, Y1; X2, Y2];
        
        int = edgeInt(X1, Y1, X2, Y2, X, c, dvdx1, dvdx2, cvalues(neighbor), delta, pic);
        gradient(x, :) = gradient(x, :) + int;
    end

end
%norm = max(sqrt(gradient(:,1).^2 + gradient(:,2).^2));
norm=1000;
gradient = gradient/norm;
energy = energy/norm;

end



function int = edgeInt(X1, Y1, X2, Y2, X, c, dvdx1, dvdx2, cvalues, delta, pic)
    % The length of the edge
    length = sqrt((X1-X2)^2 + (Y1-Y2)^2);
    % The normal vecor pointing outside
    normal = [Y2-Y1; X1-X2]/length;

%     % PLOT NORMAL VECTOR
%     clf;
%     myVoronoi(X,[0,99],'p');
%     hold on;
%     vNormal = [0.5*(X1 + X2),0.5*(Y1 + Y2)] + normal';
%     midpoint = [0.5*(X1 + X2),0.5*(Y1 + Y2)];
%     vEdgeX = [vNormal(1,1),midpoint(1,1)]';
%     vEdgeY = [vNormal(1,2),midpoint(1,2)]';
%     plot(vEdgeX,vEdgeY,'k-o')
%     hold on
%     text(vEdgeX(2),vEdgeY(2),'Current Normal Vector', 'FontSize',15)

    %See the handwritten
    % Divide the edge into N pieces
    N = fix(length/delta)+1;
    delta_x=(X2-X1)/N; delta_y=(Y2-Y1)/N;

    % <v_a ` normal>, <v_b ` normal>
    VEND1 = sum(dvdx1 .* normal,1); %1 x 2
    VEND2 = sum(dvdx2 .* normal,1);
    temp = zeros(1, 2);
    for j = 1:N
        % [X1, Y1, delta_x, delta_y, j]
        VHERE = (1 - (j-0.5)/N)*VEND1 + ((j-0.5)/N)*VEND2;
        f = interpolation(X1+(j-0.5)*delta_x, Y1+(j-0.5)*delta_y, pic);
        value = (c - f)^2;
        if size(cvalues, 2) == 1
           value =  value - (cvalues - f)^2;
        end
        temp = temp + value*VHERE;
    end
    int = temp*length/N;
end

% sortV(transpose(xGenV{x}), xGenB{x}, find(corner2x==x),V, B, corners, X(x, :));
function [deriIdxSort, coorIdxSort, coordinates] = sortV(vs, bs, cs, V, B, corners, center)
    deriIdxSort = [(1:size(vs, 2)), -(1:size(bs, 2)), zeros(size(cs))];
    coorIdxSort = [vs, (-bs), zeros(size(cs))];
    coordinates = [V(vs,:);B(bs,2:3);corners(cs, :)];
    ang = atan2(coordinates(:, 2) - center(2), coordinates(:, 1) - center(1));
    [~, order] = sort(ang);
    deriIdxSort = deriIdxSort(order);
    coorIdxSort = coorIdxSort(order);
    coordinates = coordinates(order, :);
end

function [avg, loss] = polygonLoss(coordinates, pic, grids)
    in = inpolygon(grids(:,1),grids(:,2),coordinates(:,1),coordinates(:,2));
    pts = grids(in,:)+1;
    pixel = [];
    num = size(pts,1);
    for i = 1:num
        pt = pts(i, :);
        pixel = [pixel; pic(pt(1), pt(2))];
    end
    avg = sum(pixel)/num;
    loss = [];
    for i = 1:num
        loss = [loss; (avg-pixel(i))^2];
    end
    loss = sum(loss);
end

function value = interpolation(x, y, pic)
try
    value = interpolation_core(x, y, pic);
catch exception
   value = interpolation_core(x, y, pic);
end
end

function value = interpolation_core(x, y, pic)
    % t = [x,y]
    x_ = fix(x); y_ = fix(y); 
    dx = x - x_; dy = y - y_;
    x0 = int16(x_+1); y0 = int16(y_+1);
    if equal(x_, x)
        if equal(y_, y)
            value = pic(x0, y0);
        else
            p0 = pic(x0, y0); p1 = pic(x0, y0+1);
            value = (1-dy)*p0 + dy*p1;
        end
    else
        if equal(y_, y)
            p0 = pic(x0, y0); p1 = pic(x0+1, y0);
            value = (1-dx)*p0 + dx*p1;
        else
            p00 = pic(x0, y0); p01 = pic(x0, y0+1);
            p10 = pic(x0+1, y0); p11 = pic(x0+1, y0+1);
            value = (1-dx)*(1-dy)*p00 + (1-dx)*dy*p01 + dx*(1-dy)*p10 + dx*dy*p11;
        end
    end
end

function cmp = equal(x, y)
cmp = abs(x-y)<1e-7;
end
