% IN SOME RARE CASES, B IS NOT OF SIZE (b,3), INSTEAD IT IS OF SIZE (b,1) or (b,2)
% NEED TO INVESTIGATE


function [V1,V2,B] = findBounds(X,V1,V2,bounds)
%
% [V1,V2,B] = findBounds(V1,V2,bounds)
%
% Generates a square bounded voronoi diagram with the bottom left corner at
% the origin. Matlab will give arbitrary vertices, this function truncates
% the edges in the voronoi diagram along a square boundary specified by
% bounds.
%
% INPUTS:
%
%   V1:     an array of x coordinates, of size (2,N), for the vertices of a voronoi diagram
%
%   V2:     an array of y coordinates, of size (2,N), for the vertices of a voronoi diagram
%
%   bounds: an array of size (1,2) that holds the lower and upper bounds
%           used to generate the boundary
%
%   Note: (N is the number of edges in the voronoi diagram)  
%
%
% OUTPUTS: 
%
%   V1: an array of size (2,N) of x coordinates for the edges in the voronoi
%       diagram. V1(1,:) is connected to V1(2,:)
%
%   V2: an array of size (2,N) of y coordinates for the edges in the voronoi
%       diagram. V2(1,:) is connected to V2(2,:)
%
%   B: an array of size (b,3), where b is the number of boundary
%       coordinates. B(:,1) is an array containing the types, T, of
%       the boundary vertices. B(:,2:3) is an array containg the
%       coordinates of the boundary vertices.
%
%   Note: the ith edge in the voronoi diagram is represented by the vector
%          [V1(1,i) - V1(2,i) , V2(1,i) - V2(2,i)]
xbox = [bounds(1) bounds(1) bounds(2) bounds(2) bounds(1)];
ybox = [bounds(1) bounds(2) bounds(2) bounds(1) bounds(1)];

b_count = 0;
B = [];
% % Extend lines for unbounded cells so they intersect the boundary box
V1n = round(V1,6);

%voronoi(X(:,1),X(:,2));
%grid on
% get the unique x coordinates in row 1 from all edges & count how many
% edges they appear in.
[row1,~,ic] = unique(V1n(1,:)');
counts = accumarray(ic,1);
value_counts1 = [row1, counts];

% get the unique x coordinates in row 2 from all edges & count how many
% edges they appear in.
[row2,~,ic] = unique(V1n(2,:)');
counts = accumarray(ic,1);
value_counts2 = [row2, counts];

% check which x coordinates from row 1 appear in row2
[~,IVC2] = ismember(value_counts1(:,1),value_counts2(:,1));

% get indices elements from row1 who show at least once in row 2
nz_IVC2 = IVC2 ~= 0;
IVC2 = IVC2(nz_IVC2);

% delete those elements that are in row 1 and row 2 from value_counts2
value_counts2(IVC2,:) = [];

% extract x coordinates from the matrix
vertices = value_counts2(:,1);

for index = 1:length(vertices)
    
    %voronoi(X(:,1),X(:,2),'b');

    % get the current x coordinate
    x_coord = vertices(index);
    
    % find the edges where this x coordinate appears
    edge_row1 = find(V1n(1,:) == x_coord);
    edge_row2 = find(V1n(2,:) == x_coord);
    
    % get the x and y coordinates for the endpoints of this edge
    edge_x = V1(:,edge_row2);
    edge_y = V2(:,edge_row2);

    %hold on;
    %plot(edge_x,edge_y,'bo')
    
    % calculate slope line that is to be extended & extend the edge
    m = ( edge_y(1) - edge_y(2) ) / ( edge_x(1) - edge_x(2) );
    s = sign(edge_x(1) - edge_x(2));
    if s == 1
        y2 = edge_y(1) + m*(bounds(1) - edge_x(1));
        V1(2,edge_row2) = bounds(1);
        V2(2,edge_row2) = y2;
        
        %B(b_count+1,1) = determineType([bounds(1),y2],bounds);
        %B(b_count+1,2:3) = [bounds(1),y2];
        %b_count = b_count + 1;
    end
    
    if s == -1
        y2 = edge_y(1) + m*(bounds(2) - edge_x(1));
        V1(2,edge_row2) = bounds(2);
        V2(2,edge_row2) = y2;
        %B(b_count+1,1) = determineType([bounds(2),y2],bounds);
        %B(b_count+1,2:3) = [bounds(2),y2];
        %b_count = b_count + 1;
    end
    
    %hold on;
    %plot(V1,V2,'r--');
    %clf;
end

% Extend edge for line segments with infinite slope
vertical_edges = find(V1(1,:) == V1(2,:));
for index = 1:length(vertical_edges)
    edge = vertical_edges(index);
    y1 = V2(1,edge);
    y2 = V2(2,edge);
    
    y1_in_row1 = sum(V2(1,:) == y1);
    y1_in_row2 = sum(V2(2,:) == y1);
    y2_in_row1 = sum(V2(1,:) == y2);
    y2_in_row2 = sum(V2(2,:) == y2);
    
    y1_appears = y1_in_row1 + y1_in_row2;
    y2_appears = y2_in_row1 + y2_in_row2;
    
    if y1_appears == 1
        s = sign(y1 - y2);
        if s == 1
            V2(1,edge) = bounds(2);
        else
            V2(1,edge) = bounds(1);

        end
    end
    if y2_appears == 1
        s = sign(y2 - y1);
        if s == 1
            V2(2,edge) = bounds(2);
        else
            V2(2,edge) = bounds(1);
        end
    end
end



% get the ith line we wish to truncate
for i = 1:size(V1,2)
    x = [V1(1,i),V1(2,i)]; %(x1,x2)
    y = [V2(1,i),V2(2,i)]; %(y1,y2)
    
    %plot line for debugging
    %plot(x,y,'k-');
    %xlim(bounds)
    %ylim(bounds)
    %grid on;
    
    xx = linspace(x(1),x(2),30000);
    yy = linspace(y(1),y(2),30000);
    
    % get rid of the edge if it is completely outside of the boundary
    %a_logical = inpolygon(xx,yy,bounds,bounds); VERY SLOW
    a_logical = xx > bounds(1) & xx < bounds(2) & yy > bounds(1) & yy < bounds(2);
    %a_sum = sum(a_logical1); SLOWER THAN NNZ
    a_nnz = nnz(a_logical);
    
    % if a_nnz == 0, then the entire line segment is contained in the
    % exterior of the boundary. remove the line segment
    if a_nnz == 0 
        V1(1,i) = 0;
        V2(1,i) = 0;
        V1(2,i) = 0;
        V2(2,i) = 0;
        continue;
    end
    
    % if endpoints of line segment are outside boundary, truncate.
    if a_logical(1) == 0 &&  a_logical(end) == 0
       

         indx1 = find(a_logical, 1, 'first');
         indx2 = find(a_logical, 1,'last');

%         OLD METHOD (LESS ACCURATE)                 
%         V1(1,i) = xx(indx1);
%         V1(2,i) = xx(indx2);
%         V2(1,i)= yy(indx1);
%         V2(2,i) = yy(indx2);
%         
%         B(b_count+1,1) = determineType([xx(indx1),yy(indx1)],bounds);
%         B(b_count+1,2:3) = [xx(indx1),yy(indx1)];
%         
%         B(b_count+2,1) = determineType([xx(indx2),yy(indx2)],bounds);
%         B(b_count+2,2:3) = [xx(indx2),yy(indx2)];
%         b_count = b_count + 2;
%         
        
        % get exact boundary coordinates for first intersection
        %(from left to right)
        x_endpts = [x(1) xx(indx1)];
        y_endpts = [y(1) yy(indx1)];
        [xi,yi] = polyxpoly(x_endpts,y_endpts,xbox,ybox);
        V1(1,i) = xi;
        V2(1,i)= yi;
        B(b_count+1,1) = determineType([xi,yi],bounds);
        B(b_count+1,2:3) = [xi,yi];
        
        % get exact boundary coordinates for second intersection
        % (from left to right)
        x_endpts = [xx(indx2),x(2)];
        y_endpts = [yy(indx2),y(2)];
        [xi,yi] = polyxpoly(x_endpts,y_endpts,xbox,ybox);
        V1(2,i) = xi;
        V2(2,i) = yi;
        B(b_count+2,1) = determineType([xi,yi],bounds);
        B(b_count+2,2:3) = [xi,yi];
        b_count = b_count + 2;
        
        continue;
    end
    
    
    % find the point of intersection between the boundary and the line
    % calculated above.
    
    % Check right, top, bottom and replace entries in 2 (good)
    if (x(1) <= bounds(2) && x(2) >= bounds(2)) 
        [z1,z2] = right(x,y,bounds);
        if isempty([z1,z2])
            [z1,z2] = top(x,y,bounds);
            if isempty([z1,z2])
                [z1,z2] = bottom(x,y,bounds);
            end
        end
        
        if ~isempty([z1,z2])
            V1(2,i) = z1;
            V2(2,i) = z2;
            B(b_count+1,1) = determineType([z1,z2],bounds);
            B(b_count+1,2:3) = [z1,z2];
            b_count = b_count + 1;
            continue;
        end
    end
    
    % Check right, top, bottom and replace entries in 1 (good)
    if x(2) <= bounds(2) && x(1) >= bounds(2) 
        [z1,z2] = right(x,y,bounds);
        if isempty([z1,z2])
            [z1,z2] = top(x,y,bounds);
            if isempty([z1,z2])
                [z1,z2] = bottom(x,y,bounds);
            end
        end
        if ~isempty([z1,z2])
            V1(1,i) = z1;
            V2(1,i) = z2;
            B(b_count+1,1) = determineType([z1,z2],bounds);
            B(b_count+1,2:3) = [z1,z2];
            b_count = b_count + 1;
            continue;
        end
    end
    
    % Check left, top, bottom and replace entries in 1 (good)
    if x(1) <= bounds(1) && x(2) >= bounds(1) 
        [z1,z2] = left(x,y,bounds);
        if isempty([z1,z2])
            [z1,z2] = top(x,y,bounds);
            if isempty([z1,z2])
                [z1,z2] = bottom(x,y,bounds);
            end
        end
        if ~isempty([z1,z2])
            V1(1,i) = z1;
            V2(1,i) = z2;
            B(b_count+1,1) = determineType([z1,z2],bounds);
            B(b_count+1,2:3) = [z1,z2];
            b_count = b_count + 1;
            continue;
        end
    end
    
    % Check left, top, bottom and replace entries in 2 (good)
    if x(1) >= bounds(1) && x(2) <= bounds(1) 
        [z1,z2] = left(x,y,bounds);
        if isempty([z1,z2])
            [z1,z2] = top(x,y,bounds);
            if isempty([z1,z2])
                [z1,z2] = bottom(x,y,bounds);
            end
        end
        if ~isempty([z1,z2])
            V1(2,i) = z1;
            V2(2,i) = z2;
            B(b_count+1,1) = determineType([z1,z2],bounds);
            B(b_count+1,2:3) = [z1,z2];
            b_count = b_count + 1;
            continue;
        end
    end
    
    % Check top, left, right and replace entries in 1 (good)
    if y(1) >= bounds(2) && y(2) <= bounds(2) 
        [z1,z2] = top(x,y,bounds);
        if isempty([z1,z2])
            [z1,z2] = left(x,y,bounds);
            if isempty([z1,z2])
                [z1,z2] = right(x,y,bounds);
            end
        end
        if ~isempty([z1,z2])
            V1(1,i) = z1;
            V2(1,i) = z2;
            B(b_count+1,1) = determineType([z1,z2],bounds);
            B(b_count+1,2:3) = [z1,z2];
            b_count = b_count + 1;
            continue;
        end
    end
    
    
    % Check top, left, right and replace entries in 2 (good)
    if y(2) >= bounds(2) && y(1) <= bounds(2) 
        [z1,z2] = top(x,y,bounds);
        if isempty([z1,z2])
            [z1,z2] = left(x,y,bounds);
            if isempty([z1,z2])
                [z1,z2] = right(x,y,bounds);
            end
        end
        if ~isempty([z1,z2])
            V1(2,i) = z1;
            V2(2,i) = z2;
            B(b_count+1,1) = determineType([z1,z2],bounds);
            B(b_count+1,2:3) = [z1,z2];
            b_count = b_count + 1;
            continue;
        end
    end
    
    % Check bottom, left, right and replace entries in 1 (good)
    if y(1) <= bounds(1) && y(2) >= bounds(1)
        [z1,z2] = bottom(x,y,bounds);
        if isempty([z1,z2])
            [z1,z2] = left(x,y,bounds);
            if isempty([z1,z2])
                [z1,z2] = right(x,y,bounds);
            end
        end
        
        if ~isempty([z1,z2])
            V1(1,i) = z1;
            V2(1,i) = z2;
            B(b_count+1,1) = determineType([z1,z2],bounds);
            B(b_count+1,2:3) = [z1,z2];
            b_count = b_count + 1;
            continue;
        end
    end
    
    % Check bottom, left, right and replace entries in 2 (good)
    if y(1) >= bounds(1) && y(2) <= bounds(1)
        [z1,z2] = bottom(x,y,bounds);
        if isempty([z1,z2])
            [z1,z2] = left(x,y,bounds);
            if isempty([z1,z2])
                [z1,z2] = right(x,y,bounds);
            end
        end
        
        if ~isempty([z1,z2])
            V1(2,i) = z1;
            V2(2,i) = z2;
            B(b_count+1,1) = determineType([z1,z2],bounds);
            B(b_count+1,2:3) = [z1,z2];
            b_count = b_count + 1;
            continue;
        end
    end
    
    
    
end


end



%%%%% AUXILLARY FUNCTIONS FOR CLEANER CODE ABOVE

function [z1,z2] = top(x,y,bounds)
    [z1,z2] = polyxpoly(x,y,[bounds(1),bounds(2)],[bounds(2),bounds(2)]);
end


function [z1,z2] = bottom(x,y,bounds)
    [z1,z2] = polyxpoly(x,y,[bounds(1),bounds(2)],[bounds(1),bounds(1)]);
end


function [z1,z2] = left(x,y,bounds)
    [z1,z2] = polyxpoly(x,y,[bounds(1),bounds(1)],[bounds(1),bounds(2)]);
end


function [z1,z2] = right(x,y,bounds)
    [z1,z2] = polyxpoly(x,y,[bounds(2),bounds(2)],[bounds(1),bounds(2)]);
end

function T = determineType(x,bounds)
    % T = 1: vertex is on left boundary
    % T = 2: vertex is on right boundary
    % T = 3: vertex is on bottom boundary
    % T = 4: vertex is on top boundary
    %%%%%%%%%%%%%%%%%%%%%%%%
    % show which coordinate it is related to
%     if norm(x(1) - bounds(1)) < 1e-2
%         T = 1;
%     end
%     if norm(x(1) - bounds(2)) < 1e-2
%         T = 1;
%     end
%     if norm(x(2) - bounds(1)) < 1e-2
%         T = 2;
%     end
%     if norm(x(2) - bounds(2)) < 1e-2
%         T = 2;
%    end
    
    d1 = norm(x(1) - bounds(1));
    d2 = norm(x(1) - bounds(2));
    d3 = norm(x(2) - bounds(1));
    d4 =  norm(x(2) - bounds(2));
    
    min_d = min([d1,d2,d3,d4]);
    
    if min_d == d1 || min_d == d2
        T=1;
    else
        T=2;
    end
    
end













