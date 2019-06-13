function X = removeCenters(X,xGenV,intV,xGenB,B,minArea,bounds)



W = bounds(2); H = bounds(2);
% find the center that is closest each corner
corners = [0, 0; W, 0; 0, H; W, H];
corner2x = [];
n = size(corners, 1);
for c = 1:n
    dis = sum((X-corners(c, :)).^2, 2);
    [~, x] = min(dis);
    corner2x = [corner2x, x];
end 

centers2delete = [];
for x = 1:size(X,1)
    [~, polygon] = sortV(transpose(xGenV{x}), xGenB{x}, find(corner2x==x),intV, B, corners, X(x, :));
    area = polyarea(polygon(:,1),polygon(:,2));
    %voronoi(X(:,1),X(:,2))
    %hold on
    %%grid on
    %plot(X(x,1),X(x,2),'rx')
    if area < minArea
        centers2delete = [centers2delete , x];
    end
end



% Remove the centers based on total area of its region.
X(centers2delete,:) = [];

end
function [idxSort, coordinates] = sortV(vs, bs, cs, V, B, corners, center)
    idxSort = [(1:size(vs, 2)), -(1:size(bs, 2)), zeros(size(cs))];
    coordinates = [V(vs,:);B(bs,2:3);corners(cs, :)];
    ang = atan2(coordinates(:, 2) - center(2), coordinates(:, 1) - center(1));
    [~, order] = sort(ang);
    idxSort = idxSort(order);
    coordinates = coordinates(order, :);
    coordinates = [coordinates; coordinates(1,:)];
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Start with only interior regions
% % if the X(i,:) is an interior region, xGenV{i} = [(non empty)]
% % and xGenB{i} = [(empty)]
% interior_regions_indices = [];
% for i = 1:size(X,1)
%     if isempty(xGenB{i})
%         interior_regions_indices = [interior_regions_indices, i];
%     end
% end
% 
% centers2delete = [];
% for i = 1:length(interior_regions_indices)
%     % Get the index for the current center
%     index = interior_regions_indices(i);
%     % get vector of indices into intV of vertices that depend on this
%     % center
%     dependent_vertices = xGenV{index};
%     % Center them at the current center
%     vectors = intV(dependent_vertices,:) - X(index,:);
%     % get the angles between the vertex to center and the horizontal axis
%     angles = atan2(vectors(:,2),vectors(:,1));
%     % Sort
%     [angles_sorted,sortedIndices] = sort(angles);
%     vertices = intV(dependent_vertices,:);
%     sorted_vertices = vertices(sortedIndices,:);
%     x_coords = sorted_vertices(:,1);
%     y_coords = sorted_vertices(:,2);
%     x_coords = [x_coords; x_coords(1)];
%     y_coords = [y_coords; y_coords(1)];
%     
%     region_area = polyarea(x_coords,y_coords);
%     if region_area < minArea
%         centers2delete = [centers2delete,index];
%     end
% end
% 
% % Now calculate the area of centers which are not interior regions.
% % need to identify which regions contain the corners of the
% 

