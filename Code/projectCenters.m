function X1 = projectCenters(X,Gradient,bounds)

X1 = X + Gradient;

box_x = [bounds(1), bounds(1), bounds(2), bounds(2), bounds(1)]';
box_y = [bounds(1), bounds(2), bounds(2), bounds(1), bounds(1)]';

[in,~] = inpolygon(X1(:,1),X1(:,2),box_x,box_y);

centers_outside = find(in == 0);

%%%%% Find which boundary the point should be projected onto %%%%%
for i = 1 : length(centers_outside)
    
    indx = centers_outside(i);
    
    % If the new x_coordinate is to the left of the boundary box
    if X1(indx,1) < bounds(1)
        X1(indx,1) = bounds(1);
        
        % if the new y_coordinate is also outside the boundary, move to corner
        if X1(indx,2) > bounds(2)
            X1(indx,2) = bounds(2);
        end
        if X1(indx,2) < bounds(1)
            X1(indx,2) = bounds(1);
        end
    end
    
    % If the new x_coordinate is to the right of the boundary box
    if X1(indx,1) > bounds(2)
        X1(indx,1) = bounds(2);
        
        % if the new y_coordinate is also outside the boundary, move to corner
        if X1(indx,2) > bounds(2)
            X1(indx,2) = bounds(2);
        end
        if X1(indx,2) < bounds(1)
            X1(indx,2) = bounds(1);
        end
        %    end
    end
    
    %%%%%%%%%%
    
    
    % If the new y_coordinate is below the boundary box
    if X1(indx,2) < bounds(1)
        X1(indx,2) = bounds(1);
        if X1(indx,1) < bounds(1)
            X1(indx,1) = bounds(1);
        end
        if X1(indx,1) > bounds(2)
            X1(indx,1) = bounds(2);
        end
    end
    
    % If the new y_coordinate is to above the boundary box
    if X1(indx,2) > bounds(2)
        X1(indx,2) = bounds(2);
        
        % if the new x_coordinate is also outside the boundary, move to corner
        if X1(indx,1) < bounds(1)
            X1(indx,1) = bounds(1);
        end
        if X1(indx,1) > bounds(2)
            X1(indx,1) = bounds(2);
        end
    end
    
end   
end
