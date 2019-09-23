function [rGrad,rEnergy] = repulsion(X,dist,theta)

alpha = 1;
r = dist;
%theta = 10;


% Compute distances between all points
x_coords = X(:,1);
y_coords = X(:,2);
diffx = x_coords - x_coords';
diffy = y_coords - y_coords';
diffxs = diffx.*diffx;
diffys = diffy.*diffy;
norms = (diffxs + diffys).^(1/2);

% Compute the mask for activating the repulsion function
%upperNorms = triu(norms);
upperNorms = norms;
mask1 = upperNorms > 0 & upperNorms < r;


% Compute the energy contribution
%rEnergy = norms.^(-2) - delta^(-2);

exponent = -1./ ((theta*theta)*(r - norms).^2);
rEnergy = exp(exponent);
rEnergy(isinf(rEnergy)) = 0;
rEnergy = alpha * rEnergy.*mask1;
rEnergy = sum(sum(rEnergy));


% Get the gradient of the repulsion function
%norms4 = norms.^(-4);
%RX1 = -2*diffx.*norms4;
RX1 = alpha * exp(exponent) * 2 * theta^(-2) * (r - norms).^(-3) .* norms.^(-1).*diffx;
%RX2 = -2*diffy.*norms4;
RX2 = alpha * exp(exponent) * 2 * theta^(-2) * (r - norms).^(-3) .* norms.^(-1).*diffy;
RX1(isnan(RX1)) = 0;
RX2(isnan(RX2)) = 0;

RX1 = RX1.*mask1;
RX2 = RX2.*mask1;

RX1 = sum(RX1,2);
RX2 = sum(RX2,2);
rGrad = [RX1,RX2];

for i = 1:size(rGrad,1)
    if norm(rGrad(i,:)) == 0
        continue;
    end
    rGrad(i,:) = rGrad(i,:)/norm(rGrad(i,:));
end
end































% function rGrad = repulsion(X)
%     push = zeros(size(X));
%     mask = push;
% 
%     delta = 0.01;
% 
%     x_coords = X(:,1);
%     y_coords = X(:,2);
%     plot(x_coords,y_coords,'rx')
% 
%     % COMPUTE DISTANCES BETWEEN ALL POINTS
% 
%     x_diff = x_coords - x_coords';
%     y_diff = y_coords - y_coords';
% 
%     x_diff_sq = x_diff.*x_diff;
%     y_diff_sq = y_diff.*y_diff;
% 
%     norm_sq = x_diff_sq + y_diff_sq;
%     norms = norm_sq.^(1/2);
% 
%     % since the diagonal will be zeros
%     %(points are zero distance away from themselves)
%     norms = triu(norms);
%     [I,J] = find(norms < delta & norms > 0);
%     for i = 1:size(I)
%         center1 = I(i);
%         center2 = J(i);
%     
%         direction = X(center1,:) - X(center2,:);
%         
%         push_dist = norms(center1,center2)^(-2) - delta^(-2)
%         push_dist = push_dist*0.01;
%         %push_dist = delta - norms(center1,center2);
% 
%         push(center1,:) = -push_dist * direction;
%         push(center2,:) = push_dist * direction;
%         mask(center1,:) = 1;
%         mask(center2,:) = 1;
%     end
% end
