function voronoiGD(bounds,numCenters,maxIters,pStep,rStep,dist,theta,minArea)

X = bounds(2)*rand(numCenters,2);

% compute the voronoi diagram and all the data structures
[V1,V2,IND,intV,vDepX,xGenV,B,BX,D1,D2] = myVoronoi(X,bounds,'p');
hold on;
grid on;

% xpqTuples : which three centers generate each interior vertex
xpqTuples = getXPQTuples(xGenV,vDepX);

% DVDX : change in interior vertices wrt to changes in centers
DVDX = gradientNonBoundary(X,xpqTuples);

% DBDX : changes in boundary vertices wrt to changes in centers
% xGenB : which two centers generate which boundary vertices
[DBDX, xGenB] = gradientBoudaryToCenter(X, B, BX);

% pEnergy : Perimiter Energy
% grad : differences in x ,differences in y, of the edges
[pEnergy, grad] = energyEdgeLenGradient(V1, V2);

% gradients of vertices wrt to centers
pGradient = chainRule(xGenV, xGenB, DVDX, DBDX, grad, D1, D2);
pGradient = normalize(pGradient);

% rGrad : repulsion gradient
% rEnergy: energy contributed by the repulsion function
[rGrad,rEnergy] = repulsion(X,dist,theta);

% Get Current voronoi diagram and store it
G= getframe(gcf);
F(:,:,:,1) = G.cdata;

% Get Current Energy and Store it
Energies = zeros(1,maxIters);
Energies(1) = pEnergy + rEnergy;

x_centers = X;


for iter = 2:maxIters
    
    % Use corresponding rates for gradients calculated and pass into
    % the next function. Need to add in Yi's gradient from above.
    update = - pStep*pGradient - rStep*rGrad;
    
    % Project the centers onto the boundary if the gradients move them
    % outside the boundary
    x_centers = projectCenters(x_centers,update,bounds);
    
    % Remove the centers with small area
    x_centers = removeCenters(x_centers,xGenV,intV,xGenB,B,minArea,bounds);
    
    % Calculate the new voronoi diagram and all other data structures 
    [V1,V2,IND,intV,vDepX,xGenV,B,BX,D1,D2] = myVoronoi(x_centers,bounds,'p');
    hold on;
    grid on;
    xpqTuples = getXPQTuples(xGenV,vDepX);
   
    DVDX = gradientNonBoundary(x_centers,xpqTuples);
    
    [DBDX, xGenB] = gradientBoudaryToCenter(x_centers, B, BX);
    
    [pEnergy, grad] = energyEdgeLenGradient(V1, V2);
    
    pGradient = chainRule(xGenV, xGenB, DVDX, DBDX, grad, D1, D2);
    pGradient = normalize(pGradient);
    
    [rGrad,rEnergy] = repulsion(x_centers,dist,theta);

    % Need to add in Yi's volumetric energy calculation
    Energies(iter) = pEnergy + rEnergy;

    % Get the new voronoi diagram and store it
    G= getframe(gcf);
    F(:,:,:,iter) = G.cdata;
    
    % Stopping Condition
%     if Energies(iter) - Energies(iter - 1) < 1e-8
%         break;
%     end
    clf;
    
end

% Plot the energy as a function of time after gradient descent has stopped
figure(2)
plot(1:iter,Energies(1:iter))
title('Energy per iteration');
xlabel('Iteration Number')
ylabel('Energy')

%  vidObj = VideoWriter('GD_RP','MPEG-4');
%     open(vidObj);
%  
%     % Create an animation.
%     
%     axis tight manual
%     set(gca,'nextplot','replacechildren');
%  
%     for k = 1:size(F,4)
%        imshow([F(:,:,:,k)]);
%        % Write each frame to the file.
%        currFrame = getframe(gcf);
%        writeVideo(vidObj,currFrame);
%     end
%   
%     % Close the file.
%     close(vidObj);


end





function Output = normalize(Input)

Output = zeros(size(Input));
for i = 1:size(Input,1)
    if norm(Input(i,:)) == 0
       continue; 
    end
    
    if sum(isnan(Input(i,:))) > 0
       Output(i,:) = [0 0]; 
    end
    Output(i,:) = Input(i,:) / norm(Input(i,:));
    
end

end
