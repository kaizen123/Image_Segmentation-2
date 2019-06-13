% lower and upper bounds for the boundary box of the voronoi diagram
bounds = [0,10];

% Number of centers that will be used in the voronoi diagra
numCenters = 100;

% Maximum number of iterations of gradient descent
maxIters = 500;

% Step size for the gradient of perimiter energy wrt centers
pStep = 0.003;

% Step size for the gradient of repulsion energy wrt centers
rStep = 0.21;

% Distance between centers for repulsion to activate
dist = 0.2;

% Parameter for repulsion that controls the steepness of the repulsion
% curve
theta = 10;

% Parameter for removing centers whose region has area smaller than minArea
minArea = 0.3;

% Run gradient descent
voronoiGD(bounds,numCenters,maxIters,pStep,rStep,dist,theta,minArea)