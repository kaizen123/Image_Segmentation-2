% Image to segment
pic = imread('im3.jpg');
pic = double(pic);
pic = imresize(pic,0.05);
sz = size(pic);

% lower and upper bounds for the boundary box of the voronoi diagram
bounds = [0,sz(1)-1];

% Number of centers that will be used in the voronoi diagra
numCenters = 4;

% Maximum number of iterations of gradient descent
maxIters = 1000;

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
minArea = 0.5;

% Step size for the gradient of volumetric term
vStep = 0.009;

% Step size for discretization of line integral
delta = 0.001;

% Run image segmentation
voronoiSegment(bounds,numCenters,maxIters,pStep,rStep,dist,theta,minArea,vStep,delta,pic)
