% lower and upper bounds for the boundary box of the voronoi diagram
bounds = [0,10];

% Number of centers that will be used in the voronoi diagra
numCenters = 4;

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
minArea = 0.5;

% Step size for the gradient of volumetric term
vStep = 0.009;

% Step size for discretization of line integral
delta = 0.5;

% Image to segment
pic = ones(11, 11);
pic(1:5, 1:5) = 10; %1
pic(6:11, 1:5) = 20;%4
pic(1:5, 6:11) = 30;%2
pic(6:11, 6:11) = 40;%3



% 
% imagesc(pic)
% 
% pause()

%F1 = getframe(gcf);
%F1 = F1.cdata;
%imwrite(F1,'Sample_image.png')


% Run image segmentation
voronoiSegment(bounds,numCenters,maxIters,pStep,rStep,dist,theta,minArea,vStep,delta,pic)
