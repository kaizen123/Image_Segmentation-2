addpath('..');
close all;


tic
% Set up the data
bounds = [0 99];
pic = ones(100, 100);
pic(1:50, 1:50) = 10; %1
pic(51:100, 1:50) = 20;%4
pic(1:50, 51:100) = 30;%2
pic(51:100, 51:100) = 40;%3
x_centers_init = [ 30.1 30.2; 30.3 50.4; 80.5 70.6;60.7 10.8];


%rng(7);
%rng(20)
%x_centers_init = 10*rand(4,2);

%x_centers_init = x_centers;
[V1,V2,IND,intV,vDepX,xGenV,B,BX,D1,D2] = myVoronoi(x_centers_init,bounds,'p');
xpqTuples = getXPQTuples(xGenV,vDepX);
DVDX = gradientNonBoundary(x_centers_init,xpqTuples);
[DBDX, xGenB] = gradientBoudaryToCenter(x_centers_init, B, BX);
%[energy, grad] = energyEdgeLenGradient(V1, V2);
%gradient = chainRule(xGenV, xGenB, DVDX, DBDX, grad, D1, D2);
[energy, gradient] = energyDiffGradient(xGenV, xGenB, intV, B, x_centers_init, DVDX, DBDX, vDepX, BX, pic, 0.1)

%for i = 1:size(gradient)
%   gradient(i,:) = gradient(i,:)/norm(gradient(i,:)) ;
%end

x_centers_init = [ 30.1 30.7; 30.3 50.4; 80.5 70.6;60.7 10.8];
[V1,V2,IND,intV,vDepX,xGenV,B,BX,D1,D2] = myVoronoi(x_centers_init,bounds,'p');
xpqTuples = getXPQTuples(xGenV,vDepX);
DVDX = gradientNonBoundary(x_centers_init,xpqTuples);
[DBDX, xGenB] = gradientBoudaryToCenter(x_centers_init, B, BX);
%[energy, grad] = energyEdgeLenGradient(V1, V2);
%gradient = chainRule(xGenV, xGenB, DVDX, DBDX, grad, D1, D2);
[energy_, gradient] = energyDiffGradient(xGenV, xGenB, intV, B, x_centers_init, DVDX, DBDX, vDepX, BX, pic, 0.1);
a = (energy_-energy)/0.5


% Get Current voronoi diagram and store it
F(1) = getframe(gcf);
clf;

max_iterations = 2000;
rate = 0.1;

% Get Current Energy and Store it
Energies = zeros(1,max_iterations);
Energies(1) = energy;
Gradient_norm = zeros(1,max_iterations);
Gradient_norm(1) = norm(gradient);



x_centers = x_centers_init;

for iter = 2:max_iterations
    clf;
    x_centers = x_centers - rate*gradient;
    %Compute new data
    [V1,V2,IND,intV,vDepX,xGenV,B,BX,D1,D2] = myVoronoi(x_centers,bounds,'p');
    hold on;
    plot(x_centers(:,1),x_centers(:,2),'*');
    title(['Voronoi GD w/ Learn Rate ',num2str(rate),' ( Iteration: ',num2str(iter),' )'])
    xpqTuples = getXPQTuples(xGenV,vDepX);
    DVDX = gradientNonBoundary(x_centers,xpqTuples);
    [DBDX, xGenB] = gradientBoudaryToCenter(x_centers, B, BX);
    %[energy, grad] = energyEdgeLenGradient(V1, V2);
    %gradient = chainRule(xGenV, xGenB, DVDX, DBDX, grad, D1, D2);
    [energy, gradient] = energyDiffGradient(xGenV, xGenB, intV, B, x_centers, DVDX, DBDX, vDepX, BX, pic, 0.01)
    Energies(iter) = energy;
    Gradient_norm(iter) = norm(gradient);
    
    % Get the new voronoi diagram and store it
    F(iter) = getframe(gcf);
    
    % Stopping Condition
    if norm(Gradient_norm(iter - 1) - Gradient_norm(iter)) < 1e-5
        break;
    end
    
    
energy    
    
end


figure(2)
plot(1:iter,Energies(1:iter))
title('Energy per iteration');
xlabel('Iteration Number')
ylabel('Energy')
toc