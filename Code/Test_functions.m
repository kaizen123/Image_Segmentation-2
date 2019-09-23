bounds = [0,2];
% x_centers = [9.6309 4.8890;5.4681 6.2406;5.2114 6.7914;2.3159 3.9552]/5;
%x_centers = 2*rand(9,2);
pic = [10,12,18;8,33,22;14,38,41];
%const = [20, 30, 40, 10];

[V1,V2,IND,intV,vDepX,xGenV,B,BX,D1,D2] = myVoronoi(x_centers,bounds,'p');

xpqTuples = getXPQTuples(xGenV,vDepX);

DVDX = gradientNonBoundary(x_centers,xpqTuples);

[DBDX, xGenB] = gradientBoudaryToCenter(x_centers, B, BX);
%The output here: energy is a scalar, gradient is the same shape of
%x_centers. 
% Just add them to the energy and gradent of the precious algorithm. You
% can also get some weighted sum, of course. 
[energy, gradient] = energyDiffGradient(xGenV, xGenB, intV, B, x_centers, DVDX, DBDX, pic, 0.001)

% [energy, grad1, grad2] = energyEdgeLenGradient(V1, V2);
 
