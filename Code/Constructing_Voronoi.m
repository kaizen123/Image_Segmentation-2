addpath('...')

% CONSTRUCTING A VORONOI DIAGRAM
X = [7.84427890743913,9.33759848385332;2.91570317906931,1.87460806421687;6.03533438750887,2.66178838907639;9.64422667215901,7.97830260211597;4.32484993970361,4.87603775433924;6.94752194617940,7.68958264058869;7.58099275289454,3.96006745217875;4.32642326147101,2.72938794123691;6.55498039803537,0.372346340703280;1.09755050723052,6.73294914108653];

clf

plot(X(:,1),X(:,2),'k*')
hold on;
TRI = delaunay(X(:,1),X(:,2));
TRI = [TRI, TRI(:,1)];
for i = 1:size(TRI,1)
    plot(X(TRI(i,:),1),X(TRI(i,:),2),'r-')
    hold on;
end
F1 = getframe(gcf);
F1 = F1.cdata;

pause(0.2)


voronoi(X(:,1),X(:,2),'b')
F2 = getframe(gcf);
F2 = F2.cdata;


pause(0.2)

clf
[V1,V2,~,~,~,~,~,~,~,~] = myVoronoi(X,[0,10]);
voronoi(X(:,1),X(:,2),'b');
hold on;
plot(X(:,1),X(:,2),'k*')
hold on;
F3 = getframe(gcf);
F3 = F3.cdata;

% imwrite(F1,'DT.png')
% imwrite(F2,'DTV.png')
% imwrite(F3,'Voronoi.png')
% 
% for i = 1:size(V1,2)
%     midpoint = 0.5*[V1(1,i) + V1(2,i), V2(1,i) + V2(2,i)];
%     length = norm([V1(1,i) - V1(2,i), V2(1,i) - V2(2,i)]);
%     text(midpoint(1),midpoint(2), num2str(length),'FontSize',10);
%     hold on;
% end
% F4 = getframe(gcf);
% F4 = F4.cdata;
% %imwrite(F4,'lengths.png')
p1 = X(1,:);
p2 = X(6,:);
p3 = X(4,:);
V = [V1(1,4),V2(1,4)];

plot(p1(1),p1(2),'bx') 
plot(p2(1),p2(2),'bx') 
plot(p3(1),p3(2),'bx')
plot(V(1),V(2),'bx')

text(p1(1)+0.1,p1(2)+0.1,'p_1') 
text(p2(1)+0.1,p2(2)-0.1,'p_2') 
text(p3(1)-0.1,p3(2)-0.1,'p_3')
text(V(1)-0.05,V(2)+0.1,'V')



plot([p1(1) p3(1)],[p1(2) p3(2)],'b:')
plot([p1(1) p2(1)],[p1(2) p2(2)],'b:')
plot([p2(1) p3(1)],[p2(2) p3(2)],'b:')
F5 = getframe(gcf);
F5 = F5.cdata;
imwrite(F5,'gradient1.png')

xlim([6,10])
ylim([7,10])
F6 = getframe(gcf);
F6 = F6.cdata;
imwrite(F6,'gradient2.png')

p1p2 = [p2(1) - p1(1),p2(2) - p1(2)];
p1p2 = p1p2/norm(p1p2);
p1dx = [p1(1)+0.3*p1p2(1),p1(2) + 0.3*p1p2(2)];
plot(p1dx(1),p1dx(2),'r*') 
text(p1dx(1) - 0.3,p1dx(2)+0.1,'p_1 + z_1')


X2 = X;
X2(1,:) = X(1,:) + [0.3*p1p2(1) 0.3*p1p2(2)];
voronoi(X2(:,1),X2(:,2),'r');

plot([p1dx(1), p3(1)],[p1dx(2), p3(2)],'r:')
Vdv = [8.294 7.854];
plot(Vdv(1), Vdv(2), 'r*')
text(Vdv(1) + 0.05, Vdv(2) - 0.05, 'V + w_1')

plot([V(1),V(1) + 0.15*p1p2(1)],[V(2),V(2) + 0.15*p1p2(2)],'b:')
F7 = getframe(gcf);
F7 = F7.cdata;
imwrite(F7,'gradient3.png')






