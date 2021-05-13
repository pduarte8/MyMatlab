function CombinedContourForVectorsV2(DistanceY,DistanceX,Var1,Var2,Var3,Var4,step, cmin,cmax)
% This function produces a surface colormap with contour lines representing
% two different variables
% Function arguments:
% Columns, Lines - Horizontal and vertical dimensions of the contour plot
% Var1, Var2, Var3 and Var4 - 2D Arrays, where the first dimension is the number of lines and the second is
% the number of columns (the same for
% both variables)
% Var1 and Var2 - U and V velocity components
% Var3 and Var4 - Colormap variable (e.g. the concentration of a dissolved
% substance and depth)

%size(DistanceY)
%size(DistanceX)

x = DistanceY(1:size(DistanceY));
y = DistanceX(1:size(DistanceX));
%size(x)
%size(y)
[m,Columns] = size(x);
[n,Lines] = size(y);
[X,Y] = meshgrid(x,y);
%contour(X,Y,Var3);
hold on;
%size(X)
%size(Y)
%size(Var4)
pcolor(X,Y,Var4);
colorbar;
hold on
caxis([cmin cmax]);
shading('interp');
axis equal;
XMIN = 0;
XMAX = x(1,Columns);
YMIN = 0;
YMAX = y(1,Lines);
axis([XMIN XMAX YMIN YMAX]);
xlabel('m');
ylabel('m');
hold on;
quiver(x(1:Columns),y(1:Lines),Var1(1:Lines,1:Columns),-Var2(1:Lines,1:Columns),3,'k');