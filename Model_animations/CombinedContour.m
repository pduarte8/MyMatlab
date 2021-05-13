x = [xfunction CombinedContour(Columns,Lines,Var1,Var2)
% This function produces a surface colormap with contour lines representing
% two different variables
% Function arguments:
% Columns, Lines - Horizontal and vertical dimensions of the contour plot
% Var1 - Colormap variable (e.g. the concentration of a dissolved
% substance)
% Var2 - Contour plot variable used to locate Var1 in space relative to Var2 (e.g. bathymetry)
x = [1:Columns];
y = [1:Lines];
[X,Y] = meshgrid(x,y);
%contour(X,Y,Var1);
pcolor(X,Y,Var1);
shading('interp');
colorbar;
hold on;
contour(X,Y,Var2,'w');