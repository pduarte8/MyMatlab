function CombinedContourForVectors(Columns,Lines,Var1,Var2,Var3,Var4,step)
% This function produces a surface colormap with contour lines representing
% two different variables
% Function arguments:
% Columns, Lines - Horizontal and vertical dimensions of the contour plot
% Var1, Var2 and Var3 - 2D Arrays, where the first dimension is the number of lines and the second is
% the number of columns (the same for
% both variables)
% Var1 - Colormap variable (e.g. the concentration of a dissolved
% substance)
% Var2 and Var3 - U and V velocity components

x = [1:Columns];
y = [1:Lines];
[X,Y] = meshgrid(x,y);
%contour(X,Y,Var3);
hold on
%pcolor(X,Y,Var4);
%shading('interp');
contour(X,Y,Var4);
caxis([-2 10]);
%colorbar;
hold on;
quiver(2:step:Columns,1:step:Lines,Var1(1:step:Lines,2:step:Columns),-Var2(1:step:Lines,2:step:Columns),4,'k');