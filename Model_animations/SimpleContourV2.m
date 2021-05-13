function Figure = SimpleContourV2(Columns,Lines,Var1)
x = [1:Columns];
y = [1:Lines];
[X,Y] = meshgrid(x,y);

Figure = figure;
pcolor(X,Y,Var1);
shading('interp');
caxis([0 5]);

%colorbar;

set(colorbar, 'ylim', [0 5])


