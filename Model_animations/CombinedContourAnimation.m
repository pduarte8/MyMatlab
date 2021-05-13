function M = CombinedContourAnimation(Columns,Lines,NumberOfFrames,Var1,Var2)
% This function produces a frame combining a surface colormap with contour lines representing
% two different variables
% Function arguments:
% Columns, Lines - Horizontal and vertical dimensions of the contour plot
% NumberOfFrames - Number of frames to include in the movie
% Var1 and Var2 - 3D Arrays, where the first dimension is time, the second is
% the number of lines and the third the number of columns (the same for
% both variables)
% Var1 - Colormap variable (e.g. the concentration of a dissolved
% substance)
% Var2 - Contour plot variable used to locate Var1 in space relative to Var2 (e.g. bathymetry)

for k = 1:NumberOfFrames
   for i = 1:Lines
      for j = 1:Columns
         Z1(i,j) = Var1(k,Lines+1-i,j);
         Z2(i,j) = Var2(k,Lines+1-i,j);
      end
   end
   CombinedContour(Columns,Lines,Z1,Z2)
   Frame(k) = getframe;
   hold off
   %imwrite(Frame(k),'1.hdf');
end
M = Frame;

movie2avi(M,'Movie.avi','Compression','None');
