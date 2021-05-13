function M = CombinedContourForVectorsAnimation2DHorizFrom3DDataV2(Columns,Lines,Layer,Var1,Var2,Var3,step)
%function M = CombinedContourForVectorsAnimation2DHorizFrom3DData(Columns,Lines,Layer,NumberOfFrames,Var1,Var2,Var3,Var4,step)
% This function produces a frame combining a surface colormap with contour lines representing
% two different variables
% Function arguments:
% Columns, Lines - Horizontal and vertical dimensions of the contour plot
% NumberOfFrames - Number of frames to include in the movie
% Var1, Var2 and Var3 - 4D Arrays, where the first dimension is time, the
% second is the number of vertical layers, the third
% the number of lines and the fourth the number of columns (the same for
% all variables)
% Var1 - Vector variable (the U component)
% Var2 - Vector variable (the V component)
% Var3 - Depth of each grid cell -  It is used to compute total depth by summing its values for each line and column across all layers 
% Var4 - Surface plot water quality variable
A = size(Var3);
NumberOfLayers = A(2);

for k = 1:NumberOfFrames
   Z = zeros(Lines,Columns);  
   for i = 1:Lines
      for j = 1:Columns
         U(i,j) = Var1(k,Layer,Lines+1-i,j);
         V(i,j) = Var2(k,Layer,Lines+1-i,j);
          for m = 1:A(2)
            Z(i,j) = Z(i,j) + Var3(k,m,Lines+1-i,j);
          end
         %W(i,j) = Var4(k,Layer,Lines+1-i,j);
      end
   end
   H=figure;
   CombinedContourForVectors(Columns,Lines,U,V,Z,Z,step);
   Frame(k) = getframe(H);
   delete(H);
end
M = Frame;
movie2avi(M,'Movie1.avi','Compression','None');