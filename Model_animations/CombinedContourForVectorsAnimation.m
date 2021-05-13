function M = CombinedContourForVectorsAnimation(Columns,Lines,NumberOfFrames,Var1,Var2,Var3,Var4,step)
% This function produces a frame combining a surface colormap with contour lines representing
% two different variables
% Function arguments:
% Columns, Lines - Horizontal and vertical dimensions of the contour plot
% NumberOfFrames - Number of frames to include in the movie
% Var1 and Var2 - 3D Arrays, where the first dimension is time, the second is
% the number of lines and the third the number of columns (the same for
% both variables)
% Var1 - Vector variable (the U component)
% Var2 - Vector variable (the V component)
% Var3 - Contour plot variable used to locate Var1 in space relative to Var2 (e.g. bathymetry)
%aviobj = avifile('MyMovie.avi','Compression','Indeo3','quality',100); 
for k = 1:NumberOfFrames
   for i = 1:Lines
      for j = 1:Columns
         U(i,j) = Var1(k,Lines+1-i,j);
         V(i,j) = Var2(k,Lines+1-i,j);
         Z(i,j) = Var3(k,Lines+1-i,j);
         W(i,j) = Var4(k,Lines+1-i,j);
      end
   end
   CombinedContourForVectors(Columns,Lines,U,V,Z,W,step)
   Frame(k) = getframe(gca);
   hold off
   %%aviobj = addframe(aviobj,Frame(k));
end
M = Frame;
%movie2avi(M,'MyMovie.avi','Compression','Cinepak','quality',100);
movie2avi(M,'Movie.avi','Compression','None');
%%aviobj = close(aviobj);