function ThreeDVelocityPlot(Time,U,V,W)
%U(Time, Layers, Lines, Columns)
A = size(U);
for k = 1:A(2)
    for i = 1:A(3)
      for j = 1:A(4)
         x(k,i,j)=j;
         y(k,i,j)=i;
         z(k,i,j)=k;
      end
   end
end
for k = 1:A(2)
   for i = 1:10
       for j = 1:10
         u(k,i,j)=U(Time,k,i,j);
         v(k,A(3)+1-i,j)=V(Time,k,i,j);
         w(A(2)+1-k,i,j)=W(Time,k,i,j);
       end
   end
end
quiver3(x(:,:,2:A(4)),y(:,:,2:A(4)),z(:,:,2:A(4)),u(:,:,2:A(4)),v(:,:,2:A(4)),w(:,:,2:A(4)),1);