function SaveImages(Image,ImageNumber)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

h = Image;
set(h,'NextPlot','replacechildren');
set(h,'Renderer','zbuffer');
% Attach the right suffix to the image
if ImageNumber < 10
   numstr = ['0000',num2str(ImageNumber)];
elseif ImageNumber < 100
   numstr = ['000',num2str(ImageNumber)];
elseif ImageNumber < 1000
   numstr = ['00',num2str(ImageNumber)];
elseif ImageNumber < 10000
   numstr = ['0',num2str(ImageNumber)];
else
   numstr = num2str(ImageNumber);
end       
% Save the figure as a .png file
saveas(h, ['frm', numstr, '.png']);
          
% Clear the figure
clf; 
end

%ADAPTED FROM: https://groups.csail.mit.edu/vision/sliwiki/index.php?title=Creating_Matlab_Videos

