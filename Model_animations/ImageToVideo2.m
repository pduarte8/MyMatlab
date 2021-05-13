function ImageToVideo2()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
cd /home/pduarte/pytools_git
ImFolder=uigetdir
index = 0;
for i = 1:120
    for j = 0:23
        Day = i;
        if Day < 10
            SDay = strcat('000',int2str(Day));
        elseif Day < 100
            SDay = strcat('00',int2str(Day));
        elseif Day < 1000
            SDay = strcat('0',int2str(Day));
        else
            SDay = int2str(Day); 
        end
        Hour = j;
        if Hour < 10
            SHour = strcat('0',int2str(Hour));
        else
            SHour = int2str(Hour); 
        end    
        index = index + 1
        pngFiles(index,:) = strcat('NO3_ocean_his_',SDay,'_s10_',SHour,'.png');
    end
end
VideoFile=strcat(ImFolder,'\Video');
writeObj = VideoWriter(VideoFile);
fps= 15; 
writeObj.FrameRate = fps;
open(writeObj);
for t= 1:index
   Frame=imread(strcat(ImFolder,'/',pngFiles(t,:)));
   writeVideo(writeObj,im2frame(Frame));
end
close(writeObj);
end


%Explanation of the above code: 

%  ImFolder=uigetdir;  --> This line when executed, will open a browse window in which you can select your desired folder, where the frames or sequence of images is. The path will be stored in the variable ImFolder.

% pngFile = dir(strcat(ImFolder,'\*.png')); --> here the variable pngFile will be a structure that contain all the names of PNG files that are in the directory you have specified. If your image files are in jpeg format, so give the above command as dir(strcat(ImFolder,'\*.jpg')). Here '\*.png' will be acting as a wildcard for all the PNG files available in that folder.

% S = [pngFile(:).datenum]; 
%[S,S] = sort(S);
%--> This 2 lines of code will sort the files according to the datenum property of the pngFile structure.

% pngFiles = pngFile(S); --> creating another variable pngFiles same as pngFile , but the former contains the files list in sorted form.

% VideoFile=strcat(ImFolder,'\Video'); --> Creating a Video File name as Video

% writeObj = VideoWriter(VideoFile); --> This will create a Video Object named writeObj.

% fps= 15; --> specifying the frame rate as 15 (or whatever is required or whatever you want). Your video file length will be decided by the value of fps and the number of the image files.

% writeObj.FrameRate = fps; -->  we have set the FrameRate property of writeObj  video object as specified in the variable fps.

%  open(writeObj); --> Open to start writing the frames in the writeObj video object.

% for t= 1:length(pngFiles)
%     Frame=imread(strcat(ImFolder,'\',pngFiles(t).name));
%     writeVideo(writeObj,im2frame(Frame));
%end
%--> These 4 lines will write the frames in the newly created video file, that is initially empty.
% close(writeObj); --> close the writeObj to complete the video write process.


%tAKWN FROM: http://www.divilabs.com/2013/12/create-video-file-from-sequence-of.html
