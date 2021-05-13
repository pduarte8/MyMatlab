function [SurfacePAR,VerticallyIntSubPAR,Chlorophyll] = LeituraDeDataSetsDeFicheiroHDF
prompt = {'Prefix for hdf files:','Number of hdf files to be read:','Initial line:','Final line','Initial column:','Final column:'};
dlg_title = 'Reading of hdf files';
num_lines= 1;
FileName = 'Scenario1'
nfiles = 24;
InitLin = 40;
FinLine = 121;
InitCol = 1;
FinCol = 81;


def  = {FileName,int2str(nfiles),int2str(InitLin),int2str(FinLine),int2str(InitCol),int2str(FinCol)};
answer  = inputdlg(prompt,dlg_title,num_lines,def);
FileName = answer(1)
nfiles = str2double(answer(2));
InitialLine = str2double(answer(3));
FinalLine = str2double(answer(4));
InitialColumn = str2double(answer(5));
FinalColumn = str2double(answer(6));

time = 1;

for n = 1:nfiles
    n
    HDFFile = FileName;
    S = int2str(n-1);
    HDFFile = strcat(FileName,'_',S,'.hdf')
    HDFFile = char(HDFFile);
    if (n == 1)
        y = hdfread(HDFFile,'Phytoplankton biomass');
        [a,b,c] = size(y);
        t = a*nfiles;
        lines = FinLine-(InitLin-1)
        columns = FinCol-(InitCol-1)
        SurfacePAR = zeros(t,lines,columns); 
        VerticallyIntSubPAR = zeros(t,lines,columns);
        Chlorophyll = zeros(t,lines,columns);
    end   
    SurfacePAR(time:time+9,:,:) = hdfread(HDFFile,'PAR surface irradiance with ice','Index',{[1 InitLin InitCol],[1 1 1],[10 lines columns]});
    VerticallyIntSubPAR(time:time+9,:,:) = hdfread(HDFFile,'Mean horizontal water PAR irradiance','Index',{[1 InitLin InitCol],[1 1 1],[10 lines columns]});
    Chlorophyll(time:time+9,:,:) = hdfread(HDFFile,'Phytoplankton biomass','Index',{[1 InitLin InitCol],[1 1 1],[10 lines columns]});
    time = time + 10;
end  %files