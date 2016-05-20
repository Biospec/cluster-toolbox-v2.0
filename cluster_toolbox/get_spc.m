function output=get_spc()
% output=get_spc()
% Read Thermo-Galactic SPC format file into MATLAB 

[files, path]=uigetfile('*.spc','Select File','Multiselect','on');
[m,n]=size(files);
if ~iscell(files)
    n=1;
end
for p=1:n
    if n==1
        file=files;
    else
        file=char(files(p));
    end
    filename = [path file] ;    
    SData = tgspcread(filename);
    if ~isfield(SData,'spectra')
        output(p).name = file;
        output(p).xaxis = SData.X; 
        output(p).data = double(SData.Y);
    else 
        output(p).name = file;
        output(p).xaxis = SData.X; 
        no_spectra=length(SData.Y);
        for k=1:no_spectra
            output(p).data{k}=double(SData.spectra(k).data);
        end
    end
end
