function [data,xaxis,label]=get_spc2()
% output=get_spc()
% Read Thermo-Galactic SPC format file into MATLAB 
%updated by Paul Richardson: now outputs a data matrix, xaxis, and labels
%                            based on file names.

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
cumul(length(output)+1)=zeros;
cumul(1)=1;
for i=1:length(output)
    r(i)=size(output(i).data,2);
    cumul(i+1)=sum(r(1:i))+1;
end

for i=1:length(output)
    e=cumul(i+1)-1;
    data(:,cumul(i):e)=output(i).data;
end

data=data';
xaxis=output.xaxis;
xaxis=xaxis';
label={};

for i=1:length(output)
    n=cellstr(output(i).name);
    e=cumul(i+1)-1;
    for j=cumul(i):e
        label(j)=n;
    end
end

label=label';   
    
    
end
