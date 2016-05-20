function output=get_maldi()

[files, path]=uigetfile('*.txt','Select File','Multiselect','on');
[~,n]=size(files);
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
    fid = fopen(filename);
    output(p).name=file
    temp=textscan(fid,'%f %f','headerlines',1);
    ms=temp{1};
    [sort_ms, sort_idx]=sort(ms);
    sort_intensity=temp{2}(sort_idx);
    output(p).ms = sort_ms;
    output(p).intensity = sort_intensity;
    fclose(fid);
    clear ms sort_ms sort_idx;
end
