function output=get_cdf_dims(filename)

[files, path]=uigetfile('*.cdf','Select File','Multiselect','on');
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
    ncid = netcdf.open(filename, 'nowrite') ;
    output(p).name=file
    scan_id = netcdf.inqVarID (ncid, 'scan_index' );
    scan_idx= netcdf.getVar (ncid, scan_id);
    mz_id = netcdf.inqVarID (ncid, 'mass_values' );
    mass = netcdf.getVar (ncid, mz_id); 
    ms_id = netcdf.inqVarID (ncid, 'intensity_values' );
    intensity = netcdf.getVar (ncid, ms_id);
    unique_ms=unique(mass);
%     matrix=zeros(length(scan_idx), length(unique_ms));
    for i=1:length(scan_idx)-1
        output(p).intensity{i}=...
            intensity(scan_idx(i)+1:scan_idx(i+1));
        output(p).mass{i}=mass(scan_idx(i)+1:scan_idx(i+1));
    end
    output(p).intensity{i+1}=intensity(scan_idx(end)+1:end);
    output(p).mass{i+1}=mass(scan_idx(end)+1:end);
end
