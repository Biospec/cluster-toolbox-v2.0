function output=get_cdf_tic(filename)

[files, path]=uigetfile('*.cdf','Select File','Multiselect','on');
[m,n]=size(files);
if ~iscell(files)
    n=1;
end
for p=1:n
    if n==1
        file=files;
    else
        file=char(files(p))
    end
    filename = [path file] ;    
    ncid = netcdf.open(filename, 'nowrite') ;
    output(p).name=file
    rt_id = netcdf.inqVarID (ncid, 'scan_acquisition_time' );
    output(p).rt = netcdf.getVar (ncid, rt_id); 
    tic_id = netcdf.inqVarID (ncid, 'total_intensity' );
    output(p).tic = netcdf.getVar (ncid, tic_id);
%     mz_id = netcdf.inqVarID (ncid, 'mass_values' );
%     tmp.mass = netcdf.getVar (ncid, mz_id);
%     ms_id = netcdf.inqVarID (ncid, 'intensity_values' );
%     tmp.int = netcdf.getVar (ncid, ms_id); 
%     scan_id = netcdf.inqVarID (ncid, 'scan_index' );
%     tmp.scan = netcdf.getVar (ncid, scan_id);
%     netcdf.close(ncid);
%     tmp.time=length(tmp.scan);  
%     mass_len=length(tmp.mass);
%     tmp.scan(tmp.time+1)=mass_len+1;
%     maxMass = round(max(tmp.mass)) ;
%     minMass = round(min(tmp.mass)) ;
%     MZ.masses = minMass : maxMass ;
%     scan=1;
%     % Define the starting zero matrix
%     output(p).matrix=zeros(tmp.time,length(MZ.masses));
%     for j = 1:mass_len
%         if j<=tmp.scan(scan+1)
%             output(p).matrix(scan,(round(tmp.mass(j))-minMass+1))=...
%                 output(p).matrix(scan,(round(tmp.mass(j))-minMass+1))+tmp.int(j);
%         else
%             scan=scan+1;
%             output(p).matrix(scan,(round(tmp.mass(j))-minMass+1))=tmp.int(j);
%         end
%     end
%     output(p).matrix=sparse(output(p).matrix);
%     output(p).MaS = MZ.masses;
    clear tmp;
    netcdf.close(ncid)
end