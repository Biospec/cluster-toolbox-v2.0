function output=get_cdf()
% output=get_cdf()
% Read netcdf files into Matlab
% The data will be interpolated to a common retention time window across
% all the samples.

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
    disp([num2str(n) ' files in total, loading no. ' [num2str(p)] '...']);
    filename = [path file] ;    
    ncid = netcdf.open(filename, 'nowrite') ;
    output(p).name=file;
    rt_id = netcdf.inqVarID (ncid, 'scan_acquisition_time' );
    output(p).rt = netcdf.getVar (ncid, rt_id); 
    min_rt(p)=min(output(p).rt);
    max_rt(p)=max(output(p).rt);
    tic_id = netcdf.inqVarID (ncid, 'total_intensity' );
    output(p).tic = netcdf.getVar (ncid, tic_id);
    mz_id = netcdf.inqVarID (ncid, 'mass_values' );
    tmp.mass = netcdf.getVar (ncid, mz_id);
    ms_id = netcdf.inqVarID (ncid, 'intensity_values' );
    tmp.int = netcdf.getVar (ncid, ms_id); 
    scan_id = netcdf.inqVarID (ncid, 'scan_index' );
    tmp.scan = netcdf.getVar (ncid, scan_id);
    netcdf.close(ncid);
    tmp.time=length(tmp.scan);  
    mass_len=length(tmp.mass);
    tmp.scan(tmp.time+1)=mass_len+1;
    maxMass = round(max(tmp.mass)) ;
    minMass = round(min(tmp.mass)) ;
    MZ.masses = minMass : maxMass ;
    scan=1;
    % Define the starting zero matrix
    output(p).matrix=zeros(tmp.time,length(MZ.masses));
    for j = 1:mass_len
        if j<=tmp.scan(scan+1)
            output(p).matrix(scan,(round(tmp.mass(j))-minMass+1))=...
                output(p).matrix(scan,(round(tmp.mass(j))-minMass+1))+tmp.int(j);
        else
            scan=scan+1;
            output(p).matrix(scan,(round(tmp.mass(j))-minMass+1))=tmp.int(j);
        end
    end
    output(p).matrix=sparse(output(p).matrix);
    output(p).MaS = MZ.masses;
    clear tmp;
end

disp('loading finished, start interpolation')
ref_rt=max(min_rt):mean(diff(output(1).rt)):min(max_rt);

for i=1:n
    disp(['Processing no.' num2str(i) '...']);
    tic_tmp=interp1(output(i).rt, output(i).tic, ref_rt(:));
    no_ms=size(output(i).matrix,2);
    mat=full(output(i).matrix);mat2=[];
    for ii=1:no_ms
        mat2(:,ii)=interp1(output(i).rt, mat(:,ii), ref_rt(:));
    end
    output(i).matrix=sparse(mat2); output(i).tic=tic_tmp;
    output(i).rt=ref_rt;
end
disp('All done!');