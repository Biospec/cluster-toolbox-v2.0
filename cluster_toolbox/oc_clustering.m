function [] = oc_clustering(X,names,type,param,arguments)
%[] = oc_clustering(X,names,type,param,arguments)
% Performs hierarchical clustering using the OC program by
% Dr. Geoffrey J. Barton.
% This is a frontend to the program in matlab.
%
% X         = matrix with coordinates. Each row is an
%             object which is being clustered.
% names     = a string matrix with names for each object. If names = []
%             numbers from 1..N will automatically be provided
% type      = 'Mink' (for Minkowski), 'Cheb' (for Chebyshev). Types of distance 
%             measures. 'Sim'(for Similarity measures). Use the param value below to 
%             signify which similarity measure you want.
% param     = 1,2,3... Parameters to the distance measures
%             The exponent factor in Minkowski distance. 
%             Note that {Mink,1} = Manhattan distance
%             {Mink,2} = ordinary Euclidean distance
%             {Cheb,*} = Minkowski for param = infinity
%             * = means there are no parameters for Chebychev
%             If type == 'Sim' above, param has the following options:
%             1 = correlation           (CLUSTAN method no. 3)
%             2 = dispersion            (CLUSTAN method no. 32)
%             3 = cosine                (CLUSTAN method no. 27)
%             4 = dot product           (CLUSTAN method no. 26)
%             5 = similarity ratio      (CLUSTAN method no. 28)
% arguments = string arguments to the oc clustering program.
%             include any combination of keywords WITH space in the
%             arguments string. OC has 5 arguments:
%             arg1: sim (similarity matrix), dis (distance matrix),
%             arg2: single, complete, means
%             arg3: id (numer output)
%             arg4: ps <filename, without .ps extension>
%             arg5: cut <threshold>  allows OC only to output clusters above threshold
% additional options are log,timeclus and amps (produce .tree and .tord files)
% must be followed by a filename.
%
% Example usage of this front-end:
%
% oc_clustering(X,[],'Mink',2,'dis single ps mydendro amps tt');
% 
% will produce mydendro.ps with the dendrogram and tt.tree and tt.tord files. This
% is single linkage with Minkowski distance hierarchical clustering.
%
% Copyright (C), B.K. Alsberg, 25/2-1998
%

% Here we compute the distance file (later this can
% be extended to include any type of similarity matrix.

if strcmp(type,'Sim')
   D = loc_similarity_measure(X,param);
elseif strcmp(type,'Mink')
   D = gendst(X,type,param);
end;

%keyboard
%pause

[N,M]=size(X);

%the standard file name is here set:
filename = 'octemp.dat';
d = loc_unfold(D);
if isempty(names)
   names = vec2str((1:N)');
end;

if iscell(names)
   names = cell2mat(names);
end
loc_write_to_disk(d,names,filename)

% Checking for which computer type in order to
% invoke the right OC version.

ccc = computer;

switch ccc
case 'PCWIN' 
    ocprog = 'ocnt';
case 'PCWIN64'
    ocprog = 'ocnt';
case  'MAC2'
   error('Not available!');
case 'SUN4'
   ocprog = 'oc';
case 'SOL2'
   ocprog = 'oc';
case 'HP700'
   ocprog = 'oc';
case 'SGI'
   ocprog = 'oc';
case 'SGI64'
   ocprog = 'oc';
case 'IBM_RS'
   ocprog = 'oc';
case 'ALPHA'
   ocprog = 'oc';
case 'AXP_VMSG'
   ocprog = 'oc';
case 'AXP_VMSIEEE'
   ocprog = 'oc';
case 'LNX86'
   ocprog = 'oc';
case 'VAX_VMSG'
   ocprog = 'oc';
case 'VAX_VMSD'
   ocprog = 'oc';
end;

resultsfile = 'ocresults.dat';

string = [ocprog '  '  arguments '  <  ' filename '  > ' resultsfile];
runstring = ['! ' string];
eval(runstring);



%%%%%%%%%%%%%%% LOCAL FUNCTIONS %%%%%%%%%%%%%%%%%

function [] = loc_write_to_disk(x,names,filename)
% [] = loc_write_to_disk(x,names,filename)
% Writes to disk the unfolded distance matrix. x
% is therefore a column vector in this case
% names is a string matrix where the number of
% rows MUST be identical to the number of objects 

fid = fopen(filename,'w');
N = length(names);

fprintf(fid,'%d\n',N);
for i = 1:N
   fprintf(fid,'%s\n',names(i,:));
end;
K = length(x);
for i = 1:K
   fprintf(fid,'%6.3f\n',x(i));
end;

fclose(fid);


function d = loc_unfold(D)
% d = loc_unfold(D)
% unfolds the matrix according to how they expect it in oc_clustering

[n,m]=size(D);

k=1;

for i = 1:n-1
   for j = i+1:n
      d(k) = D(i,j);
      k=k+1;
   end;
end;


function D = loc_similarity_measure(X,option)
% D = loc_similarity_measure(X,option)
% Computes a similarity matrix D from X. Option
% follows the input for which measures are calculated:
% 1 = correlation           (CLUSTAN method no. 3)
% 2 = dispersion            (CLUSTAN method no. 32)
% 3 = cosine                (CLUSTAN method no. 27)
% 4 = dot product           (CLUSTAN method no. 26)
% 5 = similarity ratio      (CLUSTAN method no. 28)

[n,m]=size(X);

D = zeros(n,n);

switch option
case 1,
   D = corrcoef(X').^2;   
   D = lintrans(D,0,1,0,100);
case 2,
   mo = mean(X');
   Y = (X' - ones(m,1)*mo)';
   C = Y*Y';
   D = (1/m)*C;
case 3,
   S = diag(diag(X*X'));
   D = S^(-1/2)*X*X'*S^(-1/2);
case 4,
   D = X*X';
case 5,
   
   for i = 1:n
      for j = 1:n
         u1 = X(i,:)*X(j,:)';
         u2 = sum(X(i,:).^2);
         u3 = sum(X(j,:).^2);
         D(i,j) = u2  - u1 - u3;
      end;
   end;
otherwise
   error('You asked for a similarity measure that does not exist - bombing out');
end;


   

