function [ndata, mdata] =emsc(data, order, varargin);
% EMSC   : EXTENDED MULTIPLICATIVE SCATTER CORRECTION (SEE H. MARTENS)
%          correct all spectra to a single background  
% I/O    : [ndata, mdata] = emsc(data,order, index);
%        : [ndata, mdata] = emsc(data,order, fdata);
%        : [ndata, mdata] = emsc(data,order, fdata, idata);
% data   : data to be corrected (N x M: N = number of points, M = number of spectra)
% order  : order of polynomial order for the background used in fitting
% index  : index of spectrum of data (or fdata, if provided) to which all other 
%          spectra are fitted (common background spectrum) (optional),
%          if zero, not provided, or empty, the mean of data is taken as independent spectrum
% fdata  : spectrum to which all other spectra are fitted (optional)
% idata  : other known interferences of which the variance is to be elimnated (optional)
% ndata  : corrected data
% mdata  : spectrum to which all other spectra were fitted

% CODT, ErasmusMC 2005: RW, PJC, TBS


index = 0;
mdata = [];
idata = [];

% determine varargin, index and/or independent fit spectrum
if nargin < 2
    
    error('need at least two input arguments, can handle four at the most');
    
elseif nargin == 2
    
    % fit data with mean
    mdata=mean(data,2);
    
elseif nargin >= 3
    
    if length(varargin{1}) == 1 % only index
        
        index = varargin{1};
        if index > size(data,3)
            error('index and data do not match');
        end
        if index == 0
            mdata = mean(data,2);
        else
            mdata = data(:,index);
        end
        
    else % fitspectrum
        
        mdata = varargin{1};
        if isempty(mdata)
            mdata=mean(data,2);
        else
            if size(data,1) ~= size(mdata,1)
                error('data sizes do not match');
            end
        end
    end
    if nargin == 4 % 
        
        idata = varargin{2};
        if size(data,1) ~= size(idata,1)
            error('data sizes do not match');
        end
        
    end
end

if ~isempty(idata)
    mdata = [mdata idata];
end;

ns = size(mdata,2);

% actual fitting
ndata = data;
for i=1:size(data,2)
    
        [b,fit,res] = multifit(data(:,i),mdata,order);
        ndata(:,i) = res./b(1) + mdata(:,1);

end

mdata = mdata(:,1);




%===========================
% internal functions
%===========================

function [b,fit,res]=multifit(y,model,varargin)
%MULTIFIT Multiple regression fit of the columns of matrix X to vector y
% Unrestricted multiple-least-squares-regression fit with optional polynomial
% background. Fit-coefficients can be positive and negative. A third input argument
% (optional) can be used to include a polynomial of order N in the model. With N=1, for
% example, y will be fitted by an extra offset and slope. If N is not specified, y will
% be fitted by the columns of X only.
%
% I/O:  [b,fit,residual]=multifit(y,X,order)
%   y     : dependent variable (column vector [m x 1])
%   X     : independent variables (array [m x n])
%   order : (optional) order of a background polynomial, which can be
%           included in the set of independent variables. If order is
%           not specified, no background will be included.
%
%   b        : fit-coefficients, if 'order' is specified the coefficients
%              for the polynomial background are appended to b (highest
%              order terms last)  
%   fit      : fit result (column vector [m x 1])
%   residual : residual   (column vector [m x 1])

% PJC, 15-1-2001

if nargin==2 %do not include a background in fitmodel
    X=model;
else % include polynomial of order 'order' in fitmodel      
    order=varargin{1};
    m=length(y);
    background=ones(m,1);
    if order>0
        for j=1:order
            background=[background [0:1/(m-1):1]'.^j];
        end
    end
    X=[model background];
end

%actual fit (b contains fit-coefficients)
b=pinv(X'*X)*X'*y;
warning on;
fit=X*b;
res=y-fit;

