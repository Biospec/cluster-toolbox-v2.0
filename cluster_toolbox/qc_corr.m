function output = qc_corr(Xin, QC_idxin, batch_idx)
%output = qc_corr(X, QC_idx, batch_idx)
%X : data matrix with QCs
%QC_idx : QC indices
%batch_idx : optional, batch number, all ones (ie. one batch only) if ignored

[m,n]=size(Xin);
% convert 0s to NaNs
% Xin(Xin==0)=NaN;

if nargin<3
    batch_idx=ones(m,1);
end

X_QC_all=Xin(QC_idxin,:);
QCs_all=zeros(m,1);
QCs_all(QC_idxin)=1;

% remove variables with more than 1/3 missing values in the QCs
nan_count=zeros(size(X_QC_all));
nan_count(isnan(X_QC_all))=1;
nan_count=sum(nan_count);
null_var=find(nan_count >= size(X_QC_all,1)/3);

disp('calculating');
Driftout=[]; Ratioout=[];
batch_ids=unique(batch_idx);
no_batches=length(unique(batch_idx));
batch_corr=zeros(no_batches,n);
for i=1:no_batches
    X=Xin(batch_idx==batch_ids(i),:);
    QCs_inner=QCs_all(batch_idx==batch_ids(i));
    [m,n]=size(X);
    X_corr=X;
    driftlines=zeros(m,n);
    ratio=zeros(n,1);
    disp(['Batch no. = ' num2str(i)])
    h=waitbar(0,'Please wait...');
    for ii=1:n
        QC_idx=find(QCs_inner==1);
        if any(ismember(null_var,ii))
            continue;
        end
        Q_vals=X(QC_idx,ii);
        % outlier detection
        outlier_Qs=outlier_detection(Q_vals);
        % remove outliers
        if ~isempty(outlier_Qs)
            QC_idx(outlier_Qs)=[];
        end
        vec_test=nan(m,1);
        vec_test(QC_idx)=X(QC_idx,ii);
        % Median absolute deviation, robust estimation of variance
        medQ=nanmedian(vec_test);
        madQ=nanmedian(abs(vec_test-medQ));
        vec_sample=X(:,ii);
        vec_sample(QC_idx)=[];
        medS=nanmedian(vec_sample);
        madS=nanmedian(abs(vec_sample-medS));
        % ratio between the variance in QC and that in samples, ideally it
        % should be <<1
        ratio(ii)=madQ/madS;
        [vec_smoothed, ~, exitflag]=smoothn(vec_test,'robust','MaxIter',500);
        if exitflag==1 & min(vec_smoothed>=0)
        % the spline algorithm was successfully converged and the
        % smoothline did not dip below 0
            vec_corr=X(:,ii)+(max(vec_smoothed)-vec_smoothed);
            X_corr(:,ii)=vec_corr;
            driftlines(:,ii)=vec_smoothed;
        else
        % the spline algorithm had failed to converge, something must
        % be wrong in the data, skip the correction
            null_var=[null_var ii];
            continue;
        end
        waitbar(ii/n,h);
    end
    Xout{i}=X_corr;
    batch_corr(i,:)=nanmedian(X_corr(QC_idx,:));
    if any(isnan(batch_corr(i,:)))
        nan_idx=find(isnan(batch_corr(i,:)));
        disp(['Warning! Excessive amount of missing values in the QC samples in batch '...
            num2str(i) ' , variable id. = ' num2str(nan_idx)]);
        disp('No batch correct was applied to these variables');
        disp('Removing these variables is recommended');
        batch_corr(i,isnan(batch_corr(i,:)))=0;
        null_var=[null_var nan_idx];
    end
    Driftout=[Driftout; driftlines];
    Ratioout=[Ratioout; ratio];
    close(h);
end
if no_batches>1
    for i=1:n
        if any(ismember(null_var,ii))
            continue;
        end
        vec_corr_batch=max(batch_corr(:,i))-batch_corr(:,i);
        % apply batch correction factors
        for ii=1:no_batches
            Xout{ii}(:,i)=Xout{ii}(:,i)+vec_corr_batch(ii);
        end
    end
    X_corr_final=[];
    for i=1:no_batches
        X_corr_final=[X_corr_final; Xout{i}];
    end
else
    X_corr_final=Xout{1};
end

output.X_corr=X_corr_final;
output.driftlines=Driftout;
output.ratios=Ratioout;
output.null_vars=null_var;
disp('Done!');
end

function [z,s,exitflag,Wtot] = smoothn(varargin)
% %SMOOTHN Robust spline smoothing for 1-D to N-D data.
%   SMOOTHN provides a fast, automatized and robust discretized smoothing
%   spline for data of any dimension.
%
%   Z = SMOOTHN(Y) automatically smoothes the uniformly-sampled array Y. Y
%   can be any N-D noisy array (time series, images, 3D data,...). Non
%   finite data (NaN or Inf) are treated as missing values.
%
%   Z = SMOOTHN(Y,S) smoothes the array Y using the smoothing parameter S.
%   S must be a real positive scalar. The larger S is, the smoother the
%   output will be. If the smoothing parameter S is omitted (see previous
%   option) or empty (i.e. S = []), it is automatically determined using
%   the generalized cross-validation (GCV) method.
%
%   Z = SMOOTHN(Y,W) or Z = SMOOTHN(Y,W,S) specifies a weighting array W of
%   real positive values, that must have the same size as Y. Note that a
%   nil weight corresponds to a missing value.
%
%   Robust smoothing
%   ----------------
%   Z = SMOOTHN(...,'robust') carries out a robust smoothing that minimizes
%   the influence of outlying data.
%
%   [Z,S] = SMOOTHN(...) also returns the calculated value for S so that
%   you can fine-tune the smoothing subsequently if needed.
%
%   An iteration process is used in the presence of weighted and/or missing
%   values. Z = SMOOTHN(...,OPTION_NAME,OPTION_VALUE) smoothes with the
%   termination parameters specified by OPTION_NAME and OPTION_VALUE. They
%   can contain the following criteria:
%       -----------------
%       TolZ:       Termination tolerance on Z (default = 1e-3)
%                   TolZ must be in ]0,1[
%       MaxIter:    Maximum number of iterations allowed (default = 100)
%       Initial:    Initial value for the iterative process (default =
%                   original data)
%       -----------------
%   Syntax: [Z,...] = SMOOTHN(...,'MaxIter',500,'TolZ',1e-4,'Initial',Z0);
%
%   [Z,S,EXITFLAG] = SMOOTHN(...) returns a boolean value EXITFLAG that
%   describes the exit condition of SMOOTHN:
%       1       SMOOTHN converged.
%       0       Maximum number of iterations was reached.
%
%   Class Support
%   -------------
%   Input array can be numeric or logical. The returned array is of class
%   double.
%
%   Reference
%   --------- 
%   Garcia D, Robust smoothing of gridded data in one and higher dimensions
%   with missing values. Computational Statistics & Data Analysis, 2010. 
%   <a
%   href="matlab:web('http://www.biomecardio.com/pageshtm/publi/csda10.pdf')">PDF download</a>
%

% Check input arguments
error(nargchk(1,12,nargin));
%% Test & prepare the variables
%---
k = 0;
while k<nargin && ~ischar(varargin{k+1}), k = k+1; end
%---
% y = array to be smoothed
y = double(varargin{1});
sizy = size(y);
noe = prod(sizy); % number of elements
if noe<2, z = y; return, end
%---
% Smoothness parameter and weights
W = ones(sizy);
s = [];
if k==2
    if isempty(varargin{2}) || isscalar(varargin{2}) % smoothn(y,s)
        s = varargin{2}; % smoothness parameter
    else % smoothn(y,W)
        W = varargin{2}; % weight array
    end
elseif k==3 % smoothn(y,W,s)
        W = varargin{2}; % weight array
        s = varargin{3}; % smoothness parameter
end
if ~isequal(size(W),sizy)
        error('MATLAB:smoothn:SizeMismatch',...
            'Arrays for data and weights must have same size.')
elseif ~isempty(s) && (~isscalar(s) || s<0)
    error('MATLAB:smoothn:IncorrectSmoothingParameter',...
        'The smoothing parameter must be a scalar >=0')
end
%---
% "Maximal number of iterations" criterion
I = find(strcmpi(varargin,'MaxIter'),1);
if isempty(I)
    MaxIter = 100; % default value for MaxIter
else
    try
        MaxIter = varargin{I+1};
    catch ME
        error('MATLAB:smoothn:IncorrectMaxIter',...
            'MaxIter must be an integer >=1')
    end
    if ~isnumeric(MaxIter) || ~isscalar(MaxIter) ||...
            MaxIter<1 || MaxIter~=round(MaxIter)
        error('MATLAB:smoothn:IncorrectMaxIter',...
            'MaxIter must be an integer >=1')        
    end    
end
%---
% "Tolerance on smoothed output" criterion
I = find(strcmpi(varargin,'TolZ'),1);
if isempty(I)
    TolZ = 1e-3; % default value for TolZ
else
    try
        TolZ = varargin{I+1};
    catch ME
        error('MATLAB:smoothn:IncorrectTolZ',...
            'TolZ must be in ]0,1[')
    end
    if ~isnumeric(TolZ) || ~isscalar(TolZ) || TolZ<=0 || TolZ>=1 
        error('MATLAB:smoothn:IncorrectTolZ',...
            'TolZ must be in ]0,1[')
    end    
end
%---
% "Initial Guess" criterion
I = find(strcmpi(varargin,'Initial'),1);
if isempty(I)
    isinitial = false; % default value for TolZ
else
    isinitial = true;
    try
        z0 = varargin{I+1};
    catch ME
        error('MATLAB:smoothn:IncorrectInitialGuess',...
            'Z0 must be a valid initial guess for Z')
    end
    if ~isnumeric(z0) || ~isequal(size(z0),sizy) 
        error('MATLAB:smoothn:IncorrectTolZ',...
            'Z0 must be a valid initial guess for Z')
    end    
end
%---
% "Weighting function" criterion
I = find(strcmpi(varargin,'Weights'),1);
if isempty(I)
    weightstr = 'bisquare'; % default weighting function
else
    try
        weightstr = lower(varargin{I+1});
    catch ME
        error('MATLAB:smoothn:IncorrectWeights',...
            'A valid weighting function must be chosen')
    end
    if ~ischar(weightstr)
        error('MATLAB:smoothn:IncorrectWeights',...
            'A valid weighting function must be chosen')
    end    
end
%---
% Weights. Zero weights are assigned to not finite values (Inf or NaN),
% (Inf/NaN values = missing data).
IsFinite = isfinite(y);
nof = nnz(IsFinite); % number of finite elements
W = W.*IsFinite;
if any(W<0)
    error('MATLAB:smoothn:NegativeWeights',...
        'Weights must all be >=0')
else 
    W = W/max(W(:));
end
%---
% Weighted or missing data?
isweighted = any(W(:)<1);
%---
% Robust smoothing?
isrobust = any(strcmpi(varargin,'robust'));
%---
% Automatic smoothing?
isauto = isempty(s);
%---
% % DCTN and IDCTN are required
% test4DCTNandIDCTN

%% Creation of the Lambda tensor
%---
% Lambda contains the eingenvalues of the difference matrix used in this
% penalized least squares process.
d = ndims(y);
Lambda = zeros(sizy);
for i = 1:d
    siz0 = ones(1,d);
    siz0(i) = sizy(i);
    Lambda = bsxfun(@plus,Lambda,...
        cos(pi*(reshape(1:sizy(i),siz0)-1)/sizy(i)));
end
Lambda = -2*(d-Lambda);
if ~isauto, Gamma = 1./(1+s*Lambda.^2); end

%% Upper and lower bound for the smoothness parameter
% The average leverage (h) is by definition in [0 1]. Weak smoothing occurs
% if h is close to 1, while over-smoothing appears when h is near 0. Upper
% and lower bounds for h are given to avoid under- or over-smoothing. See
% equation relating h to the smoothness parameter (Equation #12 in the
% referenced CSDA paper).
N = sum(sizy~=1); % tensor rank of the y-array
hMin = 1e-6; hMax = 0.99;
sMinBnd = (((1+sqrt(1+8*hMax.^(2/N)))/4./hMax.^(2/N)).^2-1)/16;
sMaxBnd = (((1+sqrt(1+8*hMin.^(2/N)))/4./hMin.^(2/N)).^2-1)/16;

%% Initialize before iterating
%---
Wtot = W;
%--- Initial conditions for z
if isweighted
    %--- With weighted/missing data
    % An initial guess is provided to ensure faster convergence. For that
    % purpose, a nearest neighbor interpolation followed by a coarse
    % smoothing are performed.
    %---
    if isinitial % an initial guess (z0) has been provided
        z = z0;
    else
        z = InitialGuess(y,IsFinite);
    end
else
    z = zeros(sizy);
end
%---
z0 = z;
y(~IsFinite) = 0; % arbitrary values for missing y-data
%---
tol = 1;
RobustIterativeProcess = true;
RobustStep = 1;
nit = 0;
%--- Error on p. Smoothness parameter s = 10^p
errp = 0.1;
opt = optimset('TolX',errp);
%--- Relaxation factor RF: to speedup convergence
RF = 1 + 0.75*isweighted;

%% Main iterative process
%---
while RobustIterativeProcess
    %--- "amount" of weights (see the function GCVscore)
    aow = sum(Wtot(:))/noe; % 0 < aow <= 1
    %---
    while tol>TolZ && nit<MaxIter
        nit = nit+1;
        DCTy = dctn(Wtot.*(y-z)+z);
        if isauto && ~rem(log2(nit),1)
            %---
            % The generalized cross-validation (GCV) method is used.
            % We seek the smoothing parameter s that minimizes the GCV
            % score i.e. s = Argmin(GCVscore).
            % Because this process is time-consuming, it is performed from
            % time to time (when nit is a power of 2)
            %---
            fminbnd(@gcv,log10(sMinBnd),log10(sMaxBnd),opt);
        end
        z = RF*idctn(Gamma.*DCTy) + (1-RF)*z;
        
        % if no weighted/missing data => tol=0 (no iteration)
        tol = isweighted*norm(z0(:)-z(:))/norm(z(:));
       
        z0 = z; % re-initialization
    end
    exitflag = nit<MaxIter;

    if isrobust %-- Robust Smoothing: iteratively re-weighted process
        %--- average leverage
        h = sqrt(1+16*s); h = sqrt(1+h)/sqrt(2)/h; h = h^N;
        %--- take robust weights into account
        Wtot = W.*RobustWeights(y-z,IsFinite,h,weightstr);
        %--- re-initialize for another iterative weighted process
        isweighted = true; tol = 1; nit = 0; 
        %---
        RobustStep = RobustStep+1;
        RobustIterativeProcess = RobustStep<4; % 3 robust steps are enough.
    else
        RobustIterativeProcess = false; % stop the whole process
    end
end

%% Warning messages
%---
if isauto
    if abs(log10(s)-log10(sMinBnd))<errp
        warning('MATLAB:smoothn:SLowerBound',...
            ['s = ' num2str(s,'%.3e') ': the lower bound for s ',...
            'has been reached. Put s as an input variable if required.'])
    elseif abs(log10(s)-log10(sMaxBnd))<errp
        warning('MATLAB:smoothn:SUpperBound',...
            ['s = ' num2str(s,'%.3e') ': the upper bound for s ',...
            'has been reached. Put s as an input variable if required.'])
    end
end
if nargout<3 && ~exitflag
    warning('MATLAB:smoothn:MaxIter',...
        ['Maximum number of iterations (' int2str(MaxIter) ') has ',...
        'been exceeded. Increase MaxIter option or decrease TolZ value.'])
end


    %% GCV score
    %---
    function GCVscore = gcv(p)
        % Search the smoothing parameter s that minimizes the GCV score
        %---
        s = 10^p;
        Gamma = 1./(1+s*Lambda.^2);
        %--- RSS = Residual sum-of-squares
        if aow>0.9 % aow = 1 means that all of the data are equally weighted
            % very much faster: does not require any inverse DCT
            RSS = norm(DCTy(:).*(Gamma(:)-1))^2;
        else
            % take account of the weights to calculate RSS:
            yhat = idctn(Gamma.*DCTy);
            RSS = norm(sqrt(Wtot(IsFinite)).*(y(IsFinite)-yhat(IsFinite)))^2;
        end
        %---
        TrH = sum(Gamma(:));
        GCVscore = RSS/nof/(1-TrH/noe)^2;
    end



    %% Robust weights
    function W = RobustWeights(r,I,h,wstr)
        % weights for robust smoothing.
        MAD = median(abs(r(I)-median(r(I)))); % median absolute deviation
        u = abs(r/(1.4826*MAD)/sqrt(1-h)); % studentized residuals
        if strcmp(wstr,'cauchy')
            c = 2.385; W = 1./(1+(u/c).^2); % Cauchy weights
        elseif strcmp(wstr,'talworth')
            c = 2.795; W = u<c; % Talworth weights
        else
            c = 4.685; W = (1-(u/c).^2).^2.*((u/c)<1); % bisquare weights
        end
        W(isnan(W)) = 0;
    end



    %% Initial Guess with weighted/missing data
    function z = InitialGuess(y,I)
        %-- nearest neighbor interpolation (in case of missing values)
        if any(~I(:))
            if license('test','image_toolbox')
                [z,L] = bwdist(I);
                z = y;
                z(~I) = y(L(~I));
            else
            % If BWDIST does not exist, NaN values are all replaced with the
            % same scalar. The initial guess is not optimal and a warning
            % message thus appears.
                z = y;
                z(~I) = mean(y(I));
                warning('MATLAB:smoothn:InitialGuess',...
                    ['BWDIST (Image Processing Toolbox) does not exist. ',...
                    'The initial guess may not be optimal; additional',...
                    ' iterations can thus be required to ensure complete',...
                    ' convergence. Increase ''MaxIter'' criterion if necessary.'])    
            end
        else
            z = y;
        end
        %-- coarse fast smoothing using one-tenth of the DCT coefficients
        siz = size(z);
        z = dctn(z);
        for k = 1:ndims(z)
            z(ceil(siz(k)/10)+1:end,:) = 0;
            z = reshape(z,circshift(siz,[0 1-k]));
            z = shiftdim(z,1);
        end
        z = idctn(z);
    end
end


function [y,w] = dctn(y,DIM,w)

%DCTN N-D discrete cosine transform.
%   Y = DCTN(X) returns the discrete cosine transform of X. The array Y is
%   the same size as X and contains the discrete cosine transform
%   coefficients. This transform can be inverted using IDCTN.
%
%   DCTN(X,DIM) applies the DCTN operation across the dimension DIM.
%
%   Class Support
%   -------------
%   Input array can be numeric or logical. The returned array is of class
%   double.
%
%   Reference
%   ---------
%   Narasimha M. et al, On the computation of the discrete cosine
%   transform, IEEE Trans Comm, 26, 6, 1978, pp 934-936.
%
%   Example
%   -------
%       RGB = imread('autumn.tif');
%       I = rgb2gray(RGB);
%       J = dctn(I);
%       imshow(log(abs(J)),[]), colormap(jet), colorbar
%
%   The commands below set values less than magnitude 10 in the DCT matrix
%   to zero, then reconstruct the image using the inverse DCT.
%
%       J(abs(J)<10) = 0;
%       K = idctn(J);
%       figure, imshow(I)
%       figure, imshow(K,[0 255])
%
%   See also IDCTN, DCT, DCT2.
%
%   -- Damien Garcia -- 2008/06, revised 2009/11
%   website: <a
%   href="matlab:web('http://www.biomecardio.com')">www.BiomeCardio.com</a>

% ----------
%   [Y,W] = DCTN(X,DIM,W) uses and returns the weights which are used by the
%   program. If DCTN is required for several large arrays of same size, the
%   weights can be reused to make the algorithm faster. A typical syntax is
%   the following:
%      w = [];
%      for k = 1:10
%          [y{k},w] = dctn(x{k},[],w);
%      end
%   The weights (w) are calculated during the first call of DCTN then
%   reused in the next calls.
% ----------

narginchk(1,3)

y = double(y);
sizy = size(y);

% Test DIM argument
if ~exist('DIM','var'), DIM = []; end
assert(~isempty(DIM) || ~isscalar(DIM),...
    'DIM must be a scalar or an empty array')
assert(isempty(DIM) || DIM==round(DIM) && DIM>0,...
    'Dimension argument must be a positive integer scalar within indexing range.')

% If DIM is empty, a DCT is performed across each dimension

if isempty(DIM), y = squeeze(y); end % Working across singleton dimensions is useless
dimy = ndims(y);

% Some modifications are required if Y is a vector
if isvector(y)
    dimy = 1;
    if size(y,1)==1
        if DIM==1, w = []; return
        elseif DIM==2, DIM=1;
        end
        y = y.';
    elseif DIM==2, w = []; return
    end
end

% Weighting vectors
if ~exist('w','var') || isempty(w)
    w = cell(1,dimy);
    for dim = 1:dimy
        if ~isempty(DIM) && dim~=DIM, continue, end
        n = (dimy==1)*numel(y) + (dimy>1)*sizy(dim);
        w{dim} = exp(1i*(0:n-1)'*pi/2/n);
    end
end

% --- DCT algorithm ---
if ~isreal(y)
    y = complex(dctn(real(y),DIM,w),dctn(imag(y),DIM,w));
else
    for dim = 1:dimy
        if ~isempty(DIM) && dim~=DIM
            y = shiftdim(y,1);
            continue
        end
        siz = size(y);
        n = siz(1);
        y = y([1:2:n 2*floor(n/2):-2:2],:);
        y = reshape(y,n,[]);
        y = y*sqrt(2*n);
        y = ifft(y,[],1);
        y = bsxfun(@times,y,w{dim});
        y = real(y);
        y(1,:) = y(1,:)/sqrt(2);
        y = reshape(y,siz);
        y = shiftdim(y,1);
    end
end
        
y = reshape(y,sizy);
end

function [y,w] = idctn(y,DIM,w)

%IDCTN N-D inverse discrete cosine transform.
%   X = IDCTN(Y) inverts the N-D DCT transform, returning the original
%   array if Y was obtained using Y = DCTN(X).
%
%   IDCTN(X,DIM) applies the IDCTN operation across the dimension DIM.
%
%   Class Support
%   -------------
%   Input array can be numeric or logical. The returned array is of class
%   double.
%
%   Reference
%   ---------
%   Narasimha M. et al, On the computation of the discrete cosine
%   transform, IEEE Trans Comm, 26, 6, 1978, pp 934-936.
%
%   Example
%   -------
%       RGB = imread('autumn.tif');
%       I = rgb2gray(RGB);
%       J = dctn(I);
%       imshow(log(abs(J)),[]), colormap(jet), colorbar
%
%   The commands below set values less than magnitude 10 in the DCT matrix
%   to zero, then reconstruct the image using the inverse DCT.
%
%       J(abs(J)<10) = 0;
%       K = idctn(J);
%       figure, imshow(I)
%       figure, imshow(K,[0 255])
%
%   See also DCTN, IDSTN, IDCT, IDCT2.
%
%   -- Damien Garcia -- 2009/04, revised 2009/11
%   website: <a
%   href="matlab:web('http://www.biomecardio.com')">www.BiomeCardio.com</a>

% ----------
%   [Y,W] = IDCTN(X,DIM,W) uses and returns the weights which are used by
%   the program. If IDCTN is required for several large arrays of same
%   size, the weights can be reused to make the algorithm faster. A typical
%   syntax is the following:
%      w = [];
%      for k = 1:10
%          [y{k},w] = idctn(x{k},[],w);
%      end
%   The weights (w) are calculated during the first call of IDCTN then
%   reused in the next calls.
% ----------

narginchk(1,3)

y = double(y);
sizy = size(y);

% Test DIM argument
if ~exist('DIM','var'), DIM = []; end
assert(~isempty(DIM) || ~isscalar(DIM),...
    'DIM must be a scalar or an empty array')
assert(isempty(DIM) || DIM==round(DIM) && DIM>0,...
    'Dimension argument must be a positive integer scalar within indexing range.')

% If DIM is empty, a DCT is performed across each dimension

if isempty(DIM), y = squeeze(y); end % Working across singleton dimensions is useless
dimy = ndims(y);

% Some modifications are required if Y is a vector
if isvector(y)
    dimy = 1;
    if size(y,1)==1
        if DIM==1, w = []; return
        elseif DIM==2, DIM=1;
        end
        y = y.';
    elseif DIM==2, w = []; return
    end
end

% Weighing vectors
if ~exist('w','var') || isempty(w)
    w = cell(1,dimy);
    for dim = 1:dimy
        if ~isempty(DIM) && dim~=DIM, continue, end
        n = (dimy==1)*numel(y) + (dimy>1)*sizy(dim);
        w{dim} = exp(1i*(0:n-1)'*pi/2/n);
    end
end

% --- IDCT algorithm ---
if ~isreal(y)
    y = complex(idctn(real(y),DIM,w),idctn(imag(y),DIM,w));
else
    for dim = 1:dimy
        if ~isempty(DIM) && dim~=DIM
            y = shiftdim(y,1);
            continue
        end        
        siz = size(y);
        n = siz(1);
        y = reshape(y,n,[]);
        y = bsxfun(@times,y,w{dim});
        y(1,:) = y(1,:)/sqrt(2);
        y = ifft(y,[],1);
        y = real(y*sqrt(2*n));
        I = (1:n)*0.5+0.5;
        I(2:2:end) = n-I(1:2:end-1)+1;
        y = y(I,:);
        y = reshape(y,siz);
        y = shiftdim(y,1);            
    end
end
        
y = reshape(y,sizy);
end

function outlier_idx=outlier_detection(values)
% Outlier detection based on box-whisker plot / IQR 
% Only extreme cases were considered, 
% i.e. > 3*IQR + Q(75%) or < Q(25%) - 3*IQR


n=length(values);
sort_vals=sort(values);
low_IQR=sort_vals(round(n*.25));
high_IQR=sort_vals(round(n*.75));
IQR=high_IQR-low_IQR;
upper_bound=high_IQR+IQR*3;
lower_bound=low_IQR-IQR*3;
outlier_idx=find(values>upper_bound | values<lower_bound);

end
