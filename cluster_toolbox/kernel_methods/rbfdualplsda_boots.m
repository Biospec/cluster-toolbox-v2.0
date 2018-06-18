function output = rbfdualplsda_boots(data, Y, rep, no_loops)
% output = rbfdualplsr_boots(data, conc, rep, no_pcs, no_loops)
% data: pre-processed data
% Y: the class information vector
% rep: an index vector indicating which samples are from the same bio-rep.
%      all the reps of the same biological sample shall have the same
%      number. e.g. the rep_idx of an experiment having 5 biological samples 
%      with 3 analytical reps each should be [1;1;1;2;2;2;3;3;3;4;4;4;5;5;5]. 
%      Set to 1:no. of samples if not specified, i.e. assuming each bio-rep
%      has only one sample.
% no_loops: the number of iterations of bootstrapping resampling, set to
%       1000 if not specified.
% output
%   .Cknown: the know concentrations of each bootstrapping test set
%   .Cknown_av: averaged known concetrations of all the bootstrapping test sets
%   .Cpred: the predicted concentrations of each bootstrapping test sets
%   .Cpred_av: averaged predicted concetrations of all the bootstrapping test sets 
%   .Cpred_std: the std of the predicted concentrations of all the bootstrapping test sets 
%   .PLS_loadings: The PLS loadings
%   .Q2_boots: cross-validated R2 of all the bootstrapping training sets,
%              used for selecting the optimal PLS components.
%   .Q2p: R2 of all the bootstrapping test sets.
%   .no_pcs: the optimal PLS components selected for each bootstrapping
%            iteration.
%   .r2: R2 of the full data set, using the minimal no. of PLS components
%       used by bootstrapping iterations.
%   .rmsec: Root mean square error of the full data set, using the minimal no. 
%          of PLS components used by bootstrapping iterations.
%   .rmsecv_boots: Root mean square error of cross-validation of each of the
%          bootstrapping training set.
%   .rmsep: Root mean square error of the test set of each of the
%          bootstrapping test set.
%   .vip_scores: Variable Importance in Projection for each of the
%          variables.


%Last update: 26/07/2016 - permutation test added Yun Xu

    if nargin<2
        help dualplsr_validation;
        output=NaN;
        return;
    end
    if nargin==2
        rep=[1:size(data,1)]';
        no_loops=1000;
    end

    if nargin==3
        no_loops=1000;
    end
    unique_rep=unique(rep);
    no_reps=length(unique_rep); no_samples = size(data,1);
    h=waitbar(0, 'Please wait...');
    [m,c]=size(Y);
    if c==1
        unique_class = unique(Y);
        no_class = length(unique_class);
        Y_all = zeros(m, no_class);
        for i=1:no_class
            Y_all(Y==unique_class(i), i)=1;
        end
        [~, class_known] = max(Y_all,[],2);
    else
        no_class = c;
        Y_all = Y;
        [~, class_known] = max(Y_all,[],2);
    end
    
    for i=1:no_loops
        perm_idx = randperm(size(data,1));
        data2 = data(perm_idx,:);
        trn_rep=unique_rep(unique(randi(no_reps,no_reps,1)));
        trn_idx=find(ismember(rep, trn_rep));
        Ytrn=Y_all(trn_idx,:); class_trn = class_known(trn_idx);
        while length(trn_rep)<=4 || length(trn_rep)>=no_reps-2 || min(sum(Ytrn))==0
            trn_rep=unique_rep(unique(randi(no_reps,no_reps,1)));
            trn_idx=find(ismember(rep, trn_rep));
            Ytrn=Y_all(trn_idx,:); class_trn = class_known(trn_idx);
        end
        
        data_trn=data(trn_idx,:);
        data_trn2=data2(trn_idx,:);
        data_tst=data; data_tst(trn_idx,:)=[];
        data_tst2=data2; data_tst2(trn_idx,:)=[];
        
        Ytst=Y_all; Ytst(trn_idx,:)=[]; 
        class_tst=class_known; class_tst(trn_idx)=[];
        no_rep_trn=length(trn_rep);
        trn_rep_idx=rep(trn_idx);
        disp(['Current loop = ' num2str(i)])
        tic
        opt_parameter = dualplsda_tune(data_trn, Ytrn, 'rbf');
        opt_k = opt_parameter.opt_k; 
        opt_gamma = opt_parameter.opt_gamma;
        ccrcv_boots=opt_parameter.max_ccr_cv;

        Xcentre = mean(data_trn);        
        Xtrn = data_trn - repmat(Xcentre, size(data_trn,1),1);
        Xtst = data_tst - repmat(Xcentre, size(data_tst,1),1);
        Xtst2 = data_tst2 - repmat(Xcentre, size(data_tst,1),1);
        Ktr = rbf_two(Xtrn, Xtrn, opt_gamma);
        Kte = rbf_two(Xtrn, Xtst, opt_gamma);
        Kte2 = rbf_two(Xtrn, Xtst2, opt_gamma);
        [alpha, Yhat] = dualpls(Ktr, Kte, Ytrn, opt_k);
        [pred_raw, predC] = max(Yhat,[],2);
        [alpha, Yhat2] = dualpls(Ktr, Kte2, Ytrn, opt_k);
        [pred_raw2, predC2] = max(Yhat2,[],2);
        ccr_p = length(find(predC==class_tst))/length(class_tst);
        ccr_p2 = length(find(predC2==class_tst))/length(class_tst);
        for ii=1:no_class
            for iii=1:no_class
                conf_mat(ii,iii) = length(find(class_tst==ii & predC==iii))...
                    /length(find(class_tst==ii));
            end
        end
        for ii=1:no_class
            for iii=1:no_class
                conf_mat2(ii,iii) = length(find(class_tst==ii & predC2==iii))...
                    /length(find(class_tst==ii));
            end
        end
        output.Cknown{i}=class_tst;
        output.Cpred{i}=predC;
        output.ccr_cv(i,1)=ccrcv_boots;
        output.ccr(i,1)=ccr_p;
        output.ccr_perm(i)=ccr_p2;
        output.opt_k(i)=opt_k;
        output.opt_G(i)=opt_gamma;
        output.conf_mat(:,:,i) = conf_mat;
        output.conf_mat_perm(:,:,i) = conf_mat2;
        toc
        waitbar(i/no_loops,h);
    end
    close(h);
    conf_mat_avg = nanmean(output.conf_mat,3);
    output.conf_mat_avg = conf_mat_avg;
    output.ccr_avg = mean(output.ccr);
    figure
    h=image(conf_mat_avg);
    set(h,'CDatamapping','scaled')
    colorbar
    h=title('Averaged confusion matrix');
    set(h,'FontSize',14)
    set(gca,'XTick',1:no_class)
    set(gca,'YTick',1:no_class)
    
    disp(['Averaged CCR = ' num2str(output.ccr_avg)])
    disp('Averaged confusion matrix is:')
    disp(output.conf_mat_avg);    
    output = orderfields(output);
end


function H = hatchfill(A,STYL,ANGLE,SPACING,FACECOL)
% HATCHFILL Hatching and speckling of patch objects
%   HATCHFILL(A) fills the patch(es) with handle(s) A.
%   A can be a vector of handles or a single handle.
%   If A is a vector, then all objects of A should
%   be part of the same group for predictable results.
%   The hatch consists of black lines angled at
%   45 degrees spaced 5 pixels apart, with no color
%   filling between the lines.
%
%   HATCHFILL(A,STYL) applies STYL pattern with default paramters.
%      - STYL can be 'single' for single lines (the default),
%      'cross' for a double-crossed hatch, 'speckle' for
%     speckling inside the patch boundary, and 'outspeckle' for
%      for speckling outside the boundary. 'fill' will
%      apply only a gray fill and no hatching.
%
%   HATCHFILL(A,STYL,ANGLE,SPACING) applies a hatch/speckle with
%   customized parameters:
%      - ANGLE sets the angle of hatch lines. For speckling, it 
%      controls the width of the speckling region.
%      - SPACING controls the spacing of hatch lines or the
%      density of speckle points.
%      If STYL is 'fill', then ANGLE and SPACING are ignored.
%
%   HATCHFILL(A,STYL,ANGLE,SPACING,FACECOL) allows the user
%   to specify a fill color. (The default is 'none'.)
%
%   H = HATCHFILL(...) returns handles to the line objects
%   comprising the hatch/speckle.
%
%   Examples:
%       Gray region with hatching:
%       hh = hatchfill(a,'cross',45,5,[0.5 0.5 0.5]);
%
%       Speckled region:
%       hatchfill(a,'speckle',7,1);
%
%   NOTE: This function depends on the script hatch_xy.m
%   based on the work of R. Pawlowicz, K. Pankratov, and
%   Iram Weinstein.
%
%   Neil Tandon 11 Jul 2011

    % set defaults:
    if nargin == 1
        STYL = 'single';
        ANGLE = 45;
        SPACING = 5;
        FACECOL = 'none';
    end

    % For backwards compatability:
    if strcmpi(STYL,'none')
        STYL = 'fill';
    end

    if nargin == 2
        if strcmpi(STYL,'single') || strcmpi(STYL,'cross')
            ANGLE = 45;
            SPACING = 5;
            FACECOL = 'none';
        elseif strcmpi(STYL,'speckle') || strcmpi(STYL,'outspeckle')
            ANGLE = 7;
            SPACING = 1;
            FACECOL = 'none';
        elseif strcmpi(STYL,'fill')
            FACECOL = [0.8 0.8 0.8];
        end
    end

    if nargin == 3
        error('Invalid number of input arguments');
    end

    if nargin == 4
        if strcmpi(STYL,'fill')
            FACECOL = [0.8 0.8 0.8];
        else
            FACECOL = 'none';
        end
    end

    if ( ~strcmpi(STYL,'single') && ~strcmpi(STYL,'cross') && ...
         ~strcmpi(STYL,'speckle') && ~strcmpi(STYL,'outspeckle') && ...
         ~strcmpi(STYL,'fill') )
        error(['Invalid style: ',STYL])
    end

    linec = 'k';
    linew = 0.5;
    specksize = 2;

    % axis handle is one or two hierarchical levels up:
    % (Additional check suggested by Dan K)
    hax = get(A(1),'parent');
    is_axes = strcmpi(get(hax,'type'),'axes');
    if ~is_axes
       hax = get(hax,'parent');
    end
    is_axes = strcmpi(get(hax,'type'),'axes');

    x_is_log = 0; y_is_log = 0;
    x_is_reverse = 0; y_is_reverse = 0;

    if is_axes
       axsize_in = get(hax,'position');
       y_is_log = strcmpi(get(hax,'yscale'),'log');
       if y_is_log
           ylims = get(hax,'ylim');
           dy = (ylims(2) - ylims(1))/(log10(ylims(2))-log10(ylims(1)));
           set(hax,'units','pixels');
           axsize = get(hax,'position');
           set(hax,'position',[ axsize(1:3) dy*axsize(4) ]);
           set(hax,'units','normalized')
       end

       x_is_log = strcmpi(get(hax,'xscale'),'log');
       if x_is_log
           xlims = get(hax,'xlim');
           dx = (xlims(2) - xlims(1))/(log10(xlims(2))-log10(xlims(1)));
           set(hax,'units','pixels');
           axsize = get(hax,'position');
           set(hax,'position',[ axsize(1:2) dx*axsize(3) axsize(4) ]);
           set(hax,'units','normalized')
       end

       if strcmp(STYL,'single') || strcmp(STYL,'cross')
          y_is_reverse = strcmpi(get(hax,'ydir'),'reverse');
          if y_is_reverse
              ANGLE = -ANGLE;
          end
          x_is_reverse = strcmpi(get(hax,'xdir'),'reverse');
          if x_is_reverse
              ANGLE = 180-ANGLE;
          end
       end
    end

    % Apply hatch:
    j = 1;
    for k = 1:length(A)
        set(A,'facecolor',FACECOL);
        v = get(A(k),'vertices');
        if any(v(end,:)~=v(1,:))
            v(end+1,:) = v(1,:);
        end
        x = v(:,1);
        if x_is_log
            x = log10(v(:,1));
        end
        y = v(:,2);
        if y_is_log
            y = log10(v(:,2));
        end

        if strcmp(STYL,'fill')
            H = NaN;
            continue
        end

        [xhatch,yhatch] = hatch_xy(x,y,STYL,ANGLE,SPACING);
        if x_is_log
            xhatch = 10.^xhatch;
        end
        if y_is_log
            yhatch = 10.^yhatch;
        end
        if strcmp(STYL,'speckle') || strcmp(STYL,'outspeckle')
            if any(xhatch)
                H(j) = line(xhatch,yhatch,'marker','.','linest','none', ...
                        'markersize',specksize,'color',linec);
                j = j+1;
            end
        elseif strcmp(STYL,'single') || strcmp(STYL,'cross')
            H(j) = line(xhatch,yhatch);
            set(H(j),'color',linec,'linewidth',linew);
            j = j+1;
        end
    end

    if y_is_log || x_is_log
        set(hax,'position',axsize_in);
    end
end

%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%

function [xi,yi,x,y]=hatch_xy(x,y,varargin);
%
% M_HATCH Draws hatched or speckled interiors to a patch
%       
%    M_HATCH(LON,LAT,STYL,ANGLE,STEP,...line parameters);
%
% INPUTS:
%     X,Y - vectors of points.
%     STYL - style of fill
%     ANGLE,STEP - parameters for style
%
%     E.g.
%                 
%      'single',45,5  - single cross-hatch, 45 degrees,  5 points apart 
%      'cross',40,6   - double cross-hatch at 40 and 90+40, 6 points apart
%      'speckle',7,1  - speckled (inside) boundary of width 7 points, density 1
%                               (density >0, .1 dense 1 OK, 5 sparse)
%      'outspeckle',7,1 - speckled (outside) boundary of width 7 points, density 1
%                               (density >0, .1 dense 1 OK, 5 sparse)
%
%     
%      H=M_HATCH(...) returns handles to hatches/speckles.
%
%      [XI,YI,X,Y]=MHATCH(...) does not draw lines - instead it returns
%      vectors XI,YI of the hatch/speckle info, and X,Y of the original
%      outline modified so the first point==last point (if necessary).
%
%     Note that inside and outside speckling are done quite differently
%     and 'outside' speckling on large coastlines can be very slow.

%
% Hatch Algorithm originally by K. Pankratov, with a bit stolen from 
% Iram Weinsteins 'fancification'. Speckle modifications by R. Pawlowicz.
%
% R Pawlowicz 15/Dec/2005
  
    styl='speckle';
    angle=7;
    step=1/2;

    if length(varargin)>0 & isstr(varargin{1}),
      styl=varargin{1};
      varargin(1)=[];  
    end;
    if length(varargin)>0 & ~isstr(varargin{1}),
      angle=varargin{1};
      varargin(1)=[];  
    end;
    if length(varargin)>0 & ~isstr(varargin{1}),
      step=varargin{1};
      varargin(1)=[];
    end;

    I = zeros(1,length(x));
    %[x,y,I]=m_ll2xy(lon,lat,'clip','patch');


    if x(end)~=x(1) & y(end)~=y(1),
      x=x([1:end 1]);
      y=y([1:end 1]);
      I=I([1:end 1]);
    end;

    if strcmp(styl,'speckle') | strcmp(styl,'outspeckle'),
      angle=angle*(1-I);
    end;

    if size(x,1)~=1,
     x=x(:)';
     angle=angle(:)';
    end;
    if size(y,1)~=1,
     y=y(:)';
    end;


    % Code stolen from Weinstein hatch
    oldu = get(gca,'units');
    set(gca,'units','points');
    sza = get(gca,'pos'); sza = sza(3:4);
    set(gca,'units',oldu)   % Set axes units back

    xlim = get(gca,'xlim');
    ylim = get(gca,'ylim');
    xsc = sza(1)/(xlim(2)-xlim(1)+eps);
    ysc = sza(2)/(ylim(2)-ylim(1)+eps);

    switch lower(styl),
     case 'single',
      [xi,yi]=drawhatch(x,y,angle,step,xsc,ysc,0);
      if nargout<2,
        xi=line(xi,yi,varargin{:});
      end;  
     case 'cross',
      [xi,yi]=drawhatch(x,y,angle,step,xsc,ysc,0);
      [xi2,yi2]=drawhatch(x,y,angle+90,step,xsc,ysc,0);
      xi=[xi,xi2];
      yi=[yi,yi2];
      if nargout<2,
        xi=line(xi,yi,varargin{:});
      end;  
     case 'speckle',
      [xi,yi ]  =drawhatch(x,y,45,   step,xsc,ysc,angle);
      [xi2,yi2 ]=drawhatch(x,y,45+90,step,xsc,ysc,angle);
      xi=[xi,xi2];
      yi=[yi,yi2];
      if nargout<2,
        if any(xi),
          xi=line(xi,yi,'marker','.','linest','none','markersize',2,varargin{:});
        else
          xi=NaN;
        end;    
      end; 
     case 'outspeckle',
      [xi,yi ]  =drawhatch(x,y,45,   step,xsc,ysc,-angle);
      [xi2,yi2 ]=drawhatch(x,y,45+90,step,xsc,ysc,-angle);
      xi=[xi,xi2];
      yi=[yi,yi2];
      inside=logical(inpolygon(xi,yi,x,y)); % logical needed for v6!
      xi(inside)=[];yi(inside)=[];
      if nargout<2,
        if any(xi),
          xi=line(xi,yi,'marker','.','linest','none','markersize',2,varargin{:});
        else
          xi=NaN;
        end;    
      end; 

    end;


    return
end

function [xi,yi]=drawhatch(x,y,angle,step,xsc,ysc,speckle);
%
% This is the guts. 
%

    angle=angle*pi/180;

    % Idea here appears to be to rotate everthing so lines will be
    % horizontal, and scaled so we go in integer steps in 'y' with
    % 'points' being the units in x.
    % Center it for "good behavior".
    ca = cos(angle); sa = sin(angle);
    x0 = mean(x); y0 = mean(y);   
    x = (x-x0)*xsc; y = (y-y0)*ysc;
    yi = x*ca+y*sa;              % Rotation
    y = -x*sa+y*ca;
    x = yi;
    y = y/step;    % Make steps equal to one

    % Compute the coordinates of the hatch line ...............
    yi = ceil(y);
    yd = [diff(yi) 0]; % when diff~=0 we are crossing an integer
    fnd = find(yd);    % indices of crossings
    dm = max(abs(yd)); % max possible #of integers between points


    %
    % This is going to be pretty space-inefficient if the line segments
    % going in have very different lengths. We have one column per line
    % interval and one row per hatch line within that interval.
    %
    A = cumsum( repmat(sign(yd(fnd)),dm,1), 1);

    % Here we interpolate points along all the line segments at the
    % correct intervals.
    fnd1 = find(abs(A)<=abs( repmat(yd(fnd),dm,1) ));
    A  = A+repmat(yi(fnd),dm,1)-(A>0);
    xy = (x(fnd+1)-x(fnd))./(y(fnd+1)-y(fnd));
    xi = repmat(x(fnd),dm,1)+(A-repmat(y(fnd),dm,1) ).*repmat(xy,dm,1);
    yi = A(fnd1);
    xi = xi(fnd1);


     % Sorting points of the hatch line ........................
    %%%yi0 = min(yi); yi1 = max(yi);
    % Sort them in raster order (i.e. by x, then by y)
    % Add '2' to make sure we don't have problems going from a max(xi)
    % to a min(xi) on the next line (yi incremented by one)
    xi0 = min(xi); xi1 = max(xi);
    ci = 2*yi*(xi1-xi0)+xi;
    [ci,num] = sort(ci);
    xi = xi(num); yi = yi(num);


    % if this happens an error has occurred somewhere (we have an odd
    % # of points), and the "fix" is not correct, but for speckling anyway
    % it really doesn't make a difference.
    if rem(length(xi),2)==1, 
      disp('mhatch warning');
      xi = [xi; xi(end)];
      yi = [yi; yi(end)];
    end

     % Organize to pairs and separate by  NaN's ................
    li = length(xi);
    xi = reshape(xi,2,li/2);
    yi = reshape(yi,2,li/2);

    % The speckly part - instead of taking the line we make a point some
    % random distance in.
    if length(speckle)>1 | speckle(1)~=0,

     if length(speckle)>1,
       % Now we get the speckle parameter for each line.

       % First, carry over the speckle parameter for the segment
    %   yd=[0 speckle(1:end-1)];
       yd=[speckle(1:end)];
       A=repmat(yd(fnd),dm,1);
       speckle=A(fnd1);

       % Now give it the same preconditioning as for xi/yi
       speckle=speckle(num);
       if rem(length(speckle),2)==1, 
         speckle = [speckle; speckle(end)];
       end
       speckle=reshape(speckle,2,li/2);

     else
       speckle=[speckle;speckle];
     end;

     % Thin out the points in narrow parts.
     % This keeps everything when abs(dxi)>2*speckle, and then makes
     % it increasingly sparse for smaller intervals.
     oldxi=xi;oldyi=yi;
     dxi=diff(xi);
     nottoosmall=sum(speckle,1)~=0 & rand(1,li/2)<abs(dxi)./(max(sum(speckle,1),eps));
     xi=xi(:,nottoosmall);
     yi=yi(:,nottoosmall);
     dxi=dxi(nottoosmall);
     if size(speckle,2)>1, speckle=speckle(:,nottoosmall); end;
     % Now randomly scatter points (if there any left)
     li=length(dxi);
     if any(li),
       xi(1,:)=xi(1,:)+sign(dxi).*(1-rand(1,li).^0.5).*min(speckle(1,:),abs(dxi) );
       xi(2,:)=xi(2,:)-sign(dxi).*(1-rand(1,li).^0.5).*min(speckle(2,:),abs(dxi) );
       % Remove the 'zero' speckles
       if size(speckle,2)>1,
        xi=xi(speckle~=0);
        yi=yi(speckle~=0);
       end;
      end;

    else
     xi = [xi; ones(1,li/2)*nan];  % Separate the line segments
     yi = [yi; ones(1,li/2)*nan];
    end;
    xi = xi(:)'; yi = yi(:)';

    % Transform back to the original coordinate system
    yi = yi*step;
    xy = xi*ca-yi*sa;
    yi = xi*sa+yi*ca;
    xi = xy/xsc+x0;
    yi = yi/ysc+y0;
end
