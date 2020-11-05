function [slope, levels, log2spectt, W, wddata ] = spectramix2d(data, L, wf, k1, k2, displace, ismean, q, isGraph)
%
%  [slope, levels, log2spectt, W, wddata ] = spectramix2d(data, L, wf, k1, k2, displace,
%  ismean, q, isGraph)
%  input:  data - data in time domain
%          L - coarse level.
%          wf - wavelet filter
%          k1 -  start with coarse level k1 when calculating slope, k1 >= L.
%          k2 -  end with the level  k2 when calculating slope, k2<=log2(n)-1
%          displace - the number of displacement steps from the diagonal
%                       It should satisfy k1 >= L + abs(displace).
%                       By default, 0.
%                       0 means the diagoal
%                       1 means the one-step upper diagonal
%                       -1 means the one-step lower diagonal
%          ismean - 0 for median and 1 for mean, default is mean.
%          q - the exponent of absolute value of wavelet coefficients
%                       By default, 2.
%          isGraph - 1 for drawing picture, 0 for not drawing.
%                       By default, 1.
%  output: slope - scaling slope of log2-energies.
%          levels - integers L, L+1, ..., log2(n)-1
%          log2spectt - log2 of levelwise averages of squared wavelet
%                     coefficients
%          W - the matrix for wavelet transformation
%          wddata - the result of separate wavelet transformation;
%                   W*data*W'
%       
%  example: %2D Stick
%           % parameters
%           lnn = 10; wf = MakeONFilter('Symmlet', 6); L = 3;
%           k1 = 5; k2 = 8; displace = 1;
%           ismean = 1; q = 2; isGraph = 1;
%           ImData = Make2dSignal('StickFigure',2^lnn);
%
%           [slope, levels, log2spectt, W, wddata ] = spectramix2d(ImData, L, wf, k1, k2, displace, ismean, q, isGraph);
%
%

if nargin == 1,  L=1;  wf=[sqrt(2)/2 sqrt(2)/2];  k1=1; k2=log2(length(data))-1;  displace = 0; ismean=1; q=2; isGraph=1; end
if nargin == 2,        wf=[sqrt(2)/2 sqrt(2)/2];  k1=L; k2=log2(length(data))-1 ; displace = 0; ismean=1; q=2; isGraph=1; end
if nargin == 3,                                   k1=L; k2=log2(length(data))-1 ; displace = 0; ismean=1; q=2; isGraph=1; end
if nargin == 4,                                         k2=log2(length(data))-1 ; displace = 0; ismean=1; q=2; isGraph=1; end
if nargin == 5,                                                                   displace = 0; ismean=1; q=2; isGraph=1; end
if nargin == 6,                                                                                 ismean=1; q=2; isGraph=1; end
if nargin == 7,                                                                                           q=2; isGraph=1; end
if nargin == 8,                                                                                                isGraph=1; end

if ~( k1 >= L + abs(displace))
    disp('Should satisfy k1 >= L + abs(displace)');
    return;
end

Diff = displace;
lnn = log2(length(data));
%wddata = FWT_PO(data, L, wf);
%wddata = dwtr(data, lnn - L, wf);
W  = Wavmat(wf,length(data),lnn - L);
wddata  = W*data*W';

y = [];
%for i =  L:(lnn-1)
for i =  L+abs(Diff):(lnn-1)
    if displace >= 0
        [rIdx cIdx] = dyadupper2D(i, Diff, lnn);
    else
        [rIdx cIdx] = dyadlower2D(i, abs(Diff), lnn);
    end
    %help = wddata((2^(i)+1):(2^(i+1)));
    help = wddata(rIdx, cIdx);
    help = help(:);
    if ismean == 1
        y = [y mean(abs(help).^q)];
    elseif ismean == 0
        y = [y median(abs(help).^q)];
    else error('not known average of energies. Use mean (ismean=1) or median (ismean=0)')
    end
end

%levels = L:(lnn-1);
levels = L+abs(Diff):(lnn-1);
log2spectt = log2(y);
%yy = log2spectt(k1-L+1:k2-L+1);
yy = log2spectt(k1-L-abs(Diff)+1:k2-L+1-abs(Diff));

aa =polyfit([k1:k2], yy, 1);
slope = aa(1);
cc = polyval(aa,[k1:k2]);

if isGraph == 1
    %--- set plotting parameters -------
    lw = 2;
    set(0, 'DefaultAxesFontSize', 15);
    fs = 15;
    msize = 6;
    plot(levels, log2spectt,  'linewidth', lw)
    hold on
    plot(levels, log2spectt, 'o', 'markersize', msize)
    plot(k1:k2, yy, 'r-', 'linewidth', lw);
    plot(k1:k2, cc + 1,'g:','linewidth', lw);
    log2spectt - log2(y);
    % plot(k3:k4, yy2, 'r-', 'linewidth', lw);
    % plot(k3:k4, cc2 + 1,'g:','linewidth', lw);
    % text( k3+0.2,cc2(1)+0.6, num2cell(slope2), ...
    %     'FontName', 'Times', 'FontWeight', 'bold', 'FontSize', 9);

    text( k1+0.2,cc(1)+0.6, num2cell(slope), ...
        'FontName', 'Times', 'FontWeight', 'bold', 'FontSize', fs+3);

    xlabel('dyadic level','fontweight','bold','fontsize',fs)
    ylabel('log spectrum','fontweight','bold','fontsize',fs)
    axis tight
    hold off
end
%-------------- Brani 10/06-------------------------------------------
function W = Wavmat(h, N, k0, shift)
% WavMat -- Transformation Matrix of FWT_PO
%  Usage
%    W = WavMat(h, N, k0, shift)
%  Inputs
%    h      low-pass filter corresponding to orthogonal WT
%    N      size of matrix/length of data. Should be power of 2.
%      
%    k0     depth of transformation. Ranges from 1 to J=log2(N).
%           Default is J. 
%    shift  the matrix is not unique an any integer shift gives
%           a valid transformation. Default is 2.
%  Outputs
%    W      N x N transformation matrix 
%
%  Description
%    For a quadrature mirror filter h (low pass) the wavelet
%    matrix is formed. The algorithm is described in 
%    [BV99] Vidakovic, B. (1999). Statistical Modeling By Wavelets, Wiley,
%    on pages 115-116.
%    Any shift is valid.  Size N=1024 is still managable on a standard PC.
%
%  Usage
%    We will mimic the example 4.3.1 from [BV99] page 112.
%   > dat=[1 0 -3 2 1 0 1 2];
%   > W = WavMat(MakeONFilter('Haar',99),2^3,3,2);
%   > wt = W * dat' %should be [sqrt(2)  |  -sqrt(2) |   1 -1  | ...         
%              %  1/sqrt(2) -5/sqrt(2) 1/sqrt(2) - 1/sqrt(2) ]
%   > data = W' * wt % should return you to the 'dat'
%
%    > dat = [1 0 -3 2 1 0 1 2];
%    > filter = [sqrt(2)/2 sqrt(2)/2];
%    > W  = WavMat(filter,2^3,3);
%    > wt = W * dat' %should be [sqrt(2)  |  -sqrt(2) |   1 -1  | ...         
%                1/sqrt(2) -5/sqrt(2) 1/sqrt(2) - 1/sqrt(2) ]
%    > data = W' * wt % should return you to the 'dat'
%
%  See Also
%    FWT_PO, IWT_PO, MakeONFilter
%
if nargin==3
    shift = 2;
end
J = log2(N);
if nargin==2
    shift = 2;
    k0 = J;
end
%--make QM filter G
     h=h(:)';  g = fliplr(h .* (-1).^(1:length(h)));
if (J ~= floor(J) )
    error('N has to be a power of 2.')
end
h=[h,zeros(1,N)]; %extend filter H by 0's to sample by modulus
g=[g,zeros(1,N)]; %extend filter G by 0's to sample by modulus
oldmat = eye(2^(J-k0)); 
for k= k0:-1:1
    clear gmat; clear hmat;
         ubJk = 2^(J-k); ubJk1 = 2^(J-k+1);
   for  jj= 1:ubJk
       for ii=1:ubJk1
           modulus = mod(N+ii-2*jj+shift,ubJk1);
           modulus = modulus + (modulus == 0)*ubJk1;
           hmat(ii,jj) = h(modulus);
           gmat(ii,jj) = g(modulus);
       end
   end
  W = [oldmat * hmat'; gmat' ];
   oldmat = W;
end
%
% 
% Copyright (c) 2004. Brani Vidakovic
%        
%  
% ver 1.0 Built 8/24/04; ver 1.2 Built 12/1/2004
% This is Copyrighted Material
% Comments? e-mail brani@isye.gatech.edu

function [vIdx hIdx] = dyadupper2D(j, Diff, log2len, L)

%
% Input
%   j : the resolution of index
%   Diff : difference between vertical resolution and 
%           horizontal resolution of the starting position;
%           starting at the diagonal position means Diff = 0;
%           starting at one position upper above the diagonal means Diff=1;
%   log2len : the finest resolution of the signal
%       it is set to be defined to be log2len = log2(length of signal);
%   L : the coarsest level of the wavelet transformation
%       it is set to be 1 by default.
%
% Output
%   vIdx : index for vertical direction;
%           corresponding to row index in 2-D
%   hIdx : index for horizontal direction;
%           corresponding to column index in 2-D
% Usage
%   >>[vIdx hIdx] = dyadupper2D(8, 2, 10)
%   >>[vIdx hIdx] = dyadupper2D(3, 1, 5)
%
% Version
%   Sky Lee, Dec-29-2008, ver 1 - initial implementation
%

if nargin < 4
    L = 1;
end
if ~(j >= L+Diff)
    disp('should satisfy j >= L+Diff');
    return;
end
    
if ~(log2len-1-L >= Diff)
    disp('should satisfy log2len-1-L >= Diff');
    return;
end
%k = log2len-j;
vIdx = (2^(j-Diff)+1:2^(j-Diff+1))';
hIdx = (2^(j)+1:2^(j+1))';


function [vIdx hIdx] = dyadlower2D(j, Diff, log2len, L)

%
% Input
%   j : the resolution of index
%   Diff : difference between vertical resolution and 
%           horizontal resolution of the starting position;
%           starting at the diagonal position means Diff = 0;
%           starting at one position upper above the diagonal means Diff=1;
%   log2len : the finest resolution of the signal
%       it is set to be defined to be J = log2(length of signal);
%   L : the coarsest level of the wavelet transformation
%       it is set to be 1 by default.
%
% Output
%   vIdx : index for vertical direction;
%           corresponding to row index in 2-D
%   hIdx : index for horizontal direction;
%           corresponding to column index in 2-D
% Usage
%   >>[vIdx hIdx] = dyadlower2D(8, 2, 10)
%   >>[vIdx hIdx] = dyadlower2D(3, 1, 5)
%
% See also
%   dyadupper
% Version
%   Sky Lee, Dec-29-2008, ver 1 - initial implementation
%
if nargin < 4
    L = 1;
end
[hIdx vIdx] = dyadupper2D(j, Diff, log2len, L);



