function W = WavmatND(hf, N, k, shift)
% WavmatND -- Transformation Matrix for Non-Decimated WT
%  Usage
%    W = WavmatND(h, N, k0, shift)
%  Inputs
%    hf   --   low-pass filter corresponding to orthogonal WT
%    N  --    size of matrix/length of data. Should be power of 2.      
%    k  --   depth of transformation. Range >= 1. Although  not
% limited k0 > log2(N) leads to oversmoothing and
% high influence of the filter
%      
%    shift -- the matrix is not unique and any integer shift gives
%           a valid transformation. Default is 0.
%  Outputs
%    W      k0 * N x N transformation matrix 
%
%  Description
%    For a quadrature mirror filter h (low pass) the wavelet
%    matrix is formed to perform Non-decimated wavelet transform.  
%
%  Usage
%   > dat=[1 0 -3 2 1 0 1 2];
%   > W = WavmatND(MakeONFilter('Haar',1),2^3,3,2);
%   > wt = W * dat' % to obtain a transformed signal
%   > data = W' * wt % should return you to the 'dat'
%
%  See Also
%  FWT_PO, IWT_PO, MakeONFilter
%

    hf=hf(:)';  gf = fliplr(conj(hf).* (-1).^(0:length(hf)-1 ));
    W=[]; hmatold=eye(N);
    h=[hf,zeros(1,N)]; %extend filter H by 0's to sample by modulus
    g=[gf,zeros(1,N)]; %extend filter G by 0's to sample by modulus
for i = 1:k
    clear gmat; clear hmat;
    for  jj= 1:N
       for ii=1:N
           modulus = mod(N + ii - jj - shift , N) + 1;
           modulus = modulus + (modulus == 0)*N;   
           hmat(ii,jj) =  h(modulus);
           gmat(ii,jj) = g(modulus);
       end
    end
%-------------------------------------
    W=[ gmat' * hmatold  ; W];
    smooth =hmat'* hmatold ;
    hmatold = smooth; 
    h = [dilate_filter(hf,2^(i)-1),zeros(1,N)];
    g = [dilate_filter(gf,2^(i)-1),zeros(1,N)];
end

W=[smooth; W];
end
%--------------------------------------
function  filtd = dilate_filter(filt,k)
%-------------------------------------
%  Dilate Filter by k zeros between the original 
% taps, k integer > 0.
%-------------------------------------
newlength = (k+1)*length(filt)-k;
filtd = zeros(1,newlength);
filtd(1:(k+1):newlength)=filt;
%--------------------------------------
% Brani 10/15/02
end
function out=imbed(a)
% imbeds 0 in the filter a
b=[a; repeat(0, length(a))];
out=b(:)';
end
function b = repeat(a, n)
% Repeats an array a n times
% Usage
%   b = repeat(a, n)
% Input
%   a, n
% Output
%   b
% Brani 10/20/02

b=[];
for i=1:n
b=[b a];
end
end


