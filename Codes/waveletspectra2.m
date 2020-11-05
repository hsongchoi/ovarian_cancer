function [sloped, slopeh, slopev, levels, log2specd, log2spech, log2specv ] = ...
                    waveletspectra2(data, L, wf, k1, k2)
%
%  [sloped, slopeh, slopev, levels, log2specd, log2spech, log2specv ] = ...
%                    waveletspectra2(data, L, wf, k1, k2)
%  input:  data - data in time domain
%          L - coarse level.
%          wf - wavelet filter
%          k1 -  start with coarse level k1 when calculating slope, k1 >= L.
%          k2 -  end with the level  k2 when calculating slope, k2<=log2(n)-1
%
%  output: slope - scaling slope of log2-energies.
%          levels - integers L, L+1, ..., log2(n)-1
%          log2spec - log2 of levelwise averages of squared wavelet
%                     coefficients 
%
% 
if nargin == 1,  L=1;  wf=[sqrt(2)/2 sqrt(2)/2];  k1=1; k2=log2(length(data))-1 ;  end
if nargin == 2,        wf=[sqrt(2)/2 sqrt(2)/2];  k1=L; k2=log2(length(data))-1 ;  end
if nargin == 3,                                   k1=L; k2=log2(length(data))-1 ;  end
if nargin == 4,                                         k2=log2(length(data))-1 ;  end


siz = size(data);
lnn = log2( siz(1) );
%wddata = FWT2_PO(data, L, wf); 
wddata = dwtrn(data, lnn - L, wf); 

yd = [ ]; 
yh = [ ];
yv = [ ];
for i =  L:(lnn-1)
   [cxd, cyd] = dyad2(i,'d');
   [cxh, cyh] = dyad2(i,'h');
   [cxv, cyv] = dyad2(i,'v');
   helpd = wddata(cxd,cyd);
   helph = wddata(cxh,cyh);
   helpv = wddata(cxv,cyv);
          yd = [ yd  mean(mean(helpd.^2))  ];
          yh = [ yh  mean(mean(helph.^2)) ];  
          yv = [ yv  mean(mean(helpv.^2)) ];
end
   levels = L:(lnn-1);
    log2specd = log2(yd);
    log2spech = log2(yh);
    log2specv = log2(yv);
   yyd = log2specd(k1-L+1:k2-L+1);
   yyh = log2spech(k1-L+1:k2-L+1);
   yyv = log2specv(k1-L+1:k2-L+1);
   aad =polyfit([k1:k2], yyd, 1);
    aah =polyfit([k1:k2], yyh, 1);
     aav =polyfit([k1:k2], yyv, 1);
   sloped = aad(1); slopeh = aah(1); slopev = aav(1);
   ccd = polyval(aad,[k1:k2]);
    cch = polyval(aah,[k1:k2]);
     ccv = polyval(aav,[k1:k2]);
       %--- set plotting parameters -------
        lw = 2;  
        set(0, 'DefaultAxesFontSize', 15);
        fs = 15;
        msize = 6;
        plot(levels, log2specd,  'linewidth', lw)
        hold on
        plot(levels, log2spech,  'linewidth', lw)
        plot(levels, log2specv,  'linewidth', lw)
        
        plot(levels, log2specd, 'd', 'markersize', msize)
        plot(levels, log2spech, 'o', 'markersize', msize)
        plot(levels, log2specv, '*', 'markersize', msize)
         
        plot(k1:k2, yyd, 'r-', 'linewidth', lw)
        plot(k1:k2, yyh, 'r-', 'linewidth', lw)
        plot(k1:k2, yyv, 'r-', 'linewidth', lw)
          
          
        plot(k1:k2, ccd + 1,'g:','linewidth', lw)
        plot(k1:k2, cch + 1,'g:','linewidth', lw)
        plot(k1:k2, ccv + 1,'g:','linewidth', lw)
          
        text( k1,ccd(1)+2.5, num2cell(sloped) )
        text( k1,ccd(1)+1.5, num2cell(slopeh) )
        text( k1,ccd(1)+0.5, num2cell(slopev) )
          
        xlabel('dyadic level','fontweight','bold','fontsize',fs)
        ylabel('log spectrum','fontweight','bold','fontsize',fs)
        axis tight
        hold off
       
%-------------- Brani 10/06-------------------------------------------