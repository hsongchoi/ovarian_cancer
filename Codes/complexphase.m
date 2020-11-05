function [phase modulus] = complexphase(data, wf)

N=6; confilt=MakeCONFilter(N);

if nargin == 1,  wf=confilt; end

[I J]=size(data); ln=log2(I);
Diff=0;

%2D covariance-form complex wavelet transformation
% W  = Wavmat(wf,length(data),ln-1);
% wddata  = W*data;
wddata = dwtr(data, ln-1, wf); 

%[rIdx cIdx] = dyadupper2D(ln-1, Diff, ln);
phase=angle(wddata.^2); 
modulus=abs(wddata);

%[cxd, cyd] = dyad2(ln-1,'d'); phase_finest=angle(wddata(cxd,cyd));
   
% avg_1      = mean(phase_finest(:)); 
% variance_1 = var(phase_finest(:));
% ske_1      = skewness(phase_finest(:)); 
% kurt_1     = kurtosis(phase_finest(:));
% Q1_1       = quantile(phase_finest(:),.25); 
% med_1      = quantile(phase_finest(:),.50);
% Q3_1       = quantile(phase_finest(:),.75);
% CV_1       = nanstd(phase_finest(:))/nanmean(phase_finest(:));%*100;
% MADs_1     = mad(phase_finest(:));

