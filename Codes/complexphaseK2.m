function [avg_1 variance_1 ske_1 kurt_1 Q1_1 med_1 Q3_1 CV_1 MADs_1] = complexphaseK2(data, wf)

N=6; confilt=MakeCONFilter(N);

if nargin == 1,  wf=confilt; end

[I J]=size(data); ln=log2(I);
Diff=0;

%2D covariance-form complex wavelet transformation
W  = WavMat(wf,length(data),ln-1);
wddata  = W*data*W';

[rIdx cIdx] = dyadupper2D(ln-1, Diff, ln);
phase_finest=angle(wddata(rIdx, cIdx)); 

%[cxd, cyd] = dyad2(ln-1,'d'); phase_finest=angle(wddata(cxd,cyd));
   
avg_1      = mean(phase_finest(:)); 
variance_1 = var(phase_finest(:));
ske_1      = skewness(phase_finest(:)); 
kurt_1     = kurtosis(phase_finest(:));
Q1_1       = quantile(phase_finest(:),.25); 
med_1      = quantile(phase_finest(:),.50);
Q3_1       = quantile(phase_finest(:),.75);
CV_1       = nanstd(phase_finest(:))/nanmean(phase_finest(:));%*100;
MADs_1     = mad(phase_finest(:));

