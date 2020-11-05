clear 
clc
close all

% Take all data of protein spectra from csv files
% Sample 162 cases and 91 controls.
% break interval approximately [2000,20000] into each subinterval
% [xi, xi+1023]

N=6; filt=MakeCONFilter(N); %complex filter 
%filt=[sqrt(2)/2 sqrt(2)/2]; % haar wavelet
% filt = [ -0.07576571478934  -0.02963552764595  ...
%           0.49761866763246   0.80373875180522  ...
%           0.29785779560554  -0.09921954357694  ...
%          -0.01260396726226   0.03222310060407]; 
d = dir(fullfile('*.csv'));
dataset = [];
phasespecMatirx = [];
modhurstMatirx = [];

for caseNumber = 1 : 253
    filename = fullfile(d(caseNumber).name);
    data = csvread(filename,1,0);
    dataset(:,:,caseNumber) = data;
    column = 1;
    column2 = 1;
    
    for step = 4001 : 100 : 13001
        
        y = data(step:step+1023,2);
        [phaselevels, phasespec] = complexwaveletspectraphase(y,1,filt,5,9,1);
        [modslope, modlevels, modlog2spec] = complexwaveletspectramodulus(y,1,filt,5,9,1);
        phasespecMatirx(caseNumber,column2:column2+8) = phasespec;
        h = -1.*modslope./2 - 1/2;
        modhurstMatrix(caseNumber,column) = h;
        column = column + 1;
        column2 = column2 + 9;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% 
% load('hurstMatrix.mat')
% 
% % divde control and cases
% % Get means, medians and standard deviation of 
% % Matrix A(control) and Matirx B(case)
% control = hurstMatrix(1:91,:); cases = hurstMatrix(92:253,:);
% control_data= [mean(control); median(control); std(control)];
% cases_data = [mean(cases); median(cases); std(cases)];
% diff = control_data - cases_data;
% 
% b = glmfit(hA(:,1),hB(:,1)); % b = glmfit(x,y,distr)
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % how the hurst exponent looks like between A and B.
% figure(1)
% 
% subplot(2,2,1)
% plot(control_data(1,:))
% title('Hurst exponents of Control')
% 
% subplot(2,2,2)
% plot(cases_data(1,:))
% title('Hurst exponents of Case')
% 
% subplot(2,2,3)
% plot(control_data(1,:)-cases_data(1,:))
% title('Difference betwween Control and case')
% 
% subplot(2,2,4)
% plot(abs(control_data(1,:)-cases_data(1,:)))
% title('Absolute Difference betwween Control and case')
% 
% 
% %
% figure(2)
% 
% subplot(2,2,1)
% x = 1: 901;
% y = control_data(1,:);
% err = control_data(3,:);
% errorbar(x,y,err);
% title('error bar of Control')
% 
% subplot(2,2,2)
% y = cases_data(1,:);
% err = cases_data(3,:);
% errorbar(x,y,err);
% title('error bar of Case')
% 
% subplot(2,2,3)
% y = control_data(1,:)-cases_data(1,:);
% err = control_data(3,:) + cases_data(3,:);
% errorbar(x,y,err);
% title('error bar of Difference between Control and Case')
