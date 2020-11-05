function HurstExponent

clear all

% 1. Wavelet transform
% Take all data of protein spectra from csv files
% Sample 162 cases and 91 controls.
% break interval approximately [2000,20000] into each subinterval
% [xi, xi+1023]

filt=[sqrt(2)/2 sqrt(2)/2]; % haar wavelet
d = dir(fullfile('*.csv'));
dataset = [];
hurstMatirx = [];

for caseNumber = 1 : length(d)
    filename = fullfile(d(caseNumber).name);
    data = csvread(filename,1,0);
    dataset(:,:,caseNumber) = data;
    column = 1;
    
    for step = 4001 : 10 : 13001
        
        y = data(step:step+1023,2);
        [slope, ~,~] = waveletspectra(y,1,filt,5,9,1);
        h = -1.*slope./2 - 1/2;
        hurstMatrix(caseNumber,column) = h;
        step = step + 10;
        column = column + 1;
    end
end

% divde control and cases
% Get means, medians and standard deviation of 
% Matrix A(control) and Matirx B(case)

A = hurstMatrix(1:91,:); B = hurstMatrix(92:253,:);
hA = [mean(A)', median(A)', std(A)'];
hB = [mean(B)', median(B)', std(B)'];

% how the hurst exponent looks like between A and B.

figure(1)

subplot(2,2,1)
plot(hA(:,1))

subplot(2,2,2)
plot(hB(:,1))

subplot(2,2,3)
plot(abs(hA(:,1)-hB(:,1)))

% 2. logistic regression

b = glmfit(hA(:,1),hB(:,1)); % b = glmfit(x,y,distr)
b0 = b(1); b1 = b(2); %evaluate b0 and b1%evaluate b0 and b1
tp = 0; fp = 0; fn = 0; tn = 0;

for i = 1 : length(evaluant)
    H = intensity(evaluant(i),2);
    linpart = b0 + b1*H;
    pp=exp(linpart)./(1+exp(linpart));
    filename = d(evaluant(i)).name;
    if pp >= prediction; % Test CANCER
        if 'O' == filename(1) % True Positive
            tp = tp + 1;
        else % False Positive
            fp = fp + 1;
        end
    else
        if 'O' == filename(1) % false negative
            fn = fn + 1;
        else % True Negative
            tn = tn + 1;
        end
    end
end

confusionMatrix = [tp,fp;fn,tn];

[se sp pre ppv npv lrp accuracy yi] = sesp(tp, fp, fn, tn);

% how the hurst exponent looks like between A and B.

figure(1)

subplot(2,2,1)
plot(hA(:,1))

subplot(2,2,2)
plot(hB(:,1))

subplot(2,2,3)
