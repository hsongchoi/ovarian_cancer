function [confusionMatrix accuracy se sp]=  Ovarian(prediction)

% load 'symm4.mat'
load('hurstMatrix4.mat')
%load('hurstMatrix.mat')

% divde control and cases
% Get means, medians and standard deviation of control and case

control = hurstMatrix(1:91,:); cases = hurstMatrix(92:253,:);
control_data= [mean(control); median(control); std(control)];
cases_data = [mean(cases); median(cases); std(cases)];
diff = abs(control_data(1,:) - cases_data(1,:));
[B,I] = sort(diff,'descend');
h1 = B(1); h2 = B(2); I1 = I(1); I2 = I(2);

% intensity for Interval 1 and Interval 2
intensity = [zeros(91,1); repmat(1,162,1)];
intensity = [intensity, hurstMatrix(:,I1), hurstMatrix(:,I2)];

% Randomization process
sample_length = 253;
traning = randperm(sample_length);
traning = traning(1 : (floor(2*sample_length./3)));
evaluant = [1:sample_length];
evaluant([traning]) = []; % the rest of data(1/3 of data)
x = [];
y = [];
for i = 1 : length(traning)
    sample = intensity(traning(i),:);
    y = [y; sample(1)];
    x = [x; sample(2:3)];
end

% Simulate
%b = glmfit(x,y);
%b = glmfit(x,y,'normal','link','identity'); % need to change the probability function
b = glmfit(x,y,'binomial','link','logit'); % by using generalized linear model regression
%b = glmfit(x,y,'binomial','link','probit');
%b = glmfit(x,y,'binomial','link','comploglog'); % good for specificity
b0 = b(1); b1 = b(2); b2 = b(3); %evaluate b0, b1 and b2
tp = 0; fp = 0; fn = 0; tn = 0;

for i = 1 : length(evaluant)
    sample = intensity(evaluant(i),:);
    H1 = sample(2); H2 = sample(3);
    linpart = b0 + b1*H1 + b2*H2;
    pp=exp(linpart)./(1+exp(linpart));
    if pp >= prediction; % Test CANCER
        if sample(1) == 1 % True Positive
            tp = tp + 1;
        else % False Positive
            fp = fp + 1;
        end
    else
        if sample(1) == 1 % false negative
            fn = fn + 1;
        else % True Negative
            tn = tn + 1;
        end
    end
end

confusionMatrix = [tp,fp;fn,tn];

[se sp pre ppv npv lrp accuracy yi] = sesp(tp, fp, fn, tn);


