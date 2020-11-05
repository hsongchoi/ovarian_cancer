function [confusionMatrix accuracy se sp]=  Ovarian2(prediction)

load('hurstMatrix.mat')

% divde control and cases
% Get means, medians and standard deviation of control and case

control = hurstMatrix(1:91,:); cases = hurstMatrix(92:253,:);
control_data= [mean(control); median(control); std(control)];
cases_data = [mean(cases); median(cases); std(cases)];
diff = abs(control_data(1,:) - cases_data(1,:));
[B,I] = sort(diff,'descend');
interval = I(1:5);

% intensity for Interval 1 and Interval 2
intensity = [zeros(91,1); repmat(1,162,1)]; %1111000000
intensity = [intensity, hurstMatrix(:,interval(1)), hurstMatrix(:,interval(2))];

% Randomization process
sample_length = 253;
traning = randperm(sample_length);
traning = traning(1 : (floor(2*sample_length./3)));
evaluant = [1:sample_length];
evaluant([traning]) = []; % the rest of data(1/3 of data)
x = [];
y = [];
for i = 1 : length(traning)
    sample = intensity(traning(i),:)
    y = [y; sample(1)];
    x = [x; sample(2:3)];
end

% Simulate
b = glmfit(x,y); % by using generalized linear model regression
b0 = b(1); b1 = b(2); b2 = b(3) %evaluate b0 and b1%evaluate b0 and b1
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






%% 
[f,alfa]=mfs(data,q,j1,j2,wf, L)
a = MakeFBMNew(1024,1/2);
plot(a)
harr = [sqrt(2)/2, sqrt(2)/2];
[a,b]=mfs(a,-1:0.1:6,3,8,harr, 2);
plot(b,a,'-')
xx = 0.2:0.1:1.8;
hold on
plot(xx,-0.2*ones(size(xx)),'r_')
plot(xx,-0*ones(size(xx)),'r_')
