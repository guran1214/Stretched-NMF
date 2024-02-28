function pc = pearson(data)
% calculate Pearson correlation
n = size(data,1);
x = data(:,1); 
y = data(:,2);
pc = (n*x'*y-sum(x)*sum(y))*(n*x'*x-sum(x)^2)^-0.5*(n*y'*y-sum(y)^2)^-0.5;
end