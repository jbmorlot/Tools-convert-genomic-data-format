function [matrix,threshold,Outliersmat] =  Outlier_Removal_density(matrix)
%---------------------------
%PARAMETERS: Can be changed
num_bin = 100;      %Number of intervals to classify values
size_bin = 100;
threshold_density = 10; %Density threshold
%---------------------------

[M,N] = size(matrix);

Outliersmat = zeros(M,N);
threshold = zeros(M,1);

for i=1:M
tic;
h1 = hist(matrix(i,:),num_bin);
v = min(matrix(i,:)):size_bin:max(matrix(i,:));
% v = linspace(min(matrix(i,:)),max(matrix(i,:)),num_bin);
threshold(i) = v(find(h1<threshold_density,1,'first'));
Outliersmat(i,matrix(i,:)>threshold(i)) = 1;
matrix(i,matrix(i,:)>threshold(i))=threshold(i);
disp(['-------- Line ',num2str(i),' ------------'])
disp(['threshold = ',num2str(floor(threshold(i)))])
disp(['Number of Outliers Removed / Size matrix = ',num2str(sum(Outliersmat(i,:))),' / ',num2str(N)])
toc,
disp('----------------------------')
end
Outliersmat = sparse(Outliersmat);


end