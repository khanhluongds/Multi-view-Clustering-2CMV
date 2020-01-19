clear;
addpath(genpath(pwd));
%% Load dataset and read dataset into a cell

datasetname = 'cora_formated';
load(datasetname);
fprintf('Dataset name: %s\n',datasetname);

%% initialise R is the data matrix 
nClass = length(unique(gnd)); %document cluster number
m = nviews; %number of views
R = cell(m,1);
n = size(fea{1,1},1);
nfea = zeros(m,1);
for i = 1:m
    nfea(i) = size(fea{i,1},2);
end
for i = 1:m
    temp = zeros(n,nfea(i));
    R{i,1} = temp;
end
%% parameter setting norm or no norm, and for calculation of affinity matrix and laplacian matrix
options = [];
options.alpha = 0;
options.WeightMode = 'Binary';
mode = options.WeightMode;
options.maxIter = 100;
rand('twister',5489);

normonFea = 2;
normonS = 0;
normonG = 0;%default value 0;

if normonFea == 2
    for i = 1:m
        R{i,1} = NormalizeFea(fea{i,1});
    end
end

%Clustering in the original space
rand('twister',5489); 

%First view matrix or concatenate all views to be clustered
firstView = R{1,1};
concat = [R{1,1} R{2,1}];
disp(['Clustering in the original space, concatenating of all views.']);
[best.CA best.F best.P best.R best.nmi best.AR] = performance_kmeans(concat, nClass, gnd);
disp(['    NMI and std:       ',num2str(best.nmi(1)), ' , ', num2str(best.nmi(2))]);
disp(['    Accuracy and std:  ',num2str(best.CA(1)), ' , ', num2str(best.CA(2))]);
disp(['    F-score and std:   ',num2str(best.F(1)), ' , ', num2str(best.F(2))]);
% initialise, assign intra-type relationship and calculate optimal manifold
W = cell(m,1);
L = cell(m,1);
tempW = cell(m,1);
alpha = 1;
options.alpha = alpha;
options.normW = 1;
H = cell(m,1);
for i = 1:m
    H{i,1} = zeros(n, nClass);
end

G = initializeMV2018(R, nClass, m);
rand('twister',5489);
ncols = 10;
alphas =[1];
betas = [1 0];
para_nmfs = [1];
ks = [25];
%%
options.eta = 1;
NMIs_intra_new = para_matrix(alphas, betas, para_nmfs, ks, 15);
for iNMI = 1:size(NMIs_intra_new,1)
    options.alpha = NMIs_intra_new(iNMI,1);
    options.beta = NMIs_intra_new(iNMI,2);
    options.para_nmf = NMIs_intra_new(iNMI,3);
    options.k = NMIs_intra_new(iNMI,4);
    if (options.alpha ~= 0) & (options.beta~=0)
        disp(['2CMV - Clustering result']);
    end
    if (options.alpha ~= 0) & (options.beta==0)
        disp(['2CMV-OM - Clustering result']);
    end
    if (options.para_nmf == 0.5)
        options.ceta = 0.5; % - options.para_nmf;
    else
        options.ceta = 1;
    end
    % initialise and assign intra-type relationship
    W = cell(m,1);
    tempW = cell(m,1);
    alpha = 1;
    options.normW = 1;
    
    for v = 1:m
        W{v,1} = constructW(R{v},options);
        L{v,1} = constructL(W{v}, alpha, options);
    end
    %% Calculate the Optimal affinity matrix Wopt and the Optimal Manifold Lopt
    for i = 1:m
        W{i,1} = full(W{i,1});
    end
    for i = 1:m
        Wopt(:,:,i) = W{i,1};
    end
    Wopt = min(Wopt,[],3);
    Lopt = constructL(Wopt, alpha, options);
    
    %%
    tic
    [Hv, H_star, nIter_final, W_final, objhistory_final,nIteration] = TwoCMV_function(normonS, normonG, gnd, m, n, R, H, Lopt, nClass, nfea, options);
    time2CMV = toc;
    rand('twister',5489);
    H_com = zeros(n,nClass);
    for t=1:m
        H_com = H_com + 1/m*Hv{t};
    end
    %%
    para_1 = 0.8;%para_coupled(itest);
    para_2 = 1 - para_1;
    H_final = [para_1*H_star para_2*H_com];
    label = litekmeans(H_final,nClass,'Replicates',20);%%%%
    %%
    if sum(any(isnan(H_final), 2))~=n
        [best.CA best.F best.P best.R best.nmi best.AR] = performance_kmeans(H_final, nClass, gnd);
        disp(['    NMI and std:       ',num2str(best.nmi(1)), ' , ', num2str(best.nmi(2))]);
        disp(['    Accuracy and std:  ',num2str(best.CA(1)), ' , ', num2str(best.CA(2))]);
        disp(['    F-score and std:   ',num2str(best.F(1)), ' , ', num2str(best.F(2))]);
    end
    
end
