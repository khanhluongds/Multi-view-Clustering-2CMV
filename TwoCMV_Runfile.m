%%
clear;
addpath(genpath(pwd));
% path_ref = '/home/luong5/EXPERIMENT/Datasets'
% dataset_ref = '/home/luong5/EXPERIMENT/ref_funcs'
% addpath(genpath(path_ref));
% addpath(genpath(dataset_ref));

%% Load dataset and read dataset into a cell

datasetname = '3-sources_formated';
%%% Read dataset and store into cell R
load(datasetname);
m = nviews; %number of views 

% initialise R is a matrix of mxm matrix storing inter relationships
% between object types
R = cell(m,1);
n = size(fea{1,1},1);
nfea = zeros(m,1); %nfea{1} stores the number of objects in object type 1
%nfea(1) = size(fea{1},1);

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
normonS = 2;
normonG = 0;%default value 0;

if normonFea == 2
    for i = 1:m
        R{i,1} = NormalizeFea(fea{i,1});
    end
end


 
%% 
timerun = datestr(datetime('now'));
filename_mat = strcat(datasetname, '_2CMV_Jour_2019','_normonFea', num2str(normonFea),'_normonS',num2str(normonS),'_normonG',num2str(normonG),'_', mode,'_', timerun);
filename_mat(regexp(filename_mat,'[.]'))=[];
filename_mat(regexp(filename_mat,'[ ]'))=[];
filename_mat(regexp(filename_mat,'[-]'))=[];
filename_mat(regexp(filename_mat,'[:]'))=[];
 
fprintf('Dataset name: %s\n',datasetname);
 
nClass = length(unique(gnd)); %document cluster number
 
 
%Clustering in the original space
rand('twister',5489);

%First view matrix or concatenate all views to be clustered
firstView = R{1,1};
concat = [R{1,1} R{2,1}];
label = litekmeans(concat,nClass,'Replicates',20);
NMI_Kmeans = MutualInfo(gnd,label);
disp(['Clustering in the original space. MIhat: ',num2str(NMI_Kmeans)]);

 d = [gnd, label];
%     fscore = FScr(d);
gnd1 = gnd;
labelnew = bestMap(gnd1, label);
AC_Kmeans = length(find(gnd == labelnew))/length(gnd);

% initialise and assign intra-type relationship
W = cell(m,1);
L = cell(m,1);
tempW = cell(m,1);
alpha = 1;
options.alpha = alpha;
options.normW = 1;
for i = 1:m
    W{i} = constructW(R{i},options);
    L{i} = constructL(W{i}, alpha, options);
end
 
G = cell(m,1);
for i = 1:m
    G{i,1} = zeros(n, nClass);
end
% G = initializeMulti_view_anyviews(R, nClass);
G = initializeMV2018(R, nClass, m);
rand('twister',5489);
ncols = 10; 
alphas =[0.1 1];
betas = [1]; %10; %0.663 
para_nmfs = [1]; %10;
ks = [5]; 
%%
options.eta = 1;

% NMIs_intra = para_matrix(alphas, betas, para_nmfs, ks, ncols+1);
NMIs_intra_new = para_matrix(alphas, betas, para_nmfs, ks, 15);
NMIs_Gcoupled_new = para_matrix(alphas, betas, para_nmfs, ks, 15);
NMIs_G_comple_new = para_matrix(alphas, betas, para_nmfs, ks, 15);


% NMIs_G_comple = NMIs_intra; 
% NMIs_Gcoupled = NMIs_intra;
para_coupled = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
 
    NMIparas = zeros(size(NMIs_intra_new,1),length(para_coupled));
    ACCparas = zeros(size(NMIs_intra_new,1),length(para_coupled));
    ARparas = zeros(size(NMIs_intra_new,1),length(para_coupled));
    Fscoreparas = zeros(size(NMIs_intra_new,1),length(para_coupled));
    Precisionparas = zeros(size(NMIs_intra_new,1),length(para_coupled));
    Recallparas = zeros(size(NMIs_intra_new,1),length(para_coupled));
    
    NMIparas_std = zeros(size(NMIs_intra_new,1),length(para_coupled));
    ACCparas_std = zeros(size(NMIs_intra_new,1),length(para_coupled));
    Fscoreparas_std = zeros(size(NMIs_intra_new,1),length(para_coupled));
    ARparas_std = zeros(size(NMIs_intra_new,1),length(para_coupled));
    Recallparas_std = zeros(size(NMIs_intra_new,1),length(para_coupled));
    Precisionparas_std = zeros(size(NMIs_intra_new,1),length(para_coupled));
    
    ACC_grid = zeros(size(NMIs_intra_new,1),length(para_coupled));
    NMI_grid = zeros(size(NMIs_intra_new,1),length(para_coupled));
    Fscore_grid = zeros(size(NMIs_intra_new,1),length(para_coupled));
    
for iNMI = 1:size(NMIs_intra_new,1)
    options.alpha = NMIs_intra_new(iNMI,1);
    options.beta = NMIs_intra_new(iNMI,2);
    options.para_nmf = NMIs_intra_new(iNMI,3);
    options.k = NMIs_intra_new(iNMI,4);
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
    [G_final, Gstar_final, nIter_final, S_final, objhistory_final,nIteration] = TwoCMV_function(normonS, normonG, gnd, m, n, R, G, Lopt, nClass, nfea, options);
    time2CMV = toc;
    rand('twister',5489);
%    %%
    G_comple = zeros(n,nClass);
    
    for t=1:m
        G_comple = G_comple + 1/m*G_final{t};
    end
%     label = litekmeans(Gstar,nClass,'Replicates',20); %%


%%  old calculation for NMIs_intra

    label = litekmeans([Gstar_final G_comple],nClass,'Replicates',20);
    NMIs_intra(iNMI,5) = MutualInfo(gnd,label);
    gnd1 = gnd;
    labelnew = bestMap(gnd1, label);
    AC = length(find(gnd == labelnew))/length(gnd);
    NMIs_intra(iNMI,6) = AC;
    d = [gnd, label];
    fscore = FScr(d);
    NMIs_intra(iNMI,7) = fscore;
    NMIs_intra(iNMI,10) = time2CMV;   
    G2C = [Gstar_final G_comple];
%% new measurement when two parameters for G_2C are turned off.
if sum(any(isnan(G2C), 2))~=n
    [CA F P Recall nmi AR] = performance_kmeans(G2C, nClass, gnd);

    NMIs_intra_new(iNMI, [5,6]) = nmi;
    NMIs_intra_new(iNMI, [7,8]) = CA;
    NMIs_intra_new(iNMI, [9,10]) = F;
    NMIs_intra_new(iNMI, [11,12]) = P;
    NMIs_intra_new(iNMI, 14) = time2CMV;
end

%%    
for itest = 1:length(para_coupled)
        para_1 = para_coupled(itest);
        para_2 = 1 - para_1;
        G_2C = [para_1*Gstar_final para_2*G_comple];
        
        %%
        if sum(any(isnan(G_2C), 2))~=n
            [CA F P Recall nmi AR] = performance_kmeans(G_2C, nClass, gnd);

    %         label = litekmeans(Gstar,nClass,'Replicates',20);

            NMIparas(iNMI,itest) = nmi(1);
            NMIparas_std(iNMI,itest) = nmi(2);

            ACCparas(iNMI,itest) = CA(1);
            ACCparas_std(iNMI,itest) = CA(2);

            Fscoreparas(iNMI,itest) = F(1);
            Fscoreparas_std(iNMI,itest) = F(2);

            Recallparas(iNMI,itest) = Recall(1);
            Recallparas_std(iNMI,itest) = Recall(2);

            Precisionparas(iNMI,itest) = P(1);
            Precisionparas_std(iNMI,itest) = P(2);

            ARparas(iNMI,itest) = AR(1);
            ARparas_std(iNMI,itest) = AR(2);
        end
        
        
        %% old measuarement
        label = litekmeans(G_2C,nClass,'Replicates',20);
        
        NMI_grid(iNMI,itest) = MutualInfo(gnd,label);
        gnd1 = gnd;
        d = [gnd, label];
        fscore = FScr(d);
        Fscore_grid(iNMI, itest) = fscore;
        labelnew = bestMap(gnd1, label);
        AC = length(find(gnd == labelnew))/length(gnd);
        ACC_grid(iNMI,itest) = AC;
        
  
end


    
%% New measurement result
    AResult_NMIparas = [NMIparas NMIs_intra_new NMIparas_std];
    AResult_ACCparas = [ACCparas NMIs_intra_new ACCparas_std];
    AResult_ARparas = [ARparas NMIs_intra_new ARparas_std];
    AResult_Recallparas = [Recallparas NMIs_intra_new Recallparas_std];
    AResult_Precisionparas = [Precisionparas NMIs_intra_new Precisionparas_std];
    AResult_Fscoreparas = [Fscoreparas NMIs_intra_new Fscoreparas_std];
    

%%
  
      
    %use couple matrix
    label = litekmeans(Gstar_final,nClass,'Replicates',20);
    NMI_temp = MutualInfo(gnd,label);
    NMIs_Gcoupled(iNMI, 5) = NMI_temp;
    gnd1 = gnd;
    labelnew = bestMap(gnd1, label);
    AC = length(find(gnd == labelnew))/length(gnd);
    NMIs_Gcoupled(iNMI,8) = AC;
    %new measurement
    if sum(any(isnan(Gstar_final), 2))~=n
        [CA F P Recall nmi AR] = performance_kmeans(Gstar_final, nClass, gnd);   
        NMIs_Gcoupled_new(iNMI, [5,6]) = nmi;
        NMIs_Gcoupled_new(iNMI, [7,8]) = CA;
        NMIs_Gcoupled_new(iNMI, [9,10]) = F;
        NMIs_Gcoupled_new(iNMI, [11,12]) = P;
        NMIs_Gcoupled_new(iNMI, 14) = time2CMV;
    end
    %use G_Comple
    label = litekmeans(G_comple,nClass,'Replicates',20);
    NMI_temp = MutualInfo(gnd,label);
    NMIs_G_comple(iNMI, 5) = NMI_temp;
    gnd1 = gnd;
    labelnew = bestMap(gnd1, label);
    AC = length(find(gnd == labelnew))/length(gnd);
    NMIs_G_comple(iNMI,8) = AC;
    
    %new measurement
    if sum(any(isnan(G_comple), 2))~=n
        [CA F P Recall nmi AR] = performance_kmeans(G_comple, nClass, gnd);   
        NMIs_G_comple_new(iNMI, [5,6]) = nmi;
        NMIs_G_comple_new(iNMI, [7,8]) = CA;
        NMIs_G_comple_new(iNMI, [9,10]) = F;
        NMIs_G_comple_new(iNMI, [11,12]) = P;
        NMIs_G_comple_new(iNMI, 14) = time2CMV;
    end
        A_NMI_result = [NMI_grid NMIs_intra]; 
        A_ACC_result = [ACC_grid NMIs_intra]; 
        A_Fscore_result = [Fscore_grid NMIs_intra];
    
    note_NMIs_intra_new = 'NMIs_intra_new has been ordered as NMI and NMI_std ACC and ACC_std F-score and F-score_std Precision and Precision_std -- time'
    save(filename_mat,'NMIs_G_comple_new','NMIs_Gcoupled_new','S_final', 'G_final','NMIs_intra_new', 'datasetname','fea','gnd', 'note_NMIs_intra_new','NMIs_intra', 'NMI_Kmeans', 'NMIs_G_comple','NMIs_Gcoupled', 'A_NMI_result', 'A_ACC_result', 'A_Fscore_result', 'AResult_NMIparas', 'AResult_ACCparas', 'AResult_ARparas', 'AResult_Recallparas', 'AResult_Precisionparas', 'AResult_Fscoreparas' );
end
 
