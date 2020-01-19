%Date 30 Sep 2019 copy from Journal_2CMV_Fullmeasurement_No_l2NormonW_cora.m
%add para in the function to allow choosing norm for W when running
%date: 13Jan 2019. Copy from TwoC_For_MV_ComplementaryManifold_l2normonW,
%TwoC_For_MV _learning using both coupled matrix and normal NMF, add diverse between two views, l2 norm on H_v,
%update rule similar to diverNMF
%parameters: beta for the diverse between H_v and H_star and eta for the l2_norm
%note: alreadly use l2-norm on H_star

%% add Learning L_complementary with the parameter gama


function [Hv, Hstar_final, nIter_final, W_final, objhistory_final, nIteration] = TwoCMV_function(normonW,normonH, gnd, m, n,  R, H, Lcomple, nClass, nfea, options)
if ~isfield(options,'error')
    options.error = 1e-6;
end
if ~isfield(options, 'maxIter')
    options.maxIter = [];
end

if ~isfield(options,'nRepeat')
    options.nRepeat = 10;
end

if ~isfield(options,'minIter')
    options.minIter = 30;
end

if ~isfield(options,'meanFitRatio')
    options.meanFitRatio = 0.1;
end

if ~isfield(options,'alpha')
    options.alpha = 100;
end

if ~isfield(options,'Optimization')
    options.Optimization = 'Multiplicative';
end

if ~exist('H1','var') %k
    H1 = [];
    H2 = [];
    W = [];
end

differror = options.error;
maxIter = options.maxIter;
nRepeat = options.nRepeat;
minIter = options.minIter - 1;
if ~isempty(maxIter) && maxIter < minIter
    minIter = maxIter;
end
meanFitRatio = options.meanFitRatio;

H = l1_norm(H,m);
for i = 1:m
    Hv{i} = H{i};
end
H_star = zeros(n,nClass);
for i = 1:m
    H_star = H_star + 1/m*H{i};
end
H_star = l1_norm_onematrix(H_star);
Hstar_final = H_star;

u = ones(1, 1);
u = [1];
q = size(u,2);%number of W
sumu = sum(u);
for i = 1:size(u,2)
    u(:,i) = u(:,i)/sumu;
end

%% Normalize the Optimal Manifold
Lcomple1 = cell(m,1);
Lcomple0 = cell(m,1);

Lcomple1 = zeros(size(Lcomple,1),size(Lcomple,2));
Lcomple0 = zeros(size(Lcomple,1),size(Lcomple,2));
Lcomple1 = (abs(Lcomple) + Lcomple)*0.5;
Lcomple0 = (abs(Lcomple) - Lcomple)*0.5;

%% end L

selectInit = 1;
Rd_t = R{1,1};
nSmp = size(Rd_t,1);
mFea = size(Rd_t,2);
% Initialize the data and feature matrices
H = initializeMV2018(R, nClass, m);

% initialise matrix W
W = cell(m,1);
for i=1:m
    W{i}=ones(nfea(i),nClass);
end
tryNo = 0;
nIter = 0;
nIteration = 0;

while tryNo < nRepeat
    tryNo = tryNo+1;
    maxErr = 1;
    while(maxErr > differror)
        nIteration = nIteration + 1;
        
        % ===================== update W ~~ update all W_v========================
        W = updateW(W, H, H_star, R, m, nfea, nClass, options);
        if normonW == 1
            W = l1_norm(W,m);
        else
            if normonW == 2
                for v = 1:m
                    W{v,1} = NormalizeFea(W{v,1});
                end
            end
        end
        % ===================== update H_star ~~ update H_star and update H ~~ update all H_v========================
        H_star = updateH_star(W, H, H_star,R,m, nClass, nSmp, options);
        H = updateH(W, H, H_star,R,m, Lcomple, Lcomple1, Lcomple0,  nClass, nSmp, options);        
        if normonH ==1
            H = l1_norm(H,m);
        end
        nIter = nIter + 1;
        %  When U, V run nIter times
        if nIter > minIter
            if selectInit
                objhistory =  CalculateObj(R, W, H, H_star, Lcomple, m, options);
                maxErr = 0;
            else
                if isempty(maxIter)
                    newobj =  CalculateObj(R, W, H, H_star, Lcomple, m, options);
                    objhistory = [objhistory newobj]; 
                    meanFit = meanFitRatio*meanFit + (1-meanFitRatio)*newobj;
                    maxErr = (meanFit-newobj)/meanFit;
                else
                    if isfield(options,'Converge') && options.Converge
                        newobj =  CalculateObj(R, W, H, H_star, Lcomple, m, options);
                        
                        objhistory = [objhistory newobj]; 
                        meanFit = meanFitRatio*meanFit + (1-meanFitRatio)*newobj;%k
                        maxErr = (meanFit-newobj)/meanFit;% k
                        
                    end
                    if nIter >= maxIter
                        maxErr = 0;
                        if isfield(options,'Converge') && options.Converge
                        else
                            objhistory = 0;
                        end
                        
                    end
                end
            end
        end
    end
    %     When nIter achieves minIter, run the following code segment
    if tryNo == 1
        for i = 1:m
            Hv{i} = H{i};
        end
        Hstar_final = H_star;
        for i=1:m
            W_final{i } = W{i };
        end
        nIter_final = nIter;
        objhistory_final = objhistory;
    else
        if objhistory(end) < objhistory_final(end)
            for i = 1:m
                Hv{i} = H{i};
            end
            Hstar_final = H_star;
            for i=1:m
                W_final{i } = W{i };
            end
            nIter_final = nIter;
            objhistory_final = objhistory;
        end
    end
    
    if selectInit
        if tryNo < nRepeat
            %re-start
            H = initializeMV2018(R, nClass,m);
            
            nIter = 0;
        else
            tryNo = tryNo - 1;
            nIter = minIter+1;
            selectInit = 0;
            for i = 1:m
                H{i} = Hv{i};
            end
            Hstar = Hstar_final;
            for i=1:m
                W{i } = W_final{i };
            end
            objhistory = objhistory_final;
            meanFit = objhistory*10;
        end
    end
end
%==========================================================================
function objhistory_final = CalculateObj(R,W, H,H_star, Lcomple, m, options)
obj_NMF = 0;
for i=1:m
    obj_NMF = obj_NMF + norm(R{i} - H{i}*W{i}','fro');
end

obj_NMF_coupled = 0;
for i=1:m
    obj_NMF_coupled = obj_NMF_coupled + norm(R{i} - H_star*W{i}','fro');
end

obj_diver = 0;

for v = 1:m
    obj_diver = obj_diver + trace(H_star*H{v}');
end
obj_l2norm = 0;
for v = 1:m
    obj_l2norm = obj_l2norm + trace(H{v}'*H{v});
end
%Complementary manifold for H_v
obj_compleManifold = 0;
for v = 1:m
    obj_compleManifold = obj_compleManifold + trace(H{v}'*Lcomple*H{v});
end
obj_manifold = options.alpha*obj_compleManifold;
objhistory_final = obj_NMF + obj_NMF_coupled + options.beta*obj_diver + options.eta*obj_l2norm + options.eta*trace(H_star'*H_star) + obj_manifold; % + beta*obj_norm;


%%
function H = l1_norm(H,m)
for p = 1:m
    for i = 1:size(H{p},1)
        if sum(H{p}(i,:))~= 0
            H{p}(i,:) = H{p}(i,:)/sum(H{p}(i,:));
        else
            for j = 1:size(H{p},2)
                H{p}(i,j) = 1/(size(H{p},2));
            end
        end
    end
end
%%
function H = l1_norm_onematrix(H)
for i = 1:size(H,1)
    if sum(H(i,:))~= 0
        H(i,:) = H(i,:)/sum(H(i,:));
    else
        for j = 1:size(H,2)
            H(i,j) = 1/(size(H,2));
        end
    end
end

function W = updateW(W,H,H_star, R, m, nfea, nClass, options) %update all W_v
ceta = options.ceta;
para_nmf = options.para_nmf;
for v=1:m
    tempup = zeros(nfea(v),nClass);
    tempun = zeros(nfea(v),nClass);
    
    tempup = tempup + R{v}'*H{v} + R{v}'*H_star;
    tempun = tempun + W{v}*H{v}'*H{v} + W{v}*H_star'*H_star;
    W{v} = W{v}.*power((tempup./tempun),(0.5));
    % 		W{i} = W{i}.*(VV1./max(VV2,1e-10));
end

%% Update H_star
function H_star = updateH_star(W, H, H_star, R, m, nClass, nSmp, options)
%alpha = options.alpha;
beta = options.beta;
eta = options.eta;

tempup_1 = zeros(nSmp, nClass);
tempun_1 =zeros(nClass, nClass);
for v = 1:m
    tempup_1 = tempup_1 + R{v}*W{v};
    tempun_1 = tempun_1 + W{v}'*W{v};
end
tempup = tempup_1;

sumH_v = zeros(nSmp, nClass);
for v = 1:m
    sumH_v = sumH_v + beta/m*H{v};
end
tempun = H_star*tempun_1 + sumH_v + eta*H_star; % + alpha*Lcompati1*H_star; %setting alpha = 0 to ignore Lcompati

for j = 1:size(H_star,2)
    for i = 1:size(H_star,1)
        if tempun(i,j)~=0
            H_star(i,j) = H_star(i,j)*(tempup(i,j)/tempun(i,j))^(0.5);
        else
            H_star(i,j) = 0;
        end
    end
end

%% update all H_v
function H = updateH(W, H, H_star, R, m, Lcomple, Lcomple1, Lcomple0, nClass, nSmp, options)
alpha = options.alpha;
beta = options.beta;
eta = options.eta;
para_nmf = options.para_nmf;
for v = 1:m
    tempup = zeros(nSmp,nClass);
    tempun = zeros(nSmp,nClass);
    tempup = tempup + R{v}*W{v};
    tempup = tempup + alpha*Lcomple0*H{v};
    sumH_t = zeros(nSmp, nClass);
    tempun = tempun + H{v}*W{v}'*W{v}+ 1/m*beta*H_star + eta*H{v}+ alpha*Lcomple1*H{v};
    
    for j = 1:size(H{v},2)
        for i = 1:size(H{v},1)
            if tempun(i,j)~=0
                H{v}(i,j) = H{v}(i,j)*(tempup(i,j)/tempun(i,j))^(0.5);
            else
                H{v}(i,j) = 0;
            end
        end
    end
end




