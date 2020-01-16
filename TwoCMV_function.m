%Date 30 Sep 2019 copy from Journal_2CMV_Fullmeasurement_No_l2NormonS_cora.m
%add para in the function to allow choosing norm for S when running
%date: 13Jan 2019. Copy from TwoC_For_MV_ComplementaryManifold_l2normonS,
%TwoC_For_MV _learning using both coupled matrix and normal NMF, add diverse between two views, l2 norm on H_v, 
%update rule similar to diverNMF
%parameters: beta for the diverse between G_v and G_star and eta for the l2_norm
%note: alreadly use l2-norm on G_star

%% add Learning L_complementary with the parameter gama

                                                                                                                                  
function [G_final, Gstar_final, nIter_final, S_final, objhistory_final, nIteration] = TwoCMV_nocetanmf(normonS,normonG, gnd, m, n,  R, G, Lcomple, nClass, nfea, options)
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

if ~exist('G1','var') %k
    G1 = [];
    G2 = [];
    S = [];
end

differror = options.error;
maxIter = options.maxIter;
nRepeat = options.nRepeat;
minIter = options.minIter - 1;
if ~isempty(maxIter) && maxIter < minIter
    minIter = maxIter;
end
meanFitRatio = options.meanFitRatio;

G = l1_norm(G,m);
for i = 1:m
    G_final{i} = G{i};
end
G_star = zeros(n,nClass);
for i = 1:m
    G_star = G_star + 1/m*G{i};
end;
G_star = l1_norm_onematrix(G_star);
Gstar_final = G_star;

Norm = 2;
NormV = 0;

%construct L based on u

u = ones(1, 1);
u = [1];
q = size(u,2);%number of W
sumu = sum(u);
for i = 1:size(u,2)
    u(:,i) = u(:,i)/sumu;
end

% u_final = u;
% Lnorm = Lcompati;
% for i = 1:m
%     Lnorm{i} = zeros(n,n);
%     for j = 1:size(u_final,2)
%         Lnorm{i} = L{i};
%     end
% end

% for i = 1:m
%     Lcomple = Lcomple + 1/m*Lnorm{i};
% end

%% calculate and use new Lcompatible
Lcomple1 = cell(m,1);
Lcomple0 = cell(m,1);

Lcomple1 = zeros(size(Lcomple,1),size(Lcomple,2));
Lcomple0 = zeros(size(Lcomple,1),size(Lcomple,2));
Lcomple1 = (abs(Lcomple) + Lcomple)*0.5;
Lcomple0 = (abs(Lcomple) - Lcomple)*0.5;

%% end chuan bi L

selectInit = 1;
Rd_t = R{1,1};
nSmp = size(Rd_t,1);
mFea = size(Rd_t,2);
% Initialize the data and feature matrices
G = initializeMV2018(R, nClass, m);

% initialise matrix S
S = cell(m,1);
for i=1:m
    S{i}=ones(nfea(i),nClass);
end
tryNo = 0;
nIter = 0;
nIteration = 0;
% options.ceta = 1;
% options.eta = 1;
% options.beta = 1;
ceta = options.ceta
eta = options.eta
para_nmf = options.para_nmf

while tryNo < nRepeat
    tryNo = tryNo+1;
    maxErr = 1;
    while(maxErr > differror)
        nIteration = nIteration + 1;
        
        % ===================== update S ~~ update all W_v========================
        S = updateS(S, G, G_star, R, m, nfea, nClass, options);
        
        if normonS == 1
            S = l1_norm(S,m);
        else 
            if normonS == 2
                for v = 1:m
                    S{v,1} = NormalizeFea(S{v,1});
                end
            end
        end
        

%         S = l1_norm(S,m);
        % ===================== update G ~~ update all H_v========================
        G_star = updateG_star(S, G, G_star,R,m, nClass, nSmp, options); 

%         G_star = NormalizeFea(G_star);
%         G_star = l1_norm_onematrix(G_star);
        G = updateG(S, G, G_star,R,m, Lcomple, Lcomple1, Lcomple0,  nClass, nSmp, options);
        % ===================== update G_star ~~ update H_star ========================
        
        if normonG ==1
            G = l1_norm(G,m);
        end
        
        nIter = nIter + 1;
        
        
        %  When U, V run nIter times
        if nIter > minIter
            if selectInit
                objhistory =  CalculateObj(R, S, G, G_star, Lcomple, m, options);
                maxErr = 0;
            else
                if isempty(maxIter)
                    newobj =  CalculateObj(R, S, G, G_star, Lcomple, m, options);
                    objhistory = [objhistory newobj]; %#ok<AGROW>
                    meanFit = meanFitRatio*meanFit + (1-meanFitRatio)*newobj;
                    maxErr = (meanFit-newobj)/meanFit;
                else
                    if isfield(options,'Converge') && options.Converge
                        newobj =  CalculateObj(R, S, G, G_star, Lcomple, m, options);
                        
                        objhistory = [objhistory newobj]; %#ok<AGROW>
                        meanFit = meanFitRatio*meanFit + (1-meanFitRatio)*newobj;%k
                        maxErr = (meanFit-newobj)/meanFit;% k
                        
                    end
                    %                     maxErr = 1;
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
            G_final{i} = G{i};
        end
        Gstar_final = G_star;
        for i=1:m
            S_final{i } = S{i };
        end
        nIter_final = nIter;
        objhistory_final = objhistory;
    else
        if objhistory(end) < objhistory_final(end)
            for i = 1:m
                G_final{i} = G{i};
            end
            Gstar_final = G_star;
            for i=1:m
                S_final{i } = S{i };
            end
            nIter_final = nIter;
            objhistory_final = objhistory;
        end
    end
    
    if selectInit
        if tryNo < nRepeat
            %re-start
            G = initializeMV2018(R, nClass,m);
            
            nIter = 0;
        else
            tryNo = tryNo - 1;
            nIter = minIter+1;
            selectInit = 0;
            for i = 1:m
                G{i} = G_final{i};
            end
            Gstar = Gstar_final;
            for i=1:m
                S{i } = S_final{i };
            end
            objhistory = objhistory_final;
            meanFit = objhistory*10;
        end
    end
end
    %==========================================================================
function objhistory_final = CalculateObj(R,S, G,G_star, Lcomple, m, options)
    obj_NMF = 0;
    for i=1:m
        obj_NMF = obj_NMF + norm(R{i} - G{i}*S{i}','fro');
    end

    obj_NMF_coupled = 0;
    for i=1:m
        obj_NMF_coupled = obj_NMF_coupled + norm(R{i} - G_star*S{i}','fro');
    end

    obj_diver = 0;
    
    for v = 1:m
        obj_diver = obj_diver + trace(G_star*G{v}');
    end
    obj_l2norm = 0; 
    for v = 1:m
        obj_l2norm = obj_l2norm + trace(G{v}'*G{v});
    end
    %Complementary manifold for G_v
    obj_compleManifold = 0;
    for v = 1:m
        obj_compleManifold = obj_compleManifold + trace(G{v}'*Lcomple*G{v});
    end
    obj_manifold = options.alpha*obj_compleManifold; 
    objhistory_final = obj_NMF + obj_NMF_coupled + options.beta*obj_diver + options.eta*obj_l2norm + options.eta*trace(G_star'*G_star) + obj_manifold; % + beta*obj_norm;


%%    
function G = l1_norm(G,m)
    for p = 1:m
        for i = 1:size(G{p},1)
            if sum(G{p}(i,:))~= 0
                G{p}(i,:) = G{p}(i,:)/sum(G{p}(i,:));
            else
                for j = 1:size(G{p},2)
                    G{p}(i,j) = 1/(size(G{p},2));
                end
            end
        end
    end
    %%
    function G = l1_norm_onematrix(G)
         for i = 1:size(G,1)
            if sum(G(i,:))~= 0
                G(i,:) = G(i,:)/sum(G(i,:));
            else
                for j = 1:size(G,2)
                    G(i,j) = 1/(size(G,2));
                end
            end
        end
 
    
    % function [G,S] = NormalizeGS(G, S)
    % n = size(G,1);
    % c = size(S,2);
    %
    % norms = sqrt(sum(G.^2,1));

%%Update all W_v    
% function S = updateS(S,G,G_star, R, m, nfea, nClass) %update all W_v
% for v=1:m
%     tempup = zeros(nfea(v),nClass);
%     tempun = zeros(nfea(v),nClass);
%     
%     tempup = tempup + R{v}'*G{v} + R{v}'*G_star;
%     tempun = tempun + S{v}*G{v}'*G{v} + S{v}*G_star'*G_star;
%     S{v} = S{v}.*power((tempup./tempun),(0.5));
%     % 		S{i} = S{i}.*(VV1./max(VV2,1e-10));
% end
function S = updateS(S,G,G_star, R, m, nfea, nClass, options) %update all W_v
ceta = options.ceta;
para_nmf = options.para_nmf;
for v=1:m
    tempup = zeros(nfea(v),nClass);
    tempun = zeros(nfea(v),nClass);
    
    tempup = tempup + R{v}'*G{v} + R{v}'*G_star;
    tempun = tempun + S{v}*G{v}'*G{v} + S{v}*G_star'*G_star;
    S{v} = S{v}.*power((tempup./tempun),(0.5));
    % 		S{i} = S{i}.*(VV1./max(VV2,1e-10));
end

%% Update G_star
function G_star = updateG_star(S, G, G_star, R, m, nClass, nSmp, options)
%alpha = options.alpha;
beta = options.beta;
eta = options.eta;

tempup_1 = zeros(nSmp, nClass);
tempun_1 =zeros(nClass, nClass);
for v = 1:m
    tempup_1 = tempup_1 + R{v}*S{v};
%     tempup_2 = tempup_2 + G{v}*D';
    tempun_1 = tempun_1 + S{v}'*S{v};
%     tempun_2 = tempun_2 + G{v}*D';
end
tempup = tempup_1;  

sumG_v = zeros(nSmp, nClass);
for v = 1:m
   sumG_v = sumG_v + beta/m*G{v};     
end
tempun = G_star*tempun_1 + sumG_v + eta*G_star; % + alpha*Lcompati1*G_star; %setting alpha = 0 to ignore Lcompati
 
for j = 1:size(G_star,2)
    for i = 1:size(G_star,1)
        if tempun(i,j)~=0
            G_star(i,j) = G_star(i,j)*(tempup(i,j)/tempun(i,j))^(0.5);
        else
            G_star(i,j) = 0;
        end
    end
end

%% update all G_v
function G = updateG(S, G, G_star, R, m, Lcomple, Lcomple1, Lcomple0, nClass, nSmp, options) 
alpha = options.alpha;
beta = options.beta;
eta = options.eta;
para_nmf = options.para_nmf;

for v = 1:m   
    tempup = zeros(nSmp,nClass);
    tempun = zeros(nSmp,nClass);
    tempup = tempup + R{v}*S{v}; 
    tempup = tempup + alpha*Lcomple0*G{v};
    sumG_t = zeros(nSmp, nClass);

    tempun = tempun + G{v}*S{v}'*S{v}+ 1/m*beta*G_star + eta*G{v}+ alpha*Lcomple1*G{v};

    for j = 1:size(G{v},2)
        for i = 1:size(G{v},1)
            if tempun(i,j)~=0
                G{v}(i,j) = G{v}(i,j)*(tempup(i,j)/tempun(i,j))^(0.5);
            else
                G{v}(i,j) = 0;
            end
        end
    end
end

                    
                    
                    
