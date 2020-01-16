function [G] = initializeMV2018(R, nClass, m)
label = cell(m,1);
rand('twister',5489);
% label{1} = litekmeans(R{1}, nClass); 
for i = 1:m
    rand('twister',5489);
    label{i} = litekmeans(R{i,1}, nClass); 
end
for h = 1:m
    for i = 1:nClass  % nCla == mCla
        G{h}(label{h} ==i, i) = 1;
    end
end
for i = 1:m
    G{i} = G{i}+0.2;
end    
