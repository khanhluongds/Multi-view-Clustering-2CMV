% a = [1 2 3 4];
% b = [0.5 1 2];
% c = [10 100];
% d = [1 5];

function M = para_matrix(a, b, c, d, ncols)

na = length(a); 
nb = length(b); 
nc = length(c); 
nd = length(d); 
nrows = na*nb*nc*nd;
M = zeros(nrows,ncols);

% a
nrepa = nb*nc*nd;
for k = 1:1
    for i = 1:na
        M((na*nrepa)*(k-1)+ nrepa*(i-1)+1:(na*nrepa)*(k-1)+ nrepa*(i-1)+nrepa,1) = repmat(a(i),nrepa,1);
%         them gia tri k vao cong thuc nay de thay cong thuc tong quat thoi
%         chu thuc ra o day k = 1 nen (na*nrepa)*(k-1) luon = 0
    end
end
% b
% nc*nd is the replicating times for each element of b[]
% truoc bien k la so phan tu mang dang xet * so lan lap moi element
% truoc bien i la so lan lap moi element
nrepb = nc*nd;
for k = 1:na
    for i = 1:nb
%         M(12*(k-1)+nc*nd*(i-1)+1:12*(k-1)+nc*nd*(i-1)+4,2)=repmat(b(i),nc*nd,1);
        M((nb*nrepb)*(k-1)+nrepb*(i-1)+1:(nb*nrepb)*(k-1)+nrepb*(i-1)+nrepb,2)=repmat(b(i),nrepb,1);
    end
end


% (so phan tu cua mang dang xet * so lan rep moi phan tu)(k-1) vi do gia
% tri buoc nhay cua k
% Truoc bien i la gia tri buoc nhay cua i, cung chinh la so lan lap lai cua
% moi gia tri trong mang dang xet
% nd so lan rep
nrepc = nd;
for k = 1:na*nb
    for i = 1:nc
        M((nc*nrepc)*(k-1)+nrepc*(i-1)+1:(nc*nrepc)*(k-1)+nrepc*(i-1)+ nrepc, 3)=repmat(c(i),nrepc,1);
%         M(4*(k-1)+nd*(i-1)+1:4*(k-1)+nd*(i-1)+ nd, 2)=repmat(b(i),nd,1);
    end
end

nrepd = 1;
for k = 1:na*nb*nc
    for i = 1:nd
        M((nd*nrepd)*(k-1)+nrepd*(i-1)+1:(nd)*(k-1)+nrepd*(i-1)+ nrepd, 4)=repmat(d(i),nrepd,1);
    end
end