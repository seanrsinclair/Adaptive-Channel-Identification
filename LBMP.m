function [x_hat,antiomp,lmp]=LBMP(kk,sparsity,groups,A,y)
[m,n]=size(A);
% -- initialization
idx = [];
idx_g=[];
res = y; % initialize residual

inds_A=[n-(n/groups/2):n 1:n/groups/2-1 n/groups/2:n-n/groups/2-1];
groupind=reshape(inds_A,n/groups,groups);
mi={};
% -- main loop
for k=1:sparsity

    % -- calculate correlation between the columns and residual
    for g=1:groups
        Ai=A(:,groupind(:,g));
        MF(g) = norm((Ai'*Ai)^(-1/2)*Ai'*res,2);
    end
    
    if k==1 %simple method
        [~,antiomp] = sort(MF,'ascend');
    end
    %store coefficients
    mi{k}=MF;
    %remove previously found channels
    MF(idx_g) = 0;
    %find one active channel
    [~,idxg] = max(MF);
    %store found channel
    idx=[idx groupind(:,idxg)'];
    idx_g=[idx_g idxg];


    A_omega = A(:,idx);
    x_ls = pinv(A_omega)*y;
    
    %claculate residual(makes it orthogonal)
    res = y - A_omega*x_ls;   

end
%calculate recovered signal
x_hat = zeros(1,n);
for i = 1:length(idx)
    x_hat(idx(i)) = x_ls(i);
end
MF2=[];
%last correlation
res=A*x_hat.';
for g=1:groups
    Ai=A(:,groupind(:,g));
    MF2(g) = norm((Ai'*Ai)^(-1/2)*Ai'*res,2);
end
%calculate the sum of correlation coefficients
s=MF2;
for i=1:sparsity
    s=s+mi{i};%*((2^i-1)/31);
end
%ignore detected active channels
s(idx_g) = inf;
% sort list
[~,lmp]=sort(s,'ascend');
%remove detected ctive channels
[~,ia,~] =intersect(lmp,idx_g);
lmp(ia)=[];