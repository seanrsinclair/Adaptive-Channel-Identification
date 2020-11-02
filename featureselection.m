clear
clc
%rng(1)

%number of channels
groups=20;

%signal length
n=200;

%uniform channel indices in Fourier Domain
inds_A=[n-(n/groups/2-1):n 1:n/groups/2 n/groups/2+1:n-n/groups/2];
groupind=reshape(inds_A,n/groups,groups);

%loac dictionary
matload=load('dictionary_200.mat');
R=matload.A;
%DFT matrix
F=fft(eye(n))/sqrt(n);

hist_mus=[];
Rnew=R;
indices=[];
%pick up to 200 wavelets
while length(indices)<200
    %initialize mutual incoherence
    mulist=NaN(size(R,1),1);
    %loop over alll rows
    parfor i =1:size(R,1)
        %only consider the rows that are not previously selected
        if ~ismember(i,indices)
            Rtemp=[R(indices,:); (R(i,:))];
            %effective sensing matrix(sparsity is in Fourier Domain)
            A=Rtemp*F';
            mus=zeros(groups,groups);
           
            %calculate correlation between channels(groups)
            for g=1:groups
                Ai=A(:,groupind(:,g));
                for k=1:groups
                    Ak=A(:,groupind(:,k));
                    mus(g,k) = norm((Ai'*Ai)^(-1/2)*Ai'*Ak,2);
                end
            end
            cloak=eye(groups)~=1; % ones everywhere but on the diagonal
            mus=[mus.*cloak];
            %mutual coherence
            mulist(i)=max(mus(:)); % max off-diagonal term 
        end
    end
    %select the row that has the least mutual coherence and add it to
    %selected rows
    [a,minmu]=min(mulist);
    indices=[indices, minmu];
    hist_mus(end+1)=a;
    Rnew=R(indices,:);
    sprintf('Adding row %d mu_b= %d',minmu,a )
end
%plot selected features
imagesc(Rnew)




