function [AR,RI,MI,HI,rowSplits,columnSplits,typeConfusions] = reportConfusionsAndRI(c1,c2)

% This code is slightly modified from the original by Uygar Sümbül.
% The copyright information for the original file is generated below.
% The original file can be downloaded from the following website as of Dec. 23, 2013:
% http://www.mathworks.com/matlabcentral/fileexchange/13916-simple-tool-for-estimating-the-number-of-clusters/content/valid_RandIndex.m

%RANDINDEX - calculates Rand Indices to compare two partitions
% ARI=RANDINDEX(c1,c2), where c1,c2 are vectors listing the
% class membership, returns the "Hubert & Arabie adjusted Rand index".
% [AR,RI,MI,HI]=RANDINDEX(c1,c2) returns the adjusted Rand index,
% the unadjusted Rand index, "Mirkin's" index and "Hubert's" index.
%
% See L. Hubert and P. Arabie (1985) "Comparing Partitions" Journal of
% Classification 2:193-218

%(C) David Corney (2000) D.Corney@cs.ucl.ac.uk

if nargin < 2 | min(size(c1)) > 1 | min(size(c2)) > 1
   error('RandIndex: Requires two vector arguments')
   return
end

C=Contingency(c1,c2); %form contingency matrix

n=sum(sum(C));
nis=sum(sum(C,2).^2); %sum of squares of sums of rows
njs=sum(sum(C,1).^2); %sum of squares of sums of columns

t1=nchoosek(n,2); %total number of pairs of entities
t2=sum(sum(C.^2)); %sum over rows & columnns of nij^2
t3=.5*(nis+njs);

%Expected index (for adjustment)
nc=(n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1));

A=t1+t2-t3; %no. agreements
D= -t2+t3; %no. disagreements

if t1==nc
   AR=0; %avoid division by zero; if k=1, define Rand = 0
else
   AR=(A-nc)/(t1-nc); %adjusted Rand - Hubert & Arabie 1985
end

RI=A/t1; %Rand 1971 %Probability of agreement
MI=D/t1; %Mirkin 1970 %p(disagreement)
HI=(A-D)/t1; %Hubert 1977 %p(agree)-p(disagree)

tmp1=sum(C>0,1); tmp2=sum(C>0,2);
rowSplits = sum(tmp1)-nnz(tmp1); columnSplits = sum(tmp2)-nnz(tmp2);
typeConfusions = rowSplits+columnSplits; % number of splits in both directions!
%for kk=1:size(C,2)
% [maxi,pos]=max(C(:,kk)); if pos>kk; tmp = C(pos,:); C(pos,:) = C(kk,:); C(kk,:) = tmp; end;
%end
%tmp = diag(diag(C)); tmp = [tmp; zeros(size(C,1)-size(tmp,1),size(tmp,2))]; tmp = [tmp zeros(size(tmp,1), size(C,2)-size(tmp,2))]; tmp = C-tmp;
%typeConfusions = nnz(tmp);

function Cont=Contingency(Mem1,Mem2)

if nargin < 2 | min(size(Mem1)) > 1 | min(size(Mem2)) > 1
   error('Contingency: Requires two vector arguments')
   return
end

Cont=zeros(max(Mem1),max(Mem2));

for i = 1:length(Mem1);
   Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
end
