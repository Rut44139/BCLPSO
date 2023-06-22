%/////////////////////////////////////////////////////////////////////////%
%--This function is Input Function that contains Load and Support Inputs--%
% Load&Support indexs are used to classify type of analysis as follows:   %
% 1 = original elastic analysis                                           %
% 2 = Non-sway case analysis                                              %
% 3 = Sway case analysis                                                  %
%/////////////////////////////////////////////////////////////////////////%
%   OR NS S                           
%SP=[1 1  1  << fixed dof                
%    5 5  5                             
%    9 9  9
%   13 13 13
%   17 17 17
%   21 21 21
%    0 22  0
%    0 23  0
%    0 24  0]
function [SUP,LOAD]=Support(D,G)
nNode=D.nNode; Con2=D.Con2; Ind_N=D.Ind_N; nMem=D.nMem; NDOF_store=G.NDOF_store;
nstory=D.nstory; nbay=D.nbay;
%% SUPPORT

for i=1:3 %(OR NS S)
    if i == 1 | i ==3 
       n=find(Ind_N(1,:)~=0);
       N_fx{i}=Ind_N(1,n);  n_fx=N_fx{i};
       dof_fx{i}=zeros(length(n_fx),3); %(Type of DOF (Free=1,Fix=0)of Support's nodes) 
       %(Ex: There are 2 fixed support's nodes; dof_fx=[0 0 0; 0 0 0])
    elseif i == 2
        N_roller=[max(n_fx)+1:max(n_fx)+nstory]; %Nodes which have additional roller supports fixed in x-axis
        NN_roll=length(N_roller); %No. of N_roller node
        N_fx{i}=[n_fx N_roller]; %Support Node
        dof_fx{i}=[dof_fx{i-1}; [zeros(NN_roll,1) ones(NN_roll,1) ones(NN_roll,1)]];
        n_fx=N_fx{i};
    end
    ST=ones((nstory+1)*(nbay+1),3);
    ndof=[];
    for j=1:length(n_fx)
        ndof(j,:)=NDOF_store{n_fx(j)};
    end
    ST(ndof)=dof_fx{i};
    dfx{i}= find(ST==0);
    dfr{i} = find(ST(1:(nstory+1)*(nbay+1),:)==1);
    dfr{i} = [dfr{i}; ((nstory+1)*(nbay+1)*3+1: (nstory+1)*(nbay+1)*3+(nNode-(nstory+1)*(nbay+1))*2)'];
end
SUP=struct('N_fx',{N_fx},'dof_fx',{dof_fx},'dfr',{dfr},'dfx',{dfx}); %Follow analysis types
