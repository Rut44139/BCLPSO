function [Geo,DL,Pfxm]=Geometry(D)
Node=D.Node; Con=D.Con; nMem=D.nMem; nNode=D.nNode; Ind_N=D.Ind_N; ncol=D.ncol;
nNode_MF=(D.nstory+1)*(D.nbay+1); nbe=D.nbe;
%% SET UP DOF&Coordinate
% 1.-(COORDINATE OF EACH MEBER)
for el=1:nMem
    CMxi(el,1)= Node(Con(el,1),1);
    CMyi(el,1)= Node(Con(el,1),2);
    CMxj(el,1)= Node(Con(el,2),1); 
    CMyj(el,1)= Node(Con(el,2),2);
end
[CMxi CMxj CMyi CMyj];
% 2.- (DOF OF EACH NODE& DOF SUPPORT TYPE)
a=1;b=3; 
for el=1:nNode_MF
    NDOF_store{el}= [a:b]; 
    a=a+3;b=b+3;
end
a=nNode_MF*3+1;b=a+1; 
for el=nNode_MF+1:nNode
    NDOF_store{el}= [a:b]; 
    a=a+2; b=b+2;
end
% 3.-(DOF OF EACH MEMBER)
for el=1:ncol+nbe
    MDOF_store{el,1}= [NDOF_store{Con(el,1)} NDOF_store{Con(el,2)}];    
end
for el=ncol+nbe+1:nMem
    MDOF_store{el,1}= [NDOF_store{Con(el,1)}(1:2) NDOF_store{Con(el,2)}(1:2)];    
end
% 4.-(Length, cos and sin)
for i=1:nMem
    L= (((CMxj(i)-CMxi(i))^2) + ((CMyj(i)-CMyi(i))^2))^0.5;
    Le(i,1)=L;
    cos= (CMxj(i)-CMxi(i))/L; c(i)=cos;
    sin= (CMyj(i)-CMyi(i))/L; s(i)=sin;
end

% 5.-Properties to find K (effective length factor)
%Find support nodes
SM=ismember(1:nNode,Ind_N(1,:)); SM=find(SM==1);
%Find other nodes
OM=setdiff(1:nNode,SM);
%Find GA ; if OM; sum([I/L]c / [I/L]b)
              %else; 0
N_cb=cellfun(@(x) find(Con(1:nMem,1)==x | Con(1:nMem,2)==x),num2cell(OM),'UniformOutput',false);
    
Geo=struct('CMxi',CMxi,'CMyi',CMyi,'CMxj',CMxj,'CMyj',CMyj,'NDOF_store',{NDOF_store},...
    'MDOF_store',{MDOF_store},'Le',Le,'c',c,'s',s,'N_cb',{N_cb});
%p1=Plot_STRUCT_TEST(D,Geo);
%Geo;
[DL,Pfxm]=LOAD(1,2,D,Geo);
end