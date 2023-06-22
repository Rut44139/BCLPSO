function [W,seta]=Fitness(D,R_MF,GO,it,sc) 
A=D.A; Le=GO.Le; nMem=D.nMem; ncol=D.ncol; nbe=D.nbe; nbr=D.nbr;
%Parameter
P=1; e1=2;
%global sol x
%e1=subs(sol,x,it);
%Objective function
%W=Ro*(A'*Le);
W=sum(sc(:,end).*Le);
%Constriend Function
%g <=0 ; if g >0, penalty is added 
%MOMENT FRAME
G1{1}= R_MF.Str; %I
G1{2}= R_MF.S*0; %Stress
G1{3}= R_MF.drift; %Drift
%G1{3}= R_MF.maxU; %Deflection !
%G-1 <=0
G_BRACE=-nbr+1; %If no bracing, constraint is violate
G_DUMMY=length(find(sc(:,1)> 267))*10^6; %For Real section
G=0; 
for i=1:length(G1)
    g1=0; g2=0;
    g1=G1{i};
    P_G1=max(g1,0); %Find ratio >1 : ratio > 1 will not equal to 1  
    G=G+sum(P_G1); %Penalty Term without multiply penalty coefficient
end
G=G+max(G_BRACE,0);
G=G+sum(max(G_DUMMY,0));

%G=G+max(G,0)*1000;

seta=double(W*(1+P*G)^e1);
%G=vpa(G);
%R_MF.Str;
%find(G1{3}>0);
%max(R_MF.Str+1)
%seta;




