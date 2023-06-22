function [DL,Pfxm,SL]=LOAD(AT,FT,D,G,Fs) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Analysis Type(AT): OR=1, NS=2, S=3
%OR, NS : MF & BF use the same function
%S : MF & BF use the different function
%Frame Types(FT): MF=1 (called DL LL SL), BF=2 (called only SL)
MDOF_store=G.MDOF_store; Le=G.Le; nbr=D.nbr; nstory=D.nstory; nbay=D.nbay;
%% EXTERNAL FORCE
%ANALYSIS TYPES
%1. DL & LL
if FT==2
    if AT==1 %OR NS
        DL=double(zeros((nstory+1)*(nbay+1)*3 +(nbr/4)*2,1)); 
        for i=[1:D.ncol]
            Pfxm{i}=double(zeros(6,1));
        end
        for i=[D.ncol+D.nbe+1:D.nMem]
            Pfxm{i}=double(zeros(4,1));
        end
        for i=D.ncol+1:D.ncol+D.nbe
            DL1=double(zeros((nstory+1)*(nbay+1)*3 +(nbr/4)*2,1)); 
            %DL DOF (Same direction of Load)
            MDOF=MDOF_store{i};
            DL1(MDOF(2),1)=D.DDL(i)*Le(i)/2;
            DL1(MDOF(5),1)=DL1(MDOF(2));
            DL1(MDOF(3),1)=D.DDL(i)*Le(i)^2/12; DL1(MDOF(6),1)=-DL1(MDOF(3)); 
            Pfxm{i}=-DL1(MDOF,1);% [p1 p2 p3 p4 p5 p6] dist. load on each member
            DL=DL+DL1;
        end   
        SL=D.SeL;
    end
end
%2. RS
if AT==3 %S
    NDOF_store=G.NDOF_store;
    DL=double(zeros((nstory+1)*(nbay+1)*3 +(nbr/4)*2,1));
    for i=[1:D.ncol+D.nbe]
        Pfxm{i}=double(zeros(6,1));
    end
    for i=[D.ncol+D.nbe+1:D.nMem]
        Pfxm{i}=double(zeros(4,1));
    end
    rs=zeros((nstory+1)*(nbay+1)*3 +(nbr/4)*2,1);
    for i=1:D.nstory
        NDOF=NDOF_store{(D.nbay+1)*(D.nstory+1)-D.nstory+i};
        fs=Fs;
        rs(NDOF(1),1)=-fs(end-D.nstory+i); % rs=DL LL SL
    end
    SL=rs(:,1);
end
