function [SL,Dis] = Stiffness(D,G,SUP,DL,Pfxm)
A=D.A; E=D.E; I=D.Ix;
dfx=SUP.dfx; dfr=SUP.dfr; nNode=D.nNode; nstory=D.nstory; nbay=D.nbay;
ncol=D.ncol; nbe=D.nbe; SeL=D.SeL; nMem=D.nMem;
Le=G.Le; c=G.c; s=G.s; MDOF_store=G.MDOF_store;

%% PART1 Stiffnes analysis
%STIFFNESS MATRIX
%ts=cputime;
%6 MDOF (FRAME)
T_store=cellfun(@(c,s) [ c  -s  0  0     0   0;
                         s   c  0  0     0   0;
                         0   0  1  0     0   0;
                         0   0  0  c    -s   0;
                         0   0  0  s     c   0;
                         0   0  0  0     0   1] , num2cell(c(1:ncol+nbe)),num2cell(s(1:ncol+nbe)), 'UniformOutput',false); 
%4 MDOF (TRUSS)
T_store=[T_store cellfun(@(c,s) [ c  -s    0     0   ;
                         s   c    0     0   ;

                         0   0    c    -s   ;
                         0   0    s     c   ], num2cell(c(ncol+nbe+1:end)),num2cell(s(ncol+nbe+1:end)), 'UniformOutput',false)]; 

E=E(1);
%6 MDOF (FRAME)
kL_store =cellfun(@(A,L,I)   [A*E/L              0                          0                 -A*E/L                  0                         0;
                                0      12/L^2 *E*I/(L)        6/L *E*I/(L)           0         -12/L^2 *E*I/(L)       6/L *E*I/(L);
                                0        6/L *E*I/(L)        (4) *E*I/(L)          0          -6/L *E*I/(L)        (2) *E*I/(L) ;
                              -A*E/L             0                           0                 A*E/L                  0                         0;
                                0     -12/L^2 *E*I/(L)       -6/L *E*I/(L)           0          12/L^2 *E*I/(L)      -6/L *E*I/(L);
                                0        6/L *E*I/(L)        (2) *E*I/(L)          0          -6/L  *E*I/(L)       (4) *E*I/(L)]  , num2cell(A(1:ncol+nbe)),num2cell(Le(1:ncol+nbe)),num2cell(I(1:ncol+nbe)), 'UniformOutput',false); 

%4 MDOF (TRUSS)
kL_store =[kL_store; cellfun(@(A,L) [A*E/L             0               -A*E/L                  0     ;
                                        0               0                  0                    0     ;
                               
                                     -A*E/L             0                A*E/L                  0     ;
                                        0               0                  0                    0       ]  , num2cell(A(ncol+nbe+1:end)),num2cell(Le(ncol+nbe+1:end)), 'UniformOutput',false)]; 

kG_store=cellfun(@(T,k)  T*k*T' , T_store',kL_store, 'UniformOutput',false);     
nNode_MF=(nbay+1)*(nstory+1);
K_BF=zeros(nNode_MF*3 + (nNode-nNode_MF)*2,nNode_MF*3 + (nNode-nNode_MF)*2); 
for i=1:nMem
    K_BF(MDOF_store{i},MDOF_store{i})= K_BF(MDOF_store{i},MDOF_store{i})+ kG_store{i};
end

tk=cputime;
for i=1:3 %OR NS S
    if i==3
        [DL,Pfxm,SeL]=LOAD(3,2,D,G,Rs); %called NS reaction MF
    end
    %--Stiffness Matrix
    Kff_BF= K_BF(dfr{i},dfr{i});
    Krf_BF= K_BF(dfx{i},dfr{i});
    %--Force
    %DL SeL 
    %--Displacement U
    %(Unfactor)
    Uf_BF   =   Kff_BF^-1   *  (DL(dfr{i})+SeL(dfr{i}));
    U_BF=zeros(nNode*3,1); U_BF(dfr{i})=Uf_BF;
    if i==1
        Dis=[U_BF];
    end
    Fs= double(Krf_BF*Uf_BF-(DL(dfx{i})));
    if i==2
        Rs=Fs;
    end
    %DEFINE VG OF EACH MEMBER
     vG=cellfun(@(x)double(U_BF(x)),MDOF_store, 'UniformOutput',false);  
    %CALCULATE INTERNAL FORCE
     SL{:,i}= cellfun(@(x,y,z,w) inv(w)*((x*y)+z), kG_store',vG',Pfxm,T_store, 'UniformOutput',false); 
end
