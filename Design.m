function [R_MF,Amp_Force,Drift]=Design(D,G,SL,Dis) %Design
%Dis1,Dis2: Unfactor displacements of DL LL SL by dofs
global alpha UBL_c UBL_b 
nstory=D.nstory;nbay=D.nbay; nMem=D.nMem; SeL=D.SeL;
Le=G.Le;ncol=D.ncol; nbe=D.nbe; nbr=D.nbr; N_cb=G.N_cb;
Ix=D.Ix; Iy=D.Iy; E=D.E; rx=D.rx; ry=D.ry; Sx=D.Sx; Zx=D.Zx; J=D.J;
bf=D.bf; tf=D.tf; d=D.d; tw=D.tw; Fy=D.Fy; Ind_N=D.Ind_N; A=D.A; 
Cw=D.Cw; Con=D.Con; DDL=D.DDL; nNode=D.nNode;
%% 1.1. Strength (Analysis)
% 0 --------------------------------DRIFT----------------------------------
u1=Dis((nstory+1)*3*(nbay)+1:3:(nstory+1)*(nbay+1)*3);
Uend=u1;
Ux=Uend(2:end);
Drift=Ux-Uend(1:end-1); 
% 1 ----------------------------- ANALYSIS --------------------------------
%1) Determine K for column (Assume K for beam and brace = 1.0)
ts=cputime;
tk=cputime;
    SM=ismember(1:nNode,Ind_N(1,:)); SM=find(SM==1);
    OM=setdiff(1:nNode,SM);          
    Gi(SM)=0;
    Gom=cell2mat(cellfun(@(x) sum(Ix(x(find(x<=ncol)))./Le(x(find(x<=ncol)))) / sum(Ix(x(find(x>ncol)))./Le(x(find(x>ncol))))  ,N_cb,'UniformOutput',false));
    Gi(OM)=Gom;     
    %Determine K (SMF); K(Unbraced)
    ebc=1:ncol; ebc=reshape(ebc,nstory,nbay+1);
    nbc=find(D.brfloor~=0);
    iub=reshape(ebc(setdiff(1:nstory,nbc),:),prod(size(ebc(setdiff(1:nstory,nbc),:))),1); ibc=reshape(ebc(nbc,:),prod(size(ebc(nbc,:))),1);
    K_MF(iub)=((1.6*Gi(Con(iub,1)).*Gi(Con(iub,2))+4*(Gi(Con(iub,1))+Gi(Con(iub,2)))+7.5)./(Gi(Con(iub,1))+Gi(Con(iub,2))+7.5)).^0.5 ; 
    %Determine K (SMF); K(Braced)
    K_MF(ibc)=(3*Gi(Con(ibc,1)).*Gi(Con(ibc,2))+1.4*(Gi(Con(ibc,1))+Gi(Con(ibc,2)))+0.64)./(3*Gi(Con(ibc,1)).*Gi(Con(ibc,2))+2*(Gi(Con(ibc,1))+Gi(Con(ibc,2)))+1.28) ; 
    K_MF(ncol+1:nMem)=1;
tek=cputime-tk;
%2) Pe1 &Pe2 divide
    %Determine Unbrace length
    Lub(1:ncol)=Le(1:ncol)*UBL_c;
    Lub(ncol+1:ncol+nbe)=Le(ncol+1:ncol+nbe)*UBL_b;
    Lub(ncol+nbe+1:nMem)=Le(ncol+nbe+1:nMem);
    %Determine Pe1 (SMF=BF); Member Consideration
    K1=1.0;
    Pe1=(pi^2*E(1)*Ix)./(K1*Lub').^2;
    %Determine Sum (Pe2); Story Consideration     Pe2 =Rm*sum(H)*L/ drift
        %LOADx = unfactor load
        %sum(H) = sum of story shear due to lateral loads (factor) 
        %drift index (drift_max/L)= 1/300
        %Rm = 0.85 for moment frame and combine system 
        %Note!! sum(H) & drift index is factor(LRFD) but they are cancelled
                %due to dividing
                
    Pe2_MF= 0.85*sum(SeL).*Le(1:nstory)./Drift; 
    %maximum load combination of horizontal factor load only is omg*QE
    %Pe2_BF= 1*omg*sum(LOAD2{3})*Le(1)./(0.02*Le(1)*[1:nstory]); 
%3) Determine Cm  (Cm for beams)
    %if m1/m2 is + ind=1
    ind=cellfun(@(x) ~ismember(max(x(3)/x(6),0),0),SL{1}(1:ncol+nbe),'UniformOutput',false);
    Cm_MF=(cellfun(@(x,y) 0.6-0.4*(-1)^(y+1) *(min(abs(x([3 6])))/max(abs(x([3 6])))) ,SL{1}(1:ncol+nbe),ind,'UniformOutput',false))';
%4) Determine B1
    N1=cell2mat(cellfun(@(x) x(1) ,SL{1}(1:ncol+nbe),'UniformOutput',false));
   %B1 >=1
    B1_MF=max(cell2mat(Cm_MF)./1-alpha*(N1'./Pe1(1:ncol+nbe)),1);
    B1_MF=[B1_MF;ones(nbr,1)];
%5) Determine B2
    %sum(Pnt) is factor load ; 
    V_L=0;
    for i=nstory:-1:1 %each nstory (kips)
        V_L=V_L + sum(DDL(ncol+i:nstory:ncol+nbe).*Le(ncol+i:nstory:ncol+nbe)); %total vertical load support by story i (kips)
        sum_Pnt(i,1)=abs(V_L);
    end
    B2=max(1./(ones(nstory,1)-alpha*(sum_Pnt./Pe2_MF')),1);
    %B2 for beam and brace =1.0
    a=ncol+1;b=ncol+nbay; B2_MF=[];
    for i=1:nstory
        B2_MF(i:nstory:ncol,1)=B2(i); %col
        B2_MF(a:b,1)=1; %beam
        a=b+1;b=a+nbay-1;
    end  
    B2_MF(1:nMem)=1;
%6) 2nd Order Forces 
    Pr_MF=cellfun(@(x1,x2,y) x1(1,:)+y*x2(1,:) ,SL{2},SL{3},num2cell(B2_MF)','UniformOutput',false);
    Mr_MF=cellfun(@(x1,x2,y1,y2) y1*x1([3,6],:)+y2*x2([3,6],:) ,SL{2}(1:ncol+nbe),SL{3}(1:ncol+nbe),num2cell(B1_MF(1:ncol+nbe))',num2cell(B2_MF(1:ncol+nbe))','UniformOutput',false);   
    Amp_Force(1:ncol+nbe,:)=[cell2mat(Pr_MF(1:ncol+nbe))' (cell2mat(Mr_MF))']; %[Pr Mr1 Mr2]
    Amp_Force(ncol+nbe+1:nMem,:)=[cell2mat(Pr_MF(ncol+nbe+1:nMem))' zeros(nbr,2)];
tes=cputime-ts;
%% 4.2. Strength (Design) not consider LB
%Design for only Y and LTB  not consider LB (C,NC,S section)
% 2 ------------------------------ DESIGN ---------------------------------    
%1) Axial Stregth
    kl_rx=K_MF'.*Le(1:nMem)./rx(1:nMem); kl_ry=1.*Lub(1:nMem)'./ry(1:nMem); kl_r_MF=max(kl_rx,kl_ry);    
   % kl_rx=K_BF.*Lub'./rx; kl_ry=1.*Lub'./ry; kl_r_BF=max(kl_rx,kl_ry);    
    %1.1) Compression member
        cm=find(Amp_Force(:,1)>0);
        %MF
        Fe=(pi^2*E(1))./(kl_r_MF(cm)).^2;
        k=Fy(cm)'./Fe; fy=Fy(cm)'; 
        Fcr(find(k<=2.25))=(0.658).^k(find(k<=2.25)).*fy(find(k<=2.25));
        Fcr(find(k>2.25))=0.877*Fe(find(k>2.25)); 
        Pnc_MF=A(cm).*Fcr';
        Pc_MF(cm)=0.85*Pnc_MF;

    %1.2) Tension member (assume no reaction member with tenslie strength)
        tm=find(Amp_Force(:,1)<=0);
        Pnt_MF=A(tm).*Fy(tm)' ; 
        Pc_MF(tm)=0.9*Pnt_MF;

%2) Flexural Strength
    %2.1 I-Section
    Fy=Fy';
    %Check Width to thickness ratio (Compact-Non-compact)
    pf=0.38*(E(1)./Fy(1)).^0.5;
    pw=3.76*(E(1)./Fy(1)).^0.5;
    %(Slender)
    rf=1*(E(1)./Fy(1)).^0.5;
    %Width to thickness ratio 
    f=bf(1:ncol+nbe)./(tf(1:ncol+nbe)*2); %flange
    w=d(1:ncol+nbe)./(tw(1:ncol+nbe)); %web
    
    fC=find(f <= pf);
    fNC=find(f > pf);
    
    wC=find(w <= pw);
    
    %-1. In-Plane Bending (Y)   Only 2 cases for interested section list
       %1.1.   f,w = C
       %(Y)
           n=intersect(fC,wC);
           Mn_I(n)=Fy(n).*Zx(n); 
       %1.2.   f = NC                w = C
       %(FLB)
           n=intersect(fNC,wC);
           Mn_I(n)=Fy(n).*Zx(n) - (Fy(n).*Zx(n) - 0.7.*Fy(n).*Sx(n)) .* (f(n)-pf)./(rf-pf); 
    %- 2. LTB   Only 1 case
        Fy=Fy';
        %--- Lb Lp Lr
        Lp=1.76.*ry(1:ncol+nbe).*(E(1:ncol+nbe)'./Fy(1:ncol+nbe)').^0.5;
        rts=sqrt(sqrt(Iy(1:ncol+nbe).*Cw(1:ncol+nbe))./Sx(1:ncol+nbe));
        c=1;  Rm=1; %Rm=1 for doubly symmetric members
        ho=d(1:ncol+nbe)-2.*(tf(1:ncol+nbe)./2);
        Lb=Lub(1:ncol+nbe)';
        Lr=1.95 .* rts(1:ncol+nbe) .* (E(1:ncol+nbe)'./(0.7.*Fy(1:ncol+nbe)')) .* (J(1:ncol+nbe).*c./(Sx(1:ncol+nbe).*ho(1:ncol+nbe)) + ((J(1:ncol+nbe).*c./(Sx(1:ncol+nbe).*ho(1:ncol+nbe))).^2 + 6.76 .* (0.7.*Fy(1:ncol+nbe)'./E(1:ncol+nbe)').^2).^0.5 ).^0.5;
        k1=find(Lb > Lp & Lb <= Lr);
        k2=find(Lb > Lr);
        k3=find(Lb <=Lp); %In case of compact web and flanges
        %--- Mn_LTB
        %2.1) Lb <=Lp
        Mn_LTB(k3)=10^12; %No LTB
        %2.2) (Lb > Lp & Lb <= Lr) and (Lb > Lr)
        k4=setdiff(1:ncol+nbe,k3); % [k1 k2]
        %Cb=12.5*Mmax/(2.5*Mmax+3*Ma+4*Mb+3*Mc)*Rm; Cb=min(Cb,3);
        W_MF=mat2cell(DDL,ones(1,nMem),1);
        W_MF=W_MF(k4);
        x=mat2cell([1/4 1/2 3/4]'*Le(k4)',3,ones(1,length(k4)));
        %-MF
        V1=cellfun(@(x) x(2,:) ,SL{1}(k4),'UniformOutput',false);
        M1=cellfun(@(x) x(3,:) ,SL{1}(k4),'UniformOutput',false);
        Mi=cellfun(@(v,m,w,x) abs(- 0.5*x.^2*w + 1.0*x*v - m) ,V1,M1,W_MF',x,'UniformOutput',false);
        Mmax=cellfun(@(x,y) max(max(x),abs(y)) ,Mi,M1,'UniformOutput',false);
        Cb=cellfun(@(mi,mx) 12.5*mx/(2.5*mx+3*mi(1)+4*mi(2)+3*mi(3))*Rm ,Mi,Mmax,'UniformOutput',false); 
        Cb=cell2mat(Cb'); uCb=zeros(ncol+nbe,1); uCb(k4)=Cb;
       
        %2.2.1) Lb > Lp & Lb <= Lr
        Mp=Fy(k1)'.*Zx(k1);
        Mn_LTB(k1)=uCb(k1).*(Mp-(Mp-0.7*Fy(k1)'.*Sx(k1)).*((Lb(k1)-Lp(k1))./(Lr(k1)-Lp(k1))));
      
        %2.2.2) Lb > Lr
        Fcr1_MF=uCb(k2).*pi^2*E(1)./(Lp(k2)./rts(k2)).^2.*(ones(length(k2),1)+0.078.*J(k2).*c./(Sx(k2).*ho(k2)).*(Lb(k2)./rts(k2)).^2).^0.5;
        Mn_LTB(k2)=Fcr1_MF.*Sx(k2);
        
    %- Mn Contorls (min(Mn_Y,Mn_LTB))
        Mn=min(Mn_I,Mn_LTB);
        Mc_MF=0.9*Mn; 
        Mc_MF(ncol+nbe+1:nMem)=inf;
    %2.2) H-SECTION (For compact section, No LB)
        %2.1. Yielding 
        %Mn_I=Fy(nMem-nbr+1:nMem)'.*Zx(nMem-nbr+1:nMem); %Zx=Plastic section modulas about x axis   
    %2.3) Struct all Mn Section
       % Mn_BF=[Mn_MF;Mn_I]; 
%3) Intereaction Equations

    P_Rt=abs(Amp_Force(:,1))./Pc_MF';
    Ratio_S=P_Rt -1;
    M_Rt= max(abs(Amp_Force(:,2)),abs(Amp_Force(:,3)))./Mc_MF';

    k1=find(P_Rt>=0.2);
    k2=find(P_Rt<0.2); 

    Ratio_Str(k1)=P_Rt(k1) + 8/9.*M_Rt(k1) -1;
    Ratio_Str(k2)=P_Rt(k2)./2 + M_Rt(k2) -1;
    %Ratio_Str(ncol+1:ncol+nbe)=-1;
    %str=Ratio_Str+1;
    %[find((str)' > 1) str(find((Ratio_Str+1)' > 1))']
   % max(str)
%% 1. Maximum Lateral Displacement (Roof)
%Uend !!!Ux at right end Node of frame (2nd order)
%col_end=(ncol-nstory+1:ncol);
%N1=Amp_Force(col_end,1);
%u1(1)=0;
%for i=2:nstory
%    u1(i)= -(1.0*(N1(i)*Le(col_end(i)) - 1.0*A(col_end(i))*E(1)*u1(i-1)))/(A(col_end(i))*E(1));
%end
%DISPLACMENT AT RIGHT END NODE
%for i=1:nstory
%    D2(i,:)=sym(['x' num2str(i)]);
%    eqn=Drift(i)+sum_Pnt(i).*D2(i)./Pe2_MF(i)==D2(i);
%    sol=vpasolve([eqn],D2(i));
%    Drift2(i)=sol;
%end
R=1/300; %max drift index
H=sum(Le(1:nstory)); %Height of frame (ft)
Ur=u1(end); 
Ratio_MaxU=abs(Ur/H)-R;
%% 1. Drift 
R=1/300; %max drift index
H=Le(1:nstory);
Uend=u1;
Ux=Uend(2:end);
Drift=Ux-Uend(1:end-1); 
Ratio_drift=abs(Drift./H)-R;

%% SUMMARY
R_MF=struct('drift',Ratio_drift','maxU',Ratio_MaxU,'Str',Ratio_Str','S',Ratio_S);
