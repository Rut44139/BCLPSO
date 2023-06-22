function Data=Input(br,findMem,Sc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RULES
%(X-bracing is count as 2 bracing)
%(Coding bases on these assumption:
 %1. Assume all bays and storys have equal length
 %2. Coordination always begins with (0,0)
 %3. All members arrangement are from bottom to top and left to right)
 %4. Seismic load only calculates by including SW
 %5. Load of all floors are equal
 %6. All supports are fixed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT DATA
%GEOMETRY
%//INPUT
Cd_bay=[0 30]; % ft
Cd_story=[0 15 [27:12:123]]; % ft
%//END

nbay=length(Cd_bay)-1; 
nstory=length(Cd_story)-1;
ncol=nstory*(nbay+1);
nbe=nstory*nbay;
if br~=0
    %GENERATE AVALIABLE BRACING INDEX
    if br==1 %Brace Frame
        All_mem=1:ncol+nbe+nstory*nbay*4;
        Ava_mem=intersect(findMem,All_mem);
        brace=find(ismember(ncol+nbe+1:ncol+nbe+nstory*nbay*4,Ava_mem)==1);
    else %Moment Frame
        brace=[]; 
    end
    
    %NODE COORDINATION
    xi=0; i=0;
    for x=Cd_bay %bay*ne
        xi=xi+1; yi=0;
        for y=Cd_story %story*ne
            yi=yi+1;
            i=i+1;
            Node(i,:)=[x y];
            Ind_N(yi,xi)=i;    
        end
    end
    nNode=length(Node); 
    
    %CONNECTION OF MEMBER (Index of node i, Index of node j))
    nc=0; nbe=0; nbr=-1; k=-3; p=0; brfloor=zeros(nstory,nbay);
    clear Con1 Con2 Con3
    Con3=[];
    for i=1:nbay+1
        for j=1:nstory+1
            if j~=nstory+1 %COLUMN
                nc=nc+1; 
                Con1(nc,:)=[Ind_N(j,i) Ind_N(j+1,i)];
            end
            if i~=nbay+1 & j~=1 %BEAM
                nbe=nbe+1; 
                Con2(nbe,:)=[Ind_N(j,i) Ind_N(j,i+1)];
            end
            if j~=nstory+1  %BRACE
                nbr=nbr+4;
            end
            if br==1
                if i~=nbay+1 & j~=nstory+1 & ismember(nbr,brace)==1
                    k=k+4; p=p+1; brfloor(j,i)=p;
                    Node(nNode+p,:)=[(Node(Ind_N(j,i),1)+Node(Ind_N(j+1,i+1),1))/2 ,(Node(Ind_N(j,i),2)+Node(Ind_N(j+1,i+1),2))/2];
                    Con3(k,:)=[Ind_N(j,i) nNode+p];
                    Con3(k+1,:)=[nNode+p Ind_N(j+1,i+1)];
                    Con3(k+2,:)=[Ind_N(j,i+1) nNode+p];
                    Con3(k+3,:)=[nNode+p Ind_N(j+1,i)];
                end
            end
        end
    end
    nNode=length(Node);
    Con=[Con1;Con2;Con3];
    nMem=length(Con);
    clear k p nbr
    nbr=length(Con3);
    
    %MATERIAL PROPERTIES
    %//INPUT
    E(1:nMem) = 29000/(1/12)^2; %ksi to kip/ft^2
    Ro= 4.573364201291581e+02; %kM/m^3 to kip/ft^3
    Fy(1:nMem)= 36/(1/12)^2; %ksi to kip/ft^2
    %//END

    %GEOMETRY PROPERTIES
    cv=1/12;
    A=Sc(1:nMem,2)* cv^2; Abr=A(ncol+nbe+1:2:nMem);
    d=Sc(1:nMem,3)* cv;
    tw=Sc(1:nMem,4)* cv;
    bf=Sc(1:nMem,5)* cv;
    tf=Sc(1:nMem,6)* cv;
    Ix=Sc(1:nMem,7)* cv^4;
    Zx=Sc(1:nMem,8)* cv^3;
    Sx=Sc(1:nMem,9)* cv^3;
    rx=Sc(1:nMem,10)* cv;
    Iy=Sc(1:nMem,11)* cv^4;
    Zy=Sc(1:nMem,12)* cv^3;
    Sy=Sc(1:nMem,13)* cv^3;
    ry=Sc(1:nMem,14)* cv;
    J=Sc(1:nMem,15)* cv^4;
    Cw=Sc(1:nMem,16)* cv^6;
    
    %EXTERNAL LOAD
    %(SL=Seismic Load[Joint Load], DDL=Dead Load[Dist. Load], DLL=Live Load[Dist. Load])  
    DLL=zeros(nMem,1); DDL=zeros(nMem,1); SeL=zeros((nstory+1)*(nbay+1)*3 +(nbr/4)*2,1); 
    %//INPUT  
    %(DLL&DDL are Frame Load Ex: DL on beam(Mem No.24) is 2 kN/m; DDL(24)=2;)
    SeL(4:3:30)=10; SeL(31) = 5; %kips
    DDL([21:29],1) = -6; %kips/ft
    DDL(30) = -3; %kips
    %//END

    Data=struct('brace',brace,'nbay',nbay,'nstory',nstory,'Node',Node,'nNode',nNode,...
        'brfloor',brfloor,'Con1',Con1,'Con2',Con2,'Con3',Con3,'Con',Con,'ncol',nc,...
        'nbe',nbe,'nbr',nbr,'nMem',nMem,'SeL',SeL,'DDL',DDL,...
        'DLL',DLL,'E',E,'A',A,'Abr',Abr,'Ix',Ix,'Iy',Iy,'Ro',Ro,'Fy',Fy,'rx',rx,'ry',ry,...
        'bf',bf,'tf',tf,'d',d,'tw',tw,'Ind_N',Ind_N,'Sx',Sx,'Sy',Sy,'Zx',Zx,'Zy',Zy,'J',J,'Cw',Cw);   %'nn',n
end
if br==0
    nMem=ncol+nbe;
    Data=struct('nbay',nbay,'nstory',nstory,'ncol',ncol,'nbe',nbe,'nMem',nMem); 
end