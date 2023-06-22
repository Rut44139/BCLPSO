function [best,Ind_n,ConvergenceCurve,GbestScore,Fit,Div]=BPSO(Data,nX,bit,nP,Max_iteration,Ex)
%Compute nescessary variables
nXsum=nX(1)+nX(2)+nX(3); %Number of all variable
nV(1)=nX(1)*bit(1); %%Number of total bit of column  
nV(2)=nX(2)*bit(2); %Number of total bit of beam 
nV(3)=nX(3)*bit(3); %Number of total bit of brace 
nVsum=nV(1)+nV(2)+nV(3); %Number of total bit of memeber
%Initial Parameters for PSO
W=0.98;              %Inirtia weight
%Wmax=0.4;         %Max inirtia weight
%Wmin=0.9;         %Min inirtia weight
c1=2;
c2=c1;
Vmax=6;
Velocity=(-Vmax) + rand(nP,nVsum)*(2*Vmax);
%---------Cognitive component--------
PbestScore=ones(nP)*10^16;
Pbest=zeros(nP,nVsum); %Pbest of position
%-----------Social component---------
GbestScore=inf;
Gbest=zeros(1,nVsum);
%--------------------------------------------------------------
ConvergenceCurve=zeros(1,Max_iteration); %Convergence vector
%-------------------------Initialization position-------------------------%
%Random initial position 
Position=rand(nP,nVsum) < 0.5; %0.26;%Position vector
fit=0;
F1=inf;
Ind250=[]; Ind500=[]; Ind750=[];
seta_best=inf; fitness_best=inf; 
%-------------------------------BPSO loop---------------------------------%
for it=1:Max_iteration
    num2str([it, length(unique(F1)),GbestScore, fitness_best, GbestScore-fit])
    Div(it)=length(unique(F1));
    %Calculate cost for each particle
    for i=1:size(Position,1)  
    %% STRUCTURE PROBLEM
    
        %ENCODING: transfer position of particle in 1,0 to real value of x
        D_Dec=Data(0,0,0); %Only input information
        [findMem,Sc] = Decode(Position(i,:),nX,bit,D_Dec,nXsum,Ex);
        
        clear D_Dec
        
        %INPUTTING GEOMETRY SETUP
        D=Data(1,findMem,Sc);
        
        %GEOMETRY CONSTRUCTION
        [G,DL,Pfxm] = Geometry(D);
        
        %SUPPORT DOF
        [SUP] = Support(D,G); %Constant by iteration and particles
            %For Brace frame, dfr need to adjust by time
         
        %ANALYSIS
        [SL,Dis] = Stiffness(D,G,SUP,DL,Pfxm);

        %DESIGN
        [R_MF,Amp_Force,Drift]=Design(D,G,SL,Dis); %Design

        %COST FUNCTION
        [fitness,seta]=Fitness(D,R_MF,G,it,Sc);
        if i==1
            F1=[F1;fitness];
        end

%% RECORD
        %//////Input position of particle(x1,x2,..xn) in objective function//////% 
        if(PbestScore(i)>seta) & fitness-seta==0 %if fitness < Pbest, fitness will be Pbest
            PbestScore(i)=fitness;
            Pbest(i,:)=Position(i,:);
        end
        if it==1 
            Ind1=zeros(D.nMem*(D.nstory*D.nbay*2),1);
            Ind1(findMem)=Sc(:,1);   

        end  
        if seta < GbestScore 
            if fitness-seta==0 | it==1
                GbestScore=seta;
                Gbest=Position(i,:);
                fit=fitness;
                D_best=D; G_best=G;
                Prop_best=Sc; S_best=Amp_Force; Ind_best=findMem;
                Drift_best=Drift;
                fitness_best=fitness;
                Rat=R_MF;
                format long
            end
        end
    end

    %tStartUpA =cputime;  
 %% UPDATE BPSO   
    %update the W of PSO
    %W=Wmax+it*((Wmax-Wmin)/Max_iteration);
    %Update the Velocity and Position of particles
    for i=1:size(Position,1)
        for j=1:size(Position,2) 
            %Equation (1)
            Velocity(i,j)=W*Velocity(i,j)+c1*rand()*(Pbest(i,j)-Position(i,j))+c2*rand()*(Gbest(j)-Position(i,j));
            
            if(Velocity(i,j)>Vmax)
                Velocity(i,j)=Vmax;
            end
            if(Velocity(i,j)<-Vmax)
                Velocity(i,j)=-Vmax;
            end    
            s=1/(1+exp(-Velocity(i,j))); %S1 transfer function    
            if rand < s
                Position(i,j)=1;
            else
                Position(i,j)=0;
            end
        end
    end
       
        
    ConvergenceCurve(it)=GbestScore; %seta
    Fit(it)=fitness_best; 

end

    %RESTRUCT SUMMARY SUPPORT DATA
    SP.N_fx=SUP.N_fx{1}; SP.dof_fx=SUP.dof_fx{1}; SP.dfr=SUP.dfr{1}; SP.dfx=SUP.dfx{1};
    %COMBINE UNFACTOR FORCE (dof)
    %STURCT SUMMARY DATA
    D_best.SeL=D.SeL+DL;
    Ind_n=struct('Ind1',Ind1,'Ind250',Ind250,'Ind500',Ind500,'Ind750',Ind750);
    best=struct('D',D_best,'G',G_best,'P',Prop_best,'S',S_best,'SUP',SP,'Ind',Ind_best,'Drift',Drift_best,'seta',seta_best,'Rat',Rat);
end
