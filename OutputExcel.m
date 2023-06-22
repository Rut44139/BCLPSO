function SX=OutputExcel(bt,id,nx,W,ConvergenceCurve,fit,Div)
D=bt.D; Ind=bt.Ind; P=bt.P; Drt=bt.Drift; nc=D.nstory*(D.nbay+1); nbe=D.nstory*(D.nbay); nbr=D.nstory*(D.nbay)*2;
ncol=D.ncol;nbe=D.nbe;nbr=D.nbr;nMem=D.nMem;Abr=D.Abr;brace=D.brace;A=D.A; nMem=D.nMem;  nstory=D.nstory;
Av=zeros(1,nMem); Av(Ind)=P(:,2)'; Drift=zeros(1,nstory); Drift(1:nstory)= Drt*12;
Ind1000=zeros(nc+nbe+nbr,1);Ind1000(Ind)=P(:,1)';
Rat=bt.Rat; fit(1)=10^6;  ConvergenceCurve(1)=10^6;
n=1000;
SX{1}=xlswrite('Result6.xlsx',[nx-1,W,Av],'Area',['A' num2str(nx)]);
SX{2}=xlswrite('Result6.xlsx',[nx-1,fit(1:n)],'Convergence Curve1',['H' num2str(nx)]);
      %xlswrite('Result Collection6.xlsx',[nx-1,fit(10001:20000)],'Convergence Curve2',['H' num2str(nx)]);
SX{3}=xlswrite('Result6.xlsx',[nx-1,id.Ind1'],'Beginning Index',['A' num2str(nx)]);
SX{4}=xlswrite('Result6.xlsx',[nx-1,Ind1000'],'Index1000 of Avalaible Members',['A' num2str(nx)]);
SX{5}=xlswrite('Result6.xlsx',[nx-1,Drift],'Story Drift',['A' num2str(nx)]);

SX{7}=xlswrite('Result6.xlsx',[nx-1,ConvergenceCurve(1:n)],'Objective Penalty function1',['H' num2str(nx)]);
       %xlswrite('Result Collection6.xlsx',[nx-1,ConvergenceCurve(10001:20000)],'Objective Penalty function2',['H' num2str(nx)]);
SX{8}=xlswrite('Result6.xlsx',[nx-1,Rat.Str'+1],'Interaction Ratio',['H' num2str(nx)]);
SX{9}=xlswrite('Result6.xlsx',[nx-1,Rat.S'+1],'Stress Ratio',['H' num2str(nx)]);
SX{10}=xlswrite('Result6.xlsx',[nx-1,Rat.drift+1],'Drift Ratio',['H' num2str(nx)]);
SX{11}=xlswrite('Result6.xlsx',[nx-1,Div(1:n)],'Diversity',['H' num2str(nx)]);

end