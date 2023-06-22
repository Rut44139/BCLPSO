function [Prc1_node,Mrc1_node]=Plot_PDELTA(bt,axis_edge)
D=bt.D; G=bt.G; S=bt.S;
nNode=D.nNode;Node=D.Node; nMem=D.nMem; N_br=D.N_br; Ind_N=D.Ind_N; Con=D.Con;
ncol=D.ncol; nbe=D.nbe; nbr=D.nbr; nstory=D.nstory;  nbay=D.nbay; ne=D.ne; Br_N=D.Br_N;
CMxi=G.CMxi; CMyi=G.CMyi; CMxj=G.CMxj; CMyj=G.CMyj; MDOF_store=G.MDOF_store; NDOF_store=G.NDOF_store;
Pr1=S(:,1); Mr1=S(:,2); Mr2=S(:,3); 
%-----------------------Separate type of dof------------------------------%
%-----------------Original coordinate of each  node-----------------------%
Cdx= Node(:,1);
Cdy= Node(:,2);
min_x=min(Cdx);     max_x=max(Cdx);
min_y=min(Cdy);     max_y=max(Cdy);
%----------------------New coordinate of each node------------------------%

%%----------------------Plot Original structure-----------------------%%
%------------------------Plot original structure--------------------------%
figure(1)
x1=[];y1=[];
for i=1:nMem
    X1(:,i)=[CMxi(i,1) CMxj(i,1)];
    Y1(:,i)=[CMyi(i,1) CMyj(i,1)];
    hold on
end


%Pr1 Mr1 Mr2 nMem
be2=[];
n1=ncol+1; n2=ncol+nstory;
for i=1:nbay
    be2=[be2 n1:n2];
    n1=n1+nstory*2; n2=n2+nstory*2;
end
Prc1_node=zeros(nNode,1);  Mrc1_node=zeros(nNode,1); %%COL
for i=1:nNode
    NDOF(i,:)=NDOF_store{i};
end

for i=1:ncol
    mdof=MDOF_store{i}; 
    node1=find(NDOF(:,1)==mdof(1));
    node2=find(NDOF(:,1)==mdof(4));
    Prc1_node(node1(1))=Pr1(i);
    Prc2_node(node2(1))=-Pr1(i);
    Mrc1_node(node1(1))=Mr1(i);
    Mrc2_node(node2(1))=Mr2(i);
end
Prbe1_node=zeros(nNode,1);  Mrbe1_node=zeros(nNode,1); %BEAM
Prbe2_node=zeros(nNode,1);  Mrbe2_node=zeros(nNode,1); 
for i=ncol+1:ncol+nbe
    if sum(i==be2) > 0
        n=find(be2==i);
        node1=Br_N(find(Br_N==Con(i,1)));
        node2=Br_N(find(Br_N==Con(i,1))+4);
        Prbe1_node(node1,1)=Pr1(ncol+n); 
        Prbe2_node(node2,1)=-Pr1(ncol+n); 
        Mrbe1_node(node1,1)=Mr1(ncol+n); 
        Mrbe2_node(node2,1)=Mr2(ncol+n);
    end
end
%[Prbe1_node Mrbe1_node Prbe2_node Mrbe2_node];
Prbr1_node=zeros(nNode,1);  Mrbr1_node=zeros(nNode,1); %BRACE
a=ncol+nbe/ne;
for i=ncol+nbe+1:nMem
    mdof=MDOF_store{i}; a=a+1;
    node1=find(NDOF(:,1)==mdof(1)); 
    node2=find(NDOF(:,1)==mdof(4));
    Prbr1_node(node1(1))=Pr1(a);
    Prbr2_node(node2(1))=-Pr1(a);
    Mrbr1_node(node1(1))=Mr1(a);
    Mrbr2_node(node2(1))=Mr2(a);
end
%%----------------Plot Support reaction of structure------------------%%
%1. INTRRNAL FORCE (COLUMN)
figure(1)
grid on
p1=plot(X1,Y1,'b','LineWidth',1.5);
title('COLUMN')
xlabel('Length(m)')
ylabel('Height(m)')
axis([min_x-axis_edge max_x+axis_edge min_y-axis_edge max_y+axis_edge]);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
hold on
%-----------------------------Plot data point-----------------------------%
%----------------------------Plot data point------------------------------%
for i=1:nNode
    x1=Cdx(i); y1=Cdy(i);
    x2=Cdx(i); y2=Cdy(i);
    Ss1=[Prc1_node(i) Mrc1_node(i)];
    Sx1=Ss1(1); Sy1=Ss1(2);
    Ss2=[Prc2_node(i) Mrc2_node(i)];
    Sx2=Ss2(1); Sy2=Ss2(2);
    Sm(i)=plot(x1,y1,'r.','MarkerSize', 25); %Dot
    textString = sprintf('(%0.2f, %0.2f)', Sx1, Sy1);
    text(x1,y1+0.2, ["Node i" textString],'Color','r', 'FontSize',10,...
           'HorizontalAlignment','center',...
           'VerticalAlignment', 'bottom');
    textString = sprintf('(%0.2f, %0.2f)', Sx2, Sy2);
    text(x2,y2-0.8,["Node j" textString],'Color','r', 'FontSize',10,...
           'HorizontalAlignment','center',...
           'VerticalAlignment', 'bottom');
end
%2. INTRRNAL FORCE (BEAM)
figure(2)
grid on
p1=plot(X1,Y1,'b','LineWidth',1.5);
title('BEAM')
xlabel('Length(m)')
ylabel('Height(m)')
axis([min_x-axis_edge max_x+axis_edge min_y-axis_edge max_y+axis_edge]);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
hold on
%-----------------------------Plot data point-----------------------------%
%----------------------------Plot data point------------------------------%
for i=1:nNode
    x1=Cdx(i); y1=Cdy(i);
    x2=Cdx(i); y2=Cdy(i);
    Ss1=[Prbe1_node(i) Mrbe1_node(i)];
    Sx1=Ss1(1); Sy1=Ss1(2);
    Ss2=[Prbe2_node(i) Mrbe2_node(i)];
    Sx2=Ss2(1); Sy2=Ss2(2);
    Sm(i)=plot(x1,y1,'r.','MarkerSize', 25); %Dot
    textString = sprintf('(%0.2f, %0.2f)', Sx1, Sy1);
    text(x1,y1+0.2, ["Node i" textString],'Color','r', 'FontSize',10,...
           'HorizontalAlignment','center',...
           'VerticalAlignment', 'bottom');
    textString = sprintf('(%0.2f, %0.2f)', Sx2, Sy2);
    text(x2,y2-0.8,["Node j" textString],'Color','r', 'FontSize',10,...
           'HorizontalAlignment','center',...
           'VerticalAlignment', 'bottom');
end
%3. INTRRNAL FORCE (BRACE)
figure(3)
grid on
p1=plot(X1,Y1,'b','LineWidth',1.5);
title('BRACE')
xlabel('Length(m)')
ylabel('Height(m)')
axis([min_x-axis_edge max_x+axis_edge min_y-axis_edge max_y+axis_edge]);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
hold on
%-----------------------------Plot data point-----------------------------%
%----------------------------Plot data point------------------------------%
for i=N_br
    x1=Cdx(i); y1=Cdy(i);
    x2=Cdx(i); y2=Cdy(i);
    Ss1=[Prbr1_node(i) Mrbr1_node(i)];
    Sx1=Ss1(1); Sy1=Ss1(2);
    Ss2=[Prbr2_node(i) Mrbr2_node(i)];
    Sx2=Ss2(1); Sy2=Ss2(2);
    Sm(i)=plot(x1,y1,'r.','MarkerSize', 25); %Dot
    textString = sprintf('(%0.2f, %0.2f)', Sx1, Sy1);
    text(x1,y1+0.2, ["Node i" textString],'Color','r', 'FontSize',10,...
           'HorizontalAlignment','center',...
           'VerticalAlignment', 'bottom');
    textString = sprintf('(%0.2f, %0.2f)', Sx2, Sy2);
    text(x2,y2-0.8,["Node j" textString],'Color','r', 'FontSize',10,...
           'HorizontalAlignment','center',...
           'VerticalAlignment', 'bottom');
end