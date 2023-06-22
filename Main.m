clear all;
close all;
clear clc
%% Main
mode=1; axis_edge=4; int_force=0; arrow=1.2; 
%Scale deformation
Scale=1000;
%% COMMAND
%   1     2      3    4  5   6  7   8   9   10  11  12  13  14 15  16 17 18
%[Index,Section Name, A, d, tw, bf. tf, Ix, Zx, Sx, rx, Iy, Zy,Sy, ry, J,Cw]
%tStart =cputime;
[Excel,txt]=xlsread('Section data.xlsx'); Name=txt(:,2);
Data=@(findMem,Sc,br) Input(findMem,Sc,br); 
%% Properties
global alpha ng_col ng_be UBL_c UBL_b drift_ind % e1 sol x
UBL_c=1; UBL_b=1/5; drift_ind=1/300;
alpha=1.0;

it_max=1000;
nX=[5 4 10]; %Number of unknown variable group nx=[column beam brace]
ng_col=[2 2 2 2 2]; % Members in each group
ng_be=[3 3 3 1];
bit=[7 9 7]; %2^bit possible variable/ unknown Ex col&be=93 variables=2^7; br=9 variables=2^5 (+ 1 topology control)
nP=50; %Number of particles
%tEndIM1 =cputime-tStart
for i=1:60   %Repeat times
N=i+1;
[bt,Ind_n,ConvergenceCurve,GbestScore,Fit,Div]=BPSO(Data,nX,bit,nP,it_max,Excel);

%[Prc1_node,Mrc1_node]=Plot_PDELTA(bt,axis_edge);
%p1=Plot_STRUCT(axis_edge,arrow,bt,Name);
Wend=Fit(end);
G=ConvergenceCurve-Fit;
CCV=ConvergenceCurve(find(G==0));CCV=CCV(end)
%figure
%p1=plot(ConvergenceCurve,'r-');
%p1.LineWidth=1.5;
%hold on
%grid on
%axis([0 it_max/num 60000 10^5]);
%p2=plot(64000*ones(1,it_max/num),'b-');
%p2.LineWidth=1.5;

SX=OutputExcel(bt,Ind_n,N,Wend,ConvergenceCurve,Fit,Div);
end