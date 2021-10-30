clear all;
close all;

Total_bus= 4;
zline1= [0.1414+0.5353i 0 0;0 0.1414+0.5353i 0;0 0 0.1414+0.5353i];      %line impedance matrix
zline2= zline1;                                                                 
yline= [0.09i 0 0;0 0.09i 0;0 0 0.09i];                                  %line admittance matrix

%Three phase transformer rating 12.47 kV/ 4.16 kV, 6 MVA.
ztrans= (0.01+0.06i)*(4.16^2*1000/6000);                  %transfer impedance actual value referred to the low voltage side
zt= [ztrans 0 0;0 ztrans 0; 0 0 ztrans];
V_mag= zeros(Total_bus,3);                                %storing the final voltage value to this array
Del= zeros(Total_bus,3);                                  %storing the final angle value to this array
nt= 12.47*sqrt(3)/4.16;                                   %transformer ratio


%Transformer modelling
at= -nt/3*[0 2 1;1 0 2;2 1 0];
bt= -(nt*ztrans)/3*[0 2 1;1 0 2;2 1 0];
ct= [0 0 0;0 0 0;0 0 0];
dt= (1/nt)*[1 -1 0;0 1 -1; -1 0 1];

At= (1/nt)*[1 0 -1;-1 1 0;0 -1 1];
Bt= [ztrans 0 0;0 ztrans 0;0 0 ztrans];

%line modelling
a1= [1 0 0;0 1 0;0 0 1]+ 0.5*yline*zline1;
b1= zline1;
c1= yline+0.25*yline*zline1*yline;
d1= a1;

a2= a1;
b2= zline2;
c2= c1;
d2= a2;

A1= inv(a1);
B1= inv(a1)*b1;

A2= inv(a2);
B2= inv(a2)*b2;

D= [1 -1 0;0 1 -1;-1 0 1];
s= [637.49+395.11i;637.49+395.11i;637.49+395.11i];                 % load rated power
q=1;

ELN= [7199.6;-3599.8-6235.036497i;-3599.8+6235.036497i];           %source L-N voltage
ELL= [10799.33+6235i;-12470i;-10799.33+6235i];                     % source L-L voltage which lead L-N by 30 degree
v{q+3}= [2078.461-1200i;-2078.4609-1200i;2400i];                   %voltage initialization at bus 4


k=2;                                                               %these two variables are defined to run the loop
p=3;

while k<p
    for n=1:3
        i{q+3}(n,1)= transpose(conj(s(n,1)*10^3/v{q+3}(n)));       %current at bus 4
    end

v{q+2}= a2*v{q+3}+b2*i{q+3};
i{q+2}= c2*v{q+3}+d2*i{q+3};


v{q+1}= at*v{q+2}+bt*i{q+2};
i{q+1}= ct*v{q+2}+dt*i{q+2};


v{q}= a1*v{q+1}+b1*i{q+1};
i{q}= c1*v{q+1}+d1*i{q+1};

v_ll=D*v{q};

z(:,:,1)= ELN;
z(:,:,k)=v{q};


if abs(z(:,:,k)-z(:,:,k-1))*sqrt(3)/12470<0.001                 %convergence criterion (Comparing L-N voltage at the source side)
     break
else
    v{q+1}= A1*ELN-B1*i{q+1};
    v{q+2}= At*v{q+1}-Bt*i{q+2};
    v{q+3}= A2*v{q+2}-B2*i{q+3};
    k=k+1;
    p=p+1;
end
end

for i=1:Total_bus                                               %storing the final value in this loop
    for j=1:3
    V_mag(i,j)= abs(v{1,i}(j));
    Del(i,j)= rad2deg(angle(v{1,i}(j)));
    end
end


disp('----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------');
disp('                                                       Modified Backward/forward Loadflow Analysis');
disp('----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------');
disp('==============================================================================================================================================================================================');
disp('| Bus    |Va|      |Vb|     |Vc |    |       Angle (Degree)          |');
disp('| No     |pu|      |pu|     |pu |      |Ph-a|    |Ph-b|     |Ph-c|');
disp('==============================================================================================================================================================================================');
for m=1:4
    fprintf('%4g', m); fprintf('%10.2f', V_mag(m,:)); fprintf(' %10.3f', Del(m,:));fprintf('\n');
    
end

