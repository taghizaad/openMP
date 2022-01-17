% clear all; close all; clc;
%----------------------------PARAMETERS------------------------------------

%Parameters init
Ts=10e-6;%time step

L=50e-3;
R=0.001;
r=3;
C = 1300e-6;

gc=C/Ts;
gl=Ts/(L+r*Ts);
m=L/(L+r*Ts);

g1=1/R;
g2=1/R;

gl1  = gl;
gl2  = gl;
gl3  = gl;
gc1  = gc;
gc2  = gc;

Ron  = 1e-3; %1e-6
Roff = 1e6; %1e9

gon = 1/Ron;
goff= 1/Roff;

nbdiode = 18;
gg = 0;

%---------------------ALL POSSIBLE DIODE STATE-----------------------------------------------------------

%all diode [anode cathode]
%Possibility Matrix(TT = Truth Table, TTd = Diode admittance)

TT = zeros(nbdiode,2^nbdiode);
ttd = zeros(nbdiode,2^nbdiode);
for  i = 1:nbdiode
if i==1
TT(i, 1:((2^nbdiode)/(2^i)))      = 0;
TT(i, 2^nbdiode/(2^i)+1 :2^nbdiode) = 1;
else
TT(i, (2^nbdiode)/(2^i)+1 :((2^i)-1)*(2^nbdiode)/(2^i)) = TT((i-1), (2^nbdiode)/(2^(i-1))+1 :2^nbdiode);
TT(i, :) = xor(TT(i,:), TT(i-1,:));
end
end

for  i = 1:2^nbdiode
    for  j = 1:nbdiode
        if TT(j,i) ==1
            ttd(j,i) = gon;
        else
            ttd(j,i) = goff; 
        end
    end
end

%-----------------------ALL INVERSE TO 3 LEVEL INVERTER------------------------------------------------
up  = 1;
up2 = 21;
HH = zeros(length(ttd)*21,21);
IH = zeros(length(ttd)*21,21);

for i=1:length(ttd)
s1 = ttd(1,i);
s2 = ttd(2,i);  
s3 = ttd(3,i); 
s4 = ttd(4,i);  
s5 = ttd(5,i); 
s6 = ttd(6,i);  
s7 = ttd(7,i);
s8 = ttd(8,i);  
s9 = ttd(9,i); 
s10 = ttd(10,i);  
s11 = ttd(11,i); 
s12 = ttd(12,i);  
gd1 = ttd(13,i);
gd2 = ttd(14,i);  
gd3 = ttd(15,i); 
gd4 = ttd(16,i);  
gd5 = ttd(17,i); 
gd6 = ttd(18,i);  

 H = ...
[ g1,   0,                     -g1,                        0,             0,             0,             0,             0,             0,             0,              0,               0,               0,    0,    0,    0, 1, 0, 0, 0, 0;
   0,  g2,                       0,                      -g2,             0,             0,             0,             0,             0,             0,              0,               0,               0,    0,    0,    0, 0, 1, 0, 0, 0;
 -g1,   0, g1 + gc1 + s1 + s5 + s9,                        0,           -s1,             0,             0,           -s5,             0,             0,            -s9,               0,               0,    0,    0,    0, 0, 0, 0, 0, 0;
   0, -g2,                       0, g2 + gc2 + s4 + s8 + s12,             0,             0,           -s4,             0,             0,           -s8,              0,               0,            -s12,    0,    0,    0, 0, 0, 0, 0, 0;
   0,   0,                     -s1,                        0, gd1 + s1 + s2,           -s2,             0,             0,             0,             0,              0,               0,               0,    0,    0,    0, 0, 0, 0, 0, 0;
   0,   0,                       0,                        0,           -s2, gl1 + s2 + s3,           -s3,             0,             0,             0,              0,               0,               0, -gl1,    0,    0, 0, 0, 0, 0, 0;
   0,   0,                       0,                      -s4,             0,           -s3, gd2 + s3 + s4,             0,             0,             0,              0,               0,               0,    0,    0,    0, 0, 0, 0, 0, 0;
   0,   0,                     -s5,                        0,             0,             0,             0, gd3 + s5 + s6,           -s6,             0,              0,               0,               0,    0,    0,    0, 0, 0, 0, 0, 0;
   0,   0,                       0,                        0,             0,             0,             0,           -s6, gl2 + s6 + s7,           -s7,              0,               0,               0,    0, -gl2,    0, 0, 0, 0, 0, 0;
   0,   0,                       0,                      -s8,             0,             0,             0,             0,           -s7, gd4 + s7 + s8,              0,               0,               0,    0,    0,    0, 0, 0, 0, 0, 0;
   0,   0,                     -s9,                        0,             0,             0,             0,             0,             0,             0, gd5 + s9 + s10,            -s10,               0,    0,    0,    0, 0, 0, 0, 0, 0;
   0,   0,                       0,                        0,             0,             0,             0,             0,             0,             0,           -s10, gl3 + s10 + s11,            -s11,    0,    0, -gl3, 0, 0, 0, 0, 0;
   0,   0,                       0,                     -s12,             0,             0,             0,             0,             0,             0,              0,            -s11, gd6 + s11 + s12,    0,    0,    0, 0, 0, 0, 0, 0;
   0,   0,                       0,                        0,             0,          -gl1,             0,             0,             0,             0,              0,               0,               0,  gl1,    0,    0, 0, 0, 1, 0, 0;
   0,   0,                       0,                        0,             0,             0,             0,             0,          -gl2,             0,              0,               0,               0,    0,  gl2,    0, 0, 0, 0, 1, 0;
   0,   0,                       0,                        0,             0,             0,             0,             0,             0,             0,              0,            -gl3,               0,    0,    0,  gl3, 0, 0, 0, 0, 1;
   1,   0,                       0,                        0,             0,             0,             0,             0,             0,             0,              0,               0,               0,    0,    0,    0, 0, 0, 0, 0, 0;
   0,   1,                       0,                        0,             0,             0,             0,             0,             0,             0,              0,               0,               0,    0,    0,    0, 0, 0, 0, 0, 0;
   0,   0,                       0,                        0,             0,             0,             0,             0,             0,             0,              0,               0,               0,    1,    0,    0, 0, 0, 0, 0, 0;
   0,   0,                       0,                        0,             0,             0,             0,             0,             0,             0,              0,               0,               0,    0,    1,    0, 0, 0, 0, 0, 0;
   0,   0,                       0,                        0,             0,             0,             0,             0,             0,             0,              0,               0,               0,    0,    0,    1, 0, 0, 0, 0, 0];

%Inverse computation
HH(up:up2,:) = H;
IH(up:up2,:) = inv(H);
up  = up +21;
up2 = up2+21;
end

