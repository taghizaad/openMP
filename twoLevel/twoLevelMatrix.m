% close all
% clear all
% clc

%Parameters of the model
Ts = 50e-6;
L  = 5e-3;
R  = 0.01;
gr = 1/R;

%discrete inductance
gl=Ts/L;

%parameter assignation
gl1=gl;
gl2=gl;
gl3=gl;
gr1=gr;
gr2=gr;
gr3=gr;

%switch parameter on/off
Ron  = 1e-3;
Roff = 1e9;

%switch value within the matrix H
gon = 1/Ron;
goff= 1/Roff;

nbdiode = 6;


%Possibility Matrix(TT = Truth Table, TTd = Diode admittance)

TT = zeros(nbdiode,2^nbdiode);

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


%% all H matrix possible variation
up  = 1;
up2 = 16;
for i=1:length(ttd)
    g1 = ttd(1,i);
    g2 = ttd(2,i);
    g3 = ttd(3,i);
    g4 = ttd(4,i);
    g5 = ttd(5,i);
    g6 = ttd(6,i);
    
    H=[g1 + g3 + g5,            0,           -g1,           -g3,           -g5,         0,         0,         0,    0,    0,    0, 1, 0, 0, 0, 0;
        0, g2 + g4 + g6,           -g2,           -g4,           -g6,         0,         0,         0,    0,    0,    0, 0, 1, 0, 0, 0;
        -g1,          -g2, g1 + g2 + gr1,             0,             0,      -gr1,         0,         0,    0,    0,    0, 0, 0, 0, 0, 0;
        -g3,          -g4,             0, g3 + g4 + gr2,             0,         0,      -gr2,         0,    0,    0,    0, 0, 0, 0, 0, 0;
        -g5,          -g6,             0,             0, g5 + g6 + gr3,         0,         0,      -gr3,    0,    0,    0, 0, 0, 0, 0, 0;
        0,            0,          -gr1,             0,             0, gl1 + gr1,         0,         0, -gl1,    0,    0, 0, 0, 0, 0, 0;
        0,            0,             0,          -gr2,             0,         0, gl2 + gr2,         0,    0, -gl2,    0, 0, 0, 0, 0, 0;
        0,            0,             0,             0,          -gr3,         0,         0, gl3 + gr3,    0,    0, -gl3, 0, 0, 0, 0, 0;
        0,            0,             0,             0,             0,      -gl1,         0,         0,  gl1,    0,    0, 0, 0, 1, 0, 0;
        0,            0,             0,             0,             0,         0,      -gl2,         0,    0,  gl2,    0, 0, 0, 0, 1, 0;
        0,            0,             0,             0,             0,         0,         0,      -gl3,    0,    0,  gl3, 0, 0, 0, 0, 1;
        1,            0,             0,             0,             0,         0,         0,         0,    0,    0,    0, 0, 0, 0, 0, 0;
        0,            1,             0,             0,             0,         0,         0,         0,    0,    0,    0, 0, 0, 0, 0, 0;
        0,            0,             0,             0,             0,         0,         0,         0,    1,    0,    0, 0, 0, 0, 0, 0;
        0,            0,             0,             0,             0,         0,         0,         0,    0,    1,    0, 0, 0, 0, 0, 0;
        0,            0,             0,             0,             0,         0,         0,         0,    0,    0,    1, 0, 0, 0, 0, 0];
    
    %Inverse
    HH(up:up2,:) = H;
    IH(up:up2,:) = inv(H);
    
    up  = up +16;
    up2 = up2+16;
end
save('H','HH')
save('IH','IH')