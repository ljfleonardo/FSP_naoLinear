clc, clear all
s = tf('s');

%Illustrative example - 4.1 Artigo Daniel, Bruno e Julio
Ts = 0.1;
G = 1/(3*s+1);
L = exp(-2*s);
P = G*L;
C = 4.3423 * (((s+1)*(s+0.5)^2)/(s*(s^2+0.49)));

%Discretização
z = tf('z',Ts);
Ld = c2d(L,Ts);
Gd = c2d(G,Ts);
Pd = c2d(P,Ts);
Cd = c2d(C,Ts,'tustin');
F_kp1 = (0.5387*z*(z-0.9457))/(z^2-1.741*z+0.7705);
F_kp2 = (4.9238*z*(z-0.9929)*(z^2-1.93*z+0.9329))/((z-0.9951)*(z-0.7377)*(z^2-1.619*z+0.7105));

Fd = F_kp1;

% Fd.InputName = 'dyd';
% Fd.OutputName = 'dfd';
% 
% Pd.InputName = 'ud';
% Pd.OutputName = 'y0d';
% 
% Gd.InputName = 'ud';
% Gd.OutputName = 'ypd';
% 
% Ld.InputName = 'ypd';
% Ld.OutputName = 'y1d';
% 
% Cd.InputName = 'ed';
% Cd.OutputName = 'ud';
% 
% Td = feedback(Pd*Cd,1);
% Td.InputName = 'rd';
% 
% S1d = sumblk('ed = rd - ypd + dfd');
% % S2d = sumblk('yd = y0d + dd');
% S2d = sumblk('yd = y0d');
% S3d = sumblk('dyd = yd - y1d');

% t = 0:Ts:50/Ts - Ts;
% dd = 0.5*sin(0.7*t);
% 
% planta_disc = connect(Pd,Gd,Ld,Cd,Fd,S1d,S2d,S3d,{'rd'},'yd');



