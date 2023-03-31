clear; close all; clc
s = tf('s');

import casadi.*
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x3 = SX.sym('x3');

u1 = SX.sym('u1');
u2 = SX.sym('u2');
u3 = SX.sym('u3');

qi = SX.sym('qi');
Ti = SX.sym('Ti');

d_min = 0;                            %Atraso m�nimo usado para predi��o
d = 10;                               %Atraso real
delay_real     = [d d d];             %Vetor de atraso real
delay_modelado = delay_real;          %Vetor de atraso modelado
d_max          = max(delay_real);     %M�ximo valor do atraso (caso sejam diferentes)
dmodelado_max  = max(delay_modelado); %M�ximo valor do atraso (caso sejam diferentes)
delay_total    = d_max+dmodelado_max; %Soma dos maiores atrasos para predi��o

%% --------------- CSTR MIMO inst�vel ---------------
global Ts Ac dH p cp ER k0 V
%---- Constantes ----
Ac = 0.05;        %Area do reator
dH = -50e6;       %Calor da rea��o
p  = 1e6;         %Densidade dos l�quidos
cp = 1;           %Calor espec�fico
ER = 1000;        %Termo de energia de ativa��o
k0 = 4.872997584; %Constante da taxa de rea��o
V  = 0.05;        %Volume do reator
Ts = 3;           %Per�odo de amostragem

%---- Ponto de Opera��o ----
%Entradas / Vari�veis Manipuladas
% u1 = q0     - Vaz�o de sa�da;
% u2 = Caf    - Concentra��o do produto A na alimenta��o do tanque; 
% u3 = Qh/pCp - taxa de remo��o de calor normalizada;
qo0  = 5e-3;
Caf0 = 5;
Qh0  = 0.75;

u0 = [qo0; Caf0; Qh0];

%Perturba��es
% qi = qi - Vaz�o de entrada; 
% Ti = Ti - Temperatura externa;
qi0 = 5e-3;
Ti0 = 350;

q0  = [qi0; Ti0];

% Sa�das
%Estados / Vari�veis Controladas
% x1 = h - n�vel dentro do tanque; 
% x2 = Ca - Concentra��o de sa�da do produto A;
% x3 = T  - Temperatura dentro do reator; 
h0  = 1;
Ca0 = 1;
T0 = 400;

x0  = [h0; Ca0; T0];

%---- Modelo ----

Rt = (k0*exp(-ER/x3));

dx1_discreto = x1 + Ts*((qi-u1)/Ac);
dx2_discreto = x2 + Ts*(qi*((u2-x2)/V) - Rt*x2);
dx3_discreto = x3 + Ts*(qi*((Ti-x3)/V) - (dH/p/cp)*Rt*x2 - u3/V);

fun_ax_ext = [dx1_discreto;dx2_discreto;dx3_discreto;qi;Ti];     %Sistema aumentado discreto
fun_x      = Function('fun_x',{x1,x2,x3,qi,Ti,u1,u2,u3},{fun_ax_ext});

fun_yx_ext = [x1;x2;x3];
fun_y      = Function('fun_y',{x1,x2,x3,qi,Ti,u1,u2,u3},{fun_yx_ext});

na = size(fun_ax_ext,1);                      %Dimens�o vetor de estados aumentado
m  = size(fun_yx_ext,1);                      %Dimens�o vetor de sa�das
%% --------------- Inicializa��o das vari�veis --------------------
iteracoes   = 900;

%% ---- Sa�das ----
x       = zeros(3,iteracoes);          %inicia vetor de sa�das
x(1,:)  = x0(1);
x(2,:)  = x0(2);
x(3,:)  = x0(3);

%---- Entradas ----
entradas      = zeros(3,iteracoes);    %inicia vetor de entradas
entradas(1,:) = u0(1);
entradas(2,:) = u0(2);
entradas(3,:) = u0(3); 

%---- Perturba��es ----
perturbacoes      = zeros(2,iteracoes); %inicia vetor de perturba��es
perturbacoes(1,:) = q0(1);                   
perturbacoes(2,:) = q0(2);

%---- Iniciliza��o das refer�ncias de controle ----
ref_h            = zeros(1,iteracoes);
ref_h(1:end)     = x(1,1);

ref_ca           = zeros(1,iteracoes);
ref_ca(1:end)    = x(2,1);

ref_T            = zeros(1,iteracoes);
ref_T(1:end)     = x(3,1);

%---- Vari�veis De Estima��o -----
e           = zeros(na-2,iteracoes); %Erro de estima��o
ef          = zeros(na-2,iteracoes); %Erro de estima��o filtrado
x_n         = zeros(na-2,iteracoes); %Sa�das nominais
x_n(:,1)    = x(:,1);                %Saidas nominais anteriores
x_pred_vect = x;                     %Vetor de sa�das preditas usado para plot

%---- Fun��es CasADi ----
A_jacobian = jacobian(fun_ax_ext,[x1 x2 x3 qi Ti]);  %C�lculo do jacobiano para matriz A
B_jacobian = jacobian(fun_ax_ext,[u1 u2 u3]);        %C�lculo do jacobiano para matriz B
C_jacobian = jacobian(fun_yx_ext,[x1 x2 x3 qi Ti]);  %C�lculo do jacobiano para matriz C
fun_a = Function('fun_a',{x1,x2,x3,qi,Ti,u1,u2,u3},{A_jacobian});
fun_b = Function('fun_b',{x1,x2,x3,qi,Ti,u1,u2,u3},{B_jacobian});
fun_c = Function('fun_a',{x1,x2,x3,qi,Ti,u1,u2,u3},{C_jacobian});

%% --------------- Projeto de Controle ---------------
for lqr=1
%---- Lineariza��o no ponto de opera��o ----
Aa = full(fun_a(x(1,1),x(2,1),x(3,1),perturbacoes(1,1),perturbacoes(2,1),entradas(1,1),entradas(2,1),entradas(3,1)));
Ba = full(fun_b(x(1,1),x(2,1),x(3,1),perturbacoes(1,1),perturbacoes(2,1),entradas(1,1),entradas(2,1),entradas(3,1)));
Ca = full(fun_c(x(1,1),x(2,1),x(3,1),perturbacoes(1,1),perturbacoes(2,1),entradas(1,1),entradas(2,1),entradas(3,1)));
Aa_x = Aa(1:na-2,1:3);
Ba_x = Ba(1:na-2,1:3);
Ca_x = Ca(1:na-2,1:3);

%---- Insere din�mica integradora no controlador ----
A_exp = [Aa_x zeros(size(Aa_x,1),size(Ca_x,1));
         -Ca_x         [1 0 0; 0 1 0;0 0 1]];
B_exp = [Ba_x;zeros(3,3)];
C_exp = [Ca_x zeros(3,3)];
D_exp = 0;
sistema = ss(A_exp,B_exp,C_exp,D_exp,Ts);                            %Sistema expandido

Co = ctrb(sistema.A,sistema.B);                                      %Matriz de Controlabilidade
rank(Co);
Q_lqr = diag([1,1,1/400^2,1,1,1/400^2]);                             %Matrizes de pondera��o
R_lqr = 100*diag([1/0.005^2,1/5^2,1/0.75^2]);
[K,P,poles] = dlqr(sistema.A,sistema.B,Q_lqr,R_lqr);                 %Ganho do controlador via dlqr

%---- Inicializa integrador no ponto de opera��o ---- 
integrador_anterior = -pinv(K(:,4:6))*(entradas(:,1)+K(:,1:3)*x);

integrador_anterior_1 = integrador_anterior(1);
integrador_atual_1    = 0;

% integrador_anterior_2 = (entradas(:,1)+K(1:3,1:3)*x_estimado)/(-K(1:3,5));
integrador_anterior_2 = integrador_anterior(2);
integrador_atual_2    = 0;
% end

integrador_anterior_3 = integrador_anterior(3);
integrador_atual_3    = 0;
end

%Vari�vel para fechar a malha de controle; 0 = malha aberta; 1 = malha fechada
controle = 1;

%% --------------- Defini��o das refer�ncias ---------------
if controle == 1 % --- MALHA FECHADA ---
    
%     ref_h(200:end)     = x(1,1)*1.02;          
    ref_ca(100:end)    = x(2,1)*1.02;            
%     ref_T(600:end)     = x(3,1)*1.02;
%     perturbacoes(1,200:end)      = q0(1)*1.02;    
%     perturbacoes(2,400:end)      = q0(2)*0.98;

else % --- MALHA ABERTA ---
%     entradas(1,200:end) = u0(1)*1.02;             
    entradas(2,400:end) = u0(2)*1.02;               %Degrau na entrada de 2%
%     entradas(3,600:end) = u0(3)*1.02;    
    perturbacoes(1,100:end)      = q0(1)*1.02;    
%     perturbacoes(2,700:end)      = q0(2)*1.02;
end
%% --------------- Projeto do Filtro ---------------

%melhroara / ver como projetar o filtro
lambda = 1e-5;
% lambda_d = 100;
F        = (s*lambda)/((s+lambda));
% z = tf('z',Ts);
% F = tf(zpk([-0.1],[-0.1],[1]));
Fd       = c2d(F,Ts);
% Fd = (z*lambda)/(z-lambda);
den_f    = Fd.den{1,1};
num_f    = Fd.num{1,1};
n        = length(den_f(2:end)); %Ordem do filtro

%% --------------- Simula��o --------------------
tic
for k = 2+delay_total:iteracoes
    %% ---- Simula��o do processo ----
    %Compila as entradas com atraso real
    entrada_atrasada_real   = [entradas(1,k-1-delay_real(1)) entradas(2,k-1-delay_real(2)) entradas(3,k-1-delay_real(3))];
    
    x(:,k) = modeloCSTR_naoIsotermico(x(:,k-1)', entrada_atrasada_real, perturbacoes(:,k-1)');
    
    %% ---- Predi��o ----
    %Compila as entradas com atraso pequeno
    entrada_atraso_minimo = [entradas(1,k-1-d_min) entradas(2,k-1-d_min) entradas(3,k-1-d_min)];
    
    %---- Simula��o do sistema sem atraso (atraso min) ----
    %estados nominais, entrada sem atraso e perturba��o no ponto de opera��o
    x_n(:,k) = modeloCSTR_naoIsotermico(x_n(:,k-d_max-1), entrada_atraso_minimo, q0);
    
    %---- Filtro ----
    e(:,k)   = x(:,k) - x_n(:,k);                                      %Erro de predi��o = saida predita - sa�da real
    ef(:,k)  = -den_f(2:end)*ef(:,k-1:-1:k-n)' + num_f*e(:,k:-1:k-n)'; %Filtragem do sinal do erro
    %---- Sa�da do preditor ----
    x_pred   = ef(:,k) + x_n(:,k);    
    
    %% ---- Controle ----
    if controle == 1
        %Integrador para altura do tanque - h
        erro_1                = ref_h(k) - x_pred(1);
        integrador_atual_1    = integrador_anterior_1 + erro_1;
        integrador_anterior_1 = integrador_atual_1;
            
        %Integrador para concentra��o do produto A - Ca
        erro_2                = ref_ca(k) - x_pred(2);
        integrador_atual_2    = integrador_anterior_2 + erro_2;
        integrador_anterior_2 = integrador_atual_2;
        
        %Integrador para temperatura interna - T
        erro_3                = ref_T(k) - x_pred(3);
        integrador_atual_3    = integrador_anterior_3 + erro_3;
        integrador_anterior_3 = integrador_atual_3;
        
        nova_entrada  = -K*[x_pred(1);x_pred(2);x_pred(3);integrador_atual_1;integrador_atual_2;integrador_atual_3];
        entradas(:,k) = nova_entrada;
    end
    
    %---- Atualiza��o das vari�veis ----
    x_n_passado          = x_n(:,k-d_max);     %Sa�das passadas
    x_pred_vect(:,k)     = x_pred'; %Vetor de plot;
    
end

Elapsed_time = toc;
text = ['Tempo de simula��o: ',num2str(Elapsed_time),' s'];
disp(text)

% ---- Plot de gr�ficos ----
plot_cstrNaoIsotermico
%% --------------- Fun��o do modelo ---------------
function x = modeloCSTR_naoIsotermico(estados,entradas,pert) 
    global Ts Ac dH p cp ER k0 V;
    
    Rt = (k0*exp(-ER/estados(3)));
    
    x(1) = estados(1) + Ts*((pert(1) - entradas(1))/Ac);
    x(2) = estados(2) + Ts*(pert(1)*((entradas(2)-estados(2))/V) - Rt*estados(2));
    x(3) = estados(3) + Ts*(pert(1)*((pert(2)-estados(3))/V) - (dH/p/cp)*Rt*estados(2) - entradas(3)/V);
    
end