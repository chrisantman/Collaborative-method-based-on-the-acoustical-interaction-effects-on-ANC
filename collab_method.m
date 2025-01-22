% Distribuido FxLMS incremental con matriz de acoples y por subredes
% versión caso concreto
clear;
clc;

%%
% Añadir carpetas
rootfolder = cd;
addpath('..\Datos\Señales');
addpath('..\Datos\Caminos');
% cargar sala
load ir_map;
sala=ir_map/(max(max(max(max(ir_map)))));  
fs=2000;
ITE=100000; 
ite=0;

% % red 2 nodos
% f1=250;
% f2=734.375;

% % red 6 nodos
% f1=867.1875;
% f2=507.8125;

% % ruido añadido
% noise=randn(2,1)';
% noise=noise-ones(1,2)*mean(noise);

% for f=f2
for f=0:7.8125:1000-1

% señal de referencia
t=0:1/fs:ITE*(1/fs); 
IN=sin(2*pi*f*t(1:end-1))';

% load x3;
% in=x3(1:ITE,1);
% IN=in-ones(ITE,1)*mean(in);

% red 2 nodos caso 11
% ref=69; 
% alt=[19,22];
% mic=[62,59];

% red 6 nodos
ref=69;
alt=[11,12,16,20,22,24];
mic=[69,68,64,60,58,56];

I=length(ref);   
J=length(alt);   
K=length(mic); 

% Nodos
nodos=[alt;
       mic];
num_nodos=length(nodos(1,:)); 

% Cálculo de caminos acústicos
M=length(sala(1,:,1,1));   
C_aux=sala(ref(1:I),:,mic(1:K));
C=permute(C_aux,[2 3 1]);
if I==1
 CPRI=C;
else
 CPRI=reshape(C,M,I*J);
end 
C_aux=sala(alt(1:J),:,mic(1:K));
C=permute(C_aux,[2 3 1]);
CSEC=reshape(C,M,J*K);

% Caminos para cada nodos
for n=1:num_nodos
    c=sala(alt(1:end),:,nodos(2,n))';
    csec(:,:,n)=c;
end

% Cálculo de señal deseada
d=zeros(ITE,K);
for k=1:K
    d(:,k)=filter(CPRI(:,k),1,IN);
end

%% Acople
%  [Erx,Etx,MCov,delay,m_retardo]=acople_vC(CSEC,J,K);
% [ M_acop1, M_acop2, M_acop3, MC,BEnermic,BEneralt ] = acople_vM(CSEC,J,K);
 % se pone 1 si el nodo de la fila está acoplado con el nodo de la columna.
 % 0 si no está acoplado. Ejemplo, como el nodo 3 (fila 3) no está acoplado al nodo
 % 1 (columna 1), se pone 0.
% % nodos separados

% for a=1:2
    
%% Configuración del algoritmo local que se ejecuta en cada nodo
L=150;      % Número de coeficientes de los filtros adaptativos
% mu=0.01;    % Paso de convergencia del algoritmo FxLMS
% mu=0.0005; % tono
mu=0.00005; % tono
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
w=zeros(L*J,num_nodos);	
buff_L=zeros(L,num_nodos);
error=zeros(ITE,num_nodos); 
e=zeros(1,num_nodos);	
buff_xf=zeros(M,num_nodos);
v=zeros(L*J,num_nodos); 
buff_y=zeros(M,J);

%% MASCARA Colaboración

%% red 2 nodos
% % Centralizado
% masc=ones(2,2); 

% %  No Colaborativo
% mu=mu/num_nodos;
% masc=eye(2);

% % Colaborativo
% mu=mu/num_nodos;
% masc=[1 1 ;
%       1 1 ];

%% red 6 nodos
% % Centralizado
% masc=ones(6,6); 

% %  No Colaborativo
% mu=mu/num_nodos;
% masc=eye(6);

% % Colaborativo
mu=mu/num_nodos;
masc=[1 1 1 0 1 0;
      1 1 1 0 0 0;
      1 1 1 1 1 0;
      0 0 1 1 1 1;
      1 0 1 1 1 1;
      0 0 0 1 1 1];
  
% masc=[1 1 1 0 0 0;
%       1 1 1 0 0 0;
%       1 1 1 1 1 0;
%       0 0 1 1 1 1;
%       0 0 1 1 1 1;
%       0 0 0 0 1 1];

      
%%

tic;
for cont=1:ITE

   y=zeros(1,J); % Variable para la señal a generar por cada una de las fuentes en cada iteración. 
   x=IN(cont,I);
   
   for nodo=1:num_nodos

        % ESTRATEGIA INCREMENTAL
        w_nodo=w(1+L*(nodo-1):L*nodo,nodo);
        
       %%
        if nodo==1
            w_ant=w(:,num_nodos);
        else
            w_ant=w(:,nodo-1);
        end
        
        %%
        % FILTRADO ADAPTATIVO
        buff_L(2:L,nodo)=buff_L(1:L-1,nodo); 
        buff_L(1,nodo)=x;	
        y(nodo)=y(nodo)+w_nodo'*buff_L(:,nodo); 

        % ACTUALIZACIÓN COEFICIENTES       
        w(:,nodo)=(w_ant-mu*v(:,nodo)*e(:,nodo));

        % FILTRADO ESTIMA 
        buff_xf(2:M,nodo)=buff_xf(1:M-1,nodo); 
        buff_xf(1,nodo)=x;
        vijk=(buff_xf(:,nodo)'*csec(:,:,nodo)).*masc(:,nodo)';
        vaux=reshape(v(:,nodo),L,K);
        vaux(2:L,:)=vaux(1:L-1,:);
        vaux(1,:)=vijk; 
        v(:,nodo)=vaux(:);

    end

    w=w(:,end)*ones(1,num_nodos);
  
    %%%%% PARTE ACÚSTICA %%%%%%%%%%%%%%%%%%%%%%%
    buff_y(2:M,:)=buff_y(1:M-1,:);
    buff_y(1,:)=y;
    yf=zeros(1,K);
     for k=1:K
        for j=1:J
          yf(k)=yf(k)+CSEC(:,k+K*(j-1))'*buff_y(:,j);
        end
     end

    e=d(cont,:)+yf;% +noise
    error(cont,:)=e;
    
end
%%
t=toc;

dpot=d.^2;
epot=error.^2;
dpp(1,:)=dpot(1,:).^2;
epp(1,:)=epot(1,:).^2;

for i=2:max(ITE)
   dpp(i,:)=0.9995*dpp(i-1,:)+(1-0.9995)*dpot(i,:).^2;
   epp(i,:)=0.9995*epp(i-1,:)+(1-0.9995)*epot(i,:).^2;
end
ATEp=epp./dpp;

% GRAFICA
ATEp=sum(ATEp,2)./K;
ATE_db=10*log10(ATEp);
t=1/fs*(1:length(ATEp))';
% % figure()
hold all;
plot(t,ATE_db);
axis([0,t(end),-120,10]);
% legend( 'Distributed FxLMS Node 1',...
%         'Distributed FxLMS Node 2');
xlabel('Time (s)')
ylabel('Noise reduction (dB)');

% % figure()
% % plot(t,ATE_db);hold all;
% plot(t,(sum(ATE_db,2)./K));hold all; %% MAL
% axis([0,t(end),-45,5]);
% legend( 'Centralized FxLMS',...
%         'Distributed NC FxLMS');

ite=ite+1;
MSE(ite,:)=ATEp(96000,:);

end

% save MSE_Colab_150coef.mat MSE


