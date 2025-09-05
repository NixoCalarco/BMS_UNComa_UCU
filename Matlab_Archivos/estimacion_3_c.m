%observador con datos reales 

close all
clear all

%cargar tablas 
load Lamb % carga la curva de lambert en forma de tabla
load EMF_pack2_unco;EMF=flipdim(EMF,1); % cargar la curva de la EMF como tabla
emf = griddedInterpolant(EMF(:,2),EMF(:,1),'cubic','cubic'); % genero la grilla para interpolar

%------ carga de datos
[f d I v1 v2 v3 v4]=textread('Ciclado_pack-2_vc_2.txt', '%f %f %f %f %f %f %f');

I=-I; % corrige el signo de la corriente de acuerdo al modelo
V=v4/1e3; %tensión en volts
ho=1;   % muestreo original: 1 sgeundo
ci=0.1; % condicion inicial de los estados de la batería. Este dato sale de V(1) y la tabla EMF inversa. 
% 

h=1;    % muestreo deseado
trm=h/ho;  %remuestreo segun el muestreo deseado 
nt=length(I);
I=I(1:trm:nt);
V=V(1:trm:nt);

h=h/3600;   %periodo de muestreo en horas
nt=length(I);
%-------------


%----  parametros  de la bateria  -------
%----------------------------------------

% %--- pack

Q=7.7; % capacidad en [Ah]
xo=[9.25e-02  3.8125e-02   5.125e-02]; %reales, ajuste por barrido


%---
a=xo(1); %inversa del cero del filtro RCE
b=xo(2); %inversa del  polo del filtro RCE
r=xo(3); % resistencia interna 
A=[1 0 ;1-exp(-h/b) exp(-h/b)];
B=[-h;(b-a)*(1-exp(-h/b))-h]/Q;


Vm=3.05; %tension minima de descarga del ensayo
% return

%-------- Ajuste del observador ------------

R2=1; % R2 se setea en 1 para hacer todo relativo a este valor
% R1 se define segun el observador usado
 
        
Pb=1e0*eye(2); % Covarianza de la perturbacion de los estados 
Xb=[ci;ci]; %setear los estados iniciales a la condicion inicial 


%--------
%Cálculo de  derivada de la emf usada en el filtro de kalman extendido(EKF) 

Df=(EMF(2:end,1)-EMF(1:end-1,1))./(EMF(2:end,2)-EMF(1:end-1,2)); %derivada como tabla

Df=[EMF(2:end,2),Df]; %ingresar con X
% Df=[EMF(2:end,1),Df]; %ingresar con V

figure(),plot(Df(:,1),Df(:,2)),grid on
sgm=[];
jk=1;



Vhat(1)=interp1(EMF(:,2),EMF(:,1),ci,'pchip'); %Vhat: estimacion de la salida. La primer estimación se hace con la condicion inicial

for EKF=0:1
% EKF=0; % estimacion con KF (lineal)
% EKF=0; % estimacion con EKF (no-lineal) 

xxmin=NaN(1,nt-1);
TT=NaN(1,nt-1);
TE=NaN(1,nt-1);
ne=[];  
Se=0;
XXe=0;
    
for i=2:nt-1
   %modelo 
%   Z=A*Z+B*I(i,1)+randn(size(B))*0;
% %   V(i)=emf(Z(2))-I(i)*r;
%  Vr(i,1)=interp1(EMF(:,2),EMF(:,1),Z(2),'pchip')-I(i,1)*r;
 
 
%______ Observador:  

 Xm(i)=interp1(EMF(:,1),EMF(:,2),V(i)+I(i)*r,'linear');% busqueda de la derivada
 
%  df=interp1(Df(:,1),Df(:,2),Xm(i),'pchip'); %ingreso con Xm
df=interp1(Df(:,1),Df(:,2),Xb(2),'spline'); %ingreso con X

ddf(i)=df;

%-----------------------
if EKF==0;
%______________ KF con Xm:   
%      

        R1=[1e0 0;0 1e1];       
        C=[0 1]; 
        K=Pb*C'*inv(R2+C*Pb*C');
        Xe=Xb+K*(Xm(i)-Xb(2));
        P=inv(inv(Pb)+C'*C/R2);
        Xb=A*Xe+B*I(i); % estimacion de los estados
        Pb=A*P*A'+R1;
    
        Vh(i)=interp1(EMF(:,2),EMF(:,1),Xb(2),'spline','extrap')-I(i)*r; % salida estimada 

        
end
if EKF==1
    
%______________ EKF con V:  

        R1=[1e-1 0;0 1e0]; % 

        C=[0 df];
        K=Pb*C'*inv(R2+C*Pb*C');
        Vhat(i)=interp1(EMF(:,2),EMF(:,1),Xb(2),'pchip')-I(i)*r;% salida estimada 
        Xe=Xb+K*(V(i)-Vhat(i));
        P=inv(inv(Pb)+C'*C/R2);Ps(i)=sqrt(P(1,1));Px(i)=sqrt(P(2,2));
        Xb=A*Xe+B*I(i); % estimacion de los estados
        Pb=A*P*A'+R1;
% ______________ 
end
 
        if Xb(1)>1;Xb(1)=1;end; if Xb(1)<0;Xb(1)=0;end %saturacion de SoC a valores permitidos
        if Xb(2)>1;Xb(2)=1;end; if Xb(2)<0;Xb(2)=0;end %saturacion de X a valores permitidos
 
        Se(i)=Xb(1); %guardo el SoC estimado        
        XXe(i)=Xb(2); %guardo X estimada     
        
%_______________ Calculo de Tiempo remanente _______________
  
Te=0; 
j=1;while I(i+1)>0 && V(i+j)>Vm && (i+j)<nt-1 ;Te=Te+1;j=j+1;end %tiempo remanente real

if  I(i)>0 && V(i)>Vm; % si pasa esto es porque comienzó la descarga. El TR solo se calcula si se está descargando 
        
    TE(i)=h*Te; % TE = tiempo remanente real en horas

%     if i>1 && I(i-1)<=0 && I(i)>0;  ini=i; ci= rand(1,1);Xb=[ci;ci]; end %Para iniciar cada descarga desde un valor random
%     
        
        
        xmin=interp1(EMF(:,1),EMF(:,2),Vm+I(i)*r,'pchip'); %X minimo para ec de lambert
        xxmin(i)=xmin; % guardo X minimo

%Lambert comienzo________________________  
        t1=(xmin-Se(i))*Q/I(i)-(b-a);  %tita 1  
        t2=(Se(i)-XXe(i))*Q/I(i)+(b-a);% tita 2
 
% w=lambertw(-(t2/b)*exp(t1/b)); % funcion de lambert del toolbox matlab
           arg=-(t2/b)*exp(t1/b); %argumento de la funcion de lambert
           [m1,m2]=min(abs(L(:,1)-arg)); % busqueda de valor mas cercano en la tabla de lambert precargada 
           w=L(m2,2) ; % busqueda en tabla 
           TT(i)=w*b-t1;if  TT(i)<0;TT(i)=0;end % tiempo remanente calculado con lambert + saturacion si es negativo
%Lambert fin________________________ 

         
           e(jk)=(TE(i)-TT(i));  % error de tiempo remanente en cada paso "i" de la dscarga        
           jk=jk+1;  
end  

if i>1 
if I(i)<0 && I(i-1)>0   %fin de la descarga
ne=[ne;100*norm(e)/sqrt(length(e))]; % al finalizar la desacarga se calcula la norma del error en toda la descarga
  e=0;
  jk=1;
end   
end


end
J(EKF+1)=median(ne(ne < 60))
figure(),hold on
plot([TE' TT'],'--'),title(['EKF = ' num2str(EKF)])
ylabel('tiempo remanente')
legend('Real','Prediccion');grid on
end


[minJ, idx] = min(J(:));  % minJ es el valor mínimo y idx es el índice lineal
[k, m] = ind2sub(size(J), idx);
% return
%% Graficos 
figure(9),plot(I)
title('Registro de corriente del ensayo');
figure(31),
plot([V(1:end-1) Vhat' Vh']),legend('Real','EKF','KF')
title('Tensiones estimadas');

figure(41),hold on
plot([TE' TT'],'--'),title('tiempo remanente');grid on
legend('Real','Prediccion')

% histograma de errores
% figure(),hist(ex,50),title('Histograma de todos los Errores');
% figure(),plot(ne,'.k','markersize',20),title('RMSE relativo del Error [%]');








