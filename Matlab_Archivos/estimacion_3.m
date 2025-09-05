%observador con datos reales

close all
clear all

format short e
load Lamb

% load EMF_2da_celda1;EMF=flipdim(EMF,1);emf = griddedInterpolant(EMF(:,2),EMF(:,1));
load EMF_pack2_unco;EMF=flipdim(EMF,1);%emf = griddedInterpolant(EMF(:,2),EMF(:,1),'cubic','cubic');
%{
%------ carga de datos
[f d I v1 v2 v3 v4]=textread('Ciclado_pack-2_vc_2.txt', '%f %f %f %f %f %f %f');I=-I;V=v4/1e3;ho=1;ci=0.1;
%

h=1;
trm=h/ho;  %remuestreo
nt=length(I);
I=I(1:trm:nt);
V=V(1:trm:nt);

h=h/3600;
nt=length(I);
%-------------


%----  parametros  de la bateria  -------
%----------------------------------------

% %--- pack
Q=7.7;
% xo=[7.7793e-01   3.1677e-01   4.7459e-02]; %sinteticos
xo=[9.25e-02  3.8125e-02   5.125e-02]; %reales, ajuste por barrido


%--- celda
% Q=1.58;
% xo=[1.7433e-01   1.0000e-02   1.0000e-01]; %<--
% xo=[1.0000e-01   1.0000e-02   4.2813e-01];

%---
a=xo(1);
b=xo(2);
r=xo(3);
A=[1 0 ;1-exp(-h/b) exp(-h/b)];
B=[-h;(b-a)*(1-exp(-h/b))-h]/Q;


% ci=interp1(EMF(:,1),EMF(:,2),V(1),'pchip');
ci= rand(1,1);
% ci=0;
Vm=3.05; %tension minima de descarga del ensayo
% return

%-------- Ajuste del observador ------------

R2=1;
% R1 se define segun el observador usado
%  R1=1e0*[1e-3 0;0 1e-3]; %


Pb=1e0*eye(2);
Xb=[ci;ci]


%--------
%derivada de la emf
Df=(EMF(2:end,1)-EMF(1:end-1,1))./(EMF(2:end,2)-EMF(1:end-1,2));
% figure(),plot(EMF(2:end,2),Df),grid on
Df=[EMF(2:end,2),Df]; %ingresar con X
% Df=[EMF(2:end,1),Df]; %ingresar con V

figure(),plot(Df(:,1),Df(:,2)),grid on
sgm=[];
jk=1;


% rs=[1e-9,1e-3,1e-1,1]; %k
% rx=[1e-9,1e-3,1e-1,1]; %m
%
% for m=1:length(rx)
%
%     for k=1:length(rs)
%
%     R1=[rs(k) 0;0 rx(m)];


Vhat(1)=interp1(EMF(:,2),EMF(:,1),ci,'pchip');

for EKF=0:1
% EKF=0;

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
%busqueda de la derivada
 Xm(i)=interp1(EMF(:,1),EMF(:,2),V(i)+I(i)*r,'linear');% 'pchip'
%  if Xm(i)>.9;Xm(i)=.9;end; if Xm(i)<0;Xm(i)=0;end

%  df=interp1(Df(:,1),Df(:,2),Xm(i),'pchip'); %ingreso con Xm
df=interp1(Df(:,1),Df(:,2),Xb(2),'spline','extrap'); %ingreso con X

ddf(i)=df;

%-----------------------
if EKF==0;
%______________ KF con Xm:
%         % R2=1 definido  antes

            R1=[1e0 0;0 1e1];
        af=1;
%          if df<.4; af=0; end;
        alfa(i)=af;
        C=[0 1];
        K=af*Pb*C'*inv(R2+C*Pb*C');
        Xe=Xb+K*(Xm(i)-Xb(2));
        P=inv(inv(Pb)+C'*C/R2);
        Xb=A*Xe+B*I(i);
        Pb=A*P*A'+R1;

        Vh(i)=interp1(EMF(:,2),EMF(:,1),Xb(2),'spline','extrap')-I(i)*r;


end
if EKF==1

%______________ EKF con V:
%       R2=1 definido  antes
        R1=[1e-1 0;0 1e0]; %

        C=[0 df];
        af=1;
%          if df<.4; af=0; end;
        alfa(i)=af;
        K=af*Pb*C'*inv(R2+C*Pb*C');
        Vhat(i)=interp1(EMF(:,2),EMF(:,1),Xb(2),'pchip')-I(i)*r;
        Xe=Xb+K*(V(i)-Vhat(i));
        P=inv(inv(Pb)+C'*C/R2);Ps(i)=sqrt(P(1,1));Px(i)=sqrt(P(2,2));
        Xb=A*Xe+B*I(i);
        Pb=A*P*A'+R1;
% ______________
end

        if Xb(1)>1;Xb(1)=1;end; if Xb(1)<0;Xb(1)=0;end
        if Xb(2)>1;Xb(2)=1;end; if Xb(2)<0;Xb(2)=0;end

        Se(i)=Xb(1);XXe(i)=Xb(2);

%_______________ Calculo de Tiempo remanente _______________

Te=0;
j=1;while I(i+1)>0 && V(i+j)>Vm && (i+j)<nt-1 ;Te=Te+1;j=j+1;end

if  I(i)>0 && V(i)>Vm; % comienza la descarga  h*Te <=120

%     if i>1 && I(i-1)<=0;  ini=i;  end %
    if i>1 && I(i-1)<=0 && I(i)>0;  ini=i; ci= rand(1,1);Xb=[ci;ci]; end %
    TE(i)=h*Te;
%     if I(i)>0 && E(i)>Em

        xmin=interp1(EMF(:,1),EMF(:,2),Vm+I(i)*r,'pchip');
        xxmin(i)=xmin;

%Lambert________________________
        t1=(xmin-Se(i))*Q/I(i)-(b-a);
        t2=(Se(i)-XXe(i))*Q/I(i)+(b-a);

% w=lambertw(-(t2/b)*exp(t1/b));
           arg=-(t2/b)*exp(t1/b);
           [m1,m2]=min(abs(L(:,1)-arg));
           w=L(m2,2) ;
           TT(i)=w*b-t1;if  TT(i)<0;TT(i)=0;end
           tt=TT(i);
           te=TE(i);
% calculo de errores en cada descarga
%            e(jk)=(TE(i)-TT(i));
           ex(i)=(TE(i)-TT(i));
           e(jk)=(te-tt);
           jk=jk+1;
end

if i>1
if I(i)<0 && I(i-1)>0  %termino la descarga
%   ne=[ne;100*norm(e)/sqrt(length(e))/(h*(i-ini))];%[ ini i]
ne=[ne;100*norm(e)/sqrt(length(e))];
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

%  mean(ne(ne < 10))
% J(m,k)=median(ne(ne < 60))
% % J(m,k)=var(V(1:end-1)-Vhat');
% %  J(m,k)=var(V(1:end-1)-Vhat');
%     end
% end
% %
[minJ, idx] = min(J(:));  % minJ es el valor mínimo y idx es el índice lineal
[k, m] = ind2sub(size(J), idx);
% return
%%
figure(9),plot(I)
figure(10),hold on
% plot([D(:,3) D(:,4)],'linewidth',2),
% plot([Se' XXe'],'--','linewidth',2),
%
% legend('SoC real','X real','SoC estimado','X estimado')
figure(31),
plot([V(1:end-1) Vhat' Vh']),legend('Real','EKF','KF')


figure(41),hold on
plot([TE' TT'],'--'),title('tiempo remanente');grid on
legend('Real','Prediccion')


figure(),hist(ex,50),title('Histograma de todos los Errores');
figure(),plot(ne,'.k','markersize',20),title('RMSE relativo del Error [%]');

return
%%


figure(22),plot([Se' XXe' xxmin']),title('Estimacion EKF');
legend('SoC','X')
%}






