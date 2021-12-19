    close all
    clear all
    s=tf('s')
    %Planta
    Cbulk=88372.83; 
    Csurface=82.11;
    Re=0.00375;
    Rs=0.00375;
    Rt=0.002745;

    Ad = -(Rs)/(Cbulk*(Re+Rs)^2)+(Re)/(Csurface*(Re+Rs)^2)-(Rs^2)/(Cbulk*Re*(Re+Rs)^2)+(Rs)/(Csurface*(Re+Rs)^2)
    Add = (Rs)/(Cbulk*Re*(Re+Rs))-1/(Csurface*(Re+Rs))
    Bd = (Re^2)/(Csurface*(Re+Rs)^2)-(Rs*Rt)/(Cbulk*Re*(Re+Rs))+(Rt)/(Csurface*(Re+Rs))+(Re*Rs)/(Csurface*(Re+Rs)^2)

    A=[-1/(Cbulk*(Re+Rs)) 1/(Cbulk*(Re+Rs)) 0 ; 1/(Csurface*(Re+Rs)) -1/(Csurface*(Re+Rs))  0 ; Ad 0 Add]
    B=[(Rs)/(Cbulk*(Re+Rs)) (Re)/(Csurface*(Re+Rs)) Bd ]'
    C=[0 0 1]
    D=0

    sys=ss(A,B,C,D);

    Ts=1
    Tss=Ts/1000;
    [Gd,Hd,Cd,Dd] = ssdata(c2d(ss(A,B,C,D),Ts)); %planta na frequência digital
    [Gd2,Hd2,Cd2,Dd2,Tss] = ssdata(c2d(ss(A,B,C,D),Tss)); %planta "contínua"
    sysd=ss(Gd,Hd,Cd,Dd);
    %Definições para o filtro de kalman
    Q=[0.001 0 0; 0 0.01 0; 0  0 10]; %MATRIZ Q
    R=10;%MATRIZ R
    P(:,:,2)=100*eye(3);% MATRIZ P
    %-------------------------------------------------------------
    hours = 3
    t_final =  3600*hours;            % Tempo total de simulação(segundos)
    i  = 0;                 % Contador do tempo total (clk)
    t  = 0:Ts:t_final;      % Tempo de simulação (considerando as interrupções)
    tt = 0:Tss:t_final;     % Tempo de simulação (considerando o CLOCK)%
    Ref=1.53*ones(1,t_final+1);    
    kp=2;
    ii=0; 
    k=1;
    x1=[1.94 zeros(1,length(tt)-1)];
    x2=[1.94 zeros(1,length(tt)-1)];
    x3=[1.94 zeros(1,length(tt)-1)];
    x_e1=zeros(1,length(t));
    x_e2=zeros(1,length(t));
    x_e3=zeros(1,length(t));
    u=zeros(1,length(t));
    y=zeros(1,length(t));
    K_obs=zeros(3,1,length(t));
    xe1_kalman=zeros(1,length(t));
    xe2_kalman=zeros(1,length(t));
    xe3_kalman=zeros(1,length(t));
    xe1_kalman1=zeros(1,length(t));
    xe2_kalman1=zeros(1,length(t));
    xe3_kalman1=zeros(1,length(t));
    Y = 0.01.*randn(1,length(t)); %Vetor contendo ruído branco
    Soc=[ 0 0 zeros(1,length(t)-1)]; %Vetor contendo o valor de SOC em %
    for i=1:length(tt)
        ii=ii+1;
        if ii==1000
           ii=0;
           k=k+1;
           y(k)=yd(i-1)+Y(k);
           y_real(k)=yd(i-1);
           %%%%%%kalman%%%%%%%%%%%%%%%%
           %predição
           P1(:,:,k)=Gd*P(:,:,k)*Gd'+Q; %atualização de P
           xe1_kalman(k)=Gd(1,1)*xe1_kalman1(k)+Gd(1,2)*xe2_kalman1(k)+Gd(1,3)*xe3_kalman1(k)+Hd(1,1)*u(k-1);
           xe2_kalman(k)=Gd(2,1)*xe1_kalman1(k)+Gd(2,2)*xe2_kalman1(k)+Gd(2,3)*xe3_kalman1(k)+Hd(2,1)*u(k-1);
           xe3_kalman(k)=Gd(3,1)*xe1_kalman1(k)+Gd(3,2)*xe2_kalman1(k)+Gd(3,3)*xe3_kalman1(k)+Hd(3,1)*u(k-1);
           S_k=Cd*P1(:,:,k)*Cd'+R;
           K_obs(:,:,k)=(P1(:,:,k)*Cd'*inv(S_k));%calculo do ganho;
           %%%%%%%atualização
           xe1_kalman1(k+1)=xe1_kalman(k)+K_obs(1,1,k)*(y(k)-xe3_kalman(k));
           xe2_kalman1(k+1)=xe2_kalman(k)+K_obs(2,1,k)*(y(k)-xe3_kalman(k));
           xe3_kalman1(k+1)=xe3_kalman(k)+K_obs(3,1,k)*(y(k)-xe3_kalman(k));
           P(:,:,k+1)=P1(:,:,k)- K_obs(:,:,k)*S_k*K_obs(:,:,k)';
           %%%%%%%%%%%%%% Ação de controle
           u(k)=Ref(k);           
           Soc_estimado(k)=1-(2.14-xe1_kalman(k))*4.73;           
           Soc_real(k)=1-(2.14-y_real(k))*4.73;           
           if Soc_estimado(k)<0
               Soc_estimado(k)=0;
           end
           if Soc_real(k)<0
               Soc_Real(k)=0;
           end
           Soc(k)=1-(2.14-y(k))*4.73;
           erro(k)= Soc(k)- Soc_estimado(k);

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          end

      x1(i+1)=Gd2(1,1)*x1(i)+Gd2(1,2)*x2(i)+Gd2(1,3)*x3(i)+Hd2(1,1)*u(k);
      x2(i+1)=Gd2(2,1)*x1(i)+Gd2(2,2)*x2(i)+Gd2(2,3)*x3(i)+Hd2(2,1)*u(k);
      x3(i+1)=Gd2(3,1)*x1(i)+Gd2(3,2)*x2(i)+Gd2(3,3)*x3(i)+Hd2(3,1)*u(k);
      yd(i)=Cd2(1,1)*x1(i)+Cd2(1,2)*x2(i)+Cd2(1,3)*x3(i);

    end

    figure
    xlabel('Tempo (min)')
    ylabel('SOC')    
    plot(t/60,Soc(1:length(t)),'DisplayName','SOC calculado com ruído')        
    title('SOC estimado vs SOC ruidoso')
    hold on
    plot(t/60,Soc_estimado(1:length(t)),'r','DisplayName','SOC estimado') 
    figure
    plot(t/60,y_real(1:length(t)),'r','DisplayName','SOC real')
    xlabel('Tempo (min)')
    ylabel('Tensão (V)')
    title('Saída real da bateria') 
    xlim([0 190])
    ylim([1.94 2.15])
