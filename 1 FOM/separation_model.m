%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PART Ⅰ
%In this part, you need to set the parameters before calculating.

%initialization
clear,clc
p0=10310;                                                                  %initial pressure[kPa]
T0=29.1+273.15;                                                              %initial fluid temperature[K]
d0=refpropm('D','T',T0,'P',p0,'CO2');                                      %initial fluid density[kg/m^3]
h0=refpropm('H','T',T0,'P',p0,'CO2');                                      %initla fluid enthalpy[J/kg]
TW=T0+0.1;                                                                 %initial wall temperature[K]
t=0; i=1;                                                                  %initial moment t[s]

%set container parameters
d=0.01;                                                                    %wall thickness[m]
R=110;                                                                     %radius[mm]
As=pi*(R*0.001)^2;                                                         %cross-sectional area[m^2]
V=0.05;                                                                    %volume[m^3]
cpw=460;                                                                   %specific heat capacity[J/(kg K)]
denw=7850;                                                                 %density[kg/m^3]

%set nozzle parameters
r=0.5;                                                                     %radius[mm]
as=pi*(r*0.001)^2;                                                         %area[m^2]
x0=0.4378;                                                              %relative height of nozzle center[#]
x_up=x0+r/2/R;                                                             %maximum relative height of the nozzle[#]
x_down=x0-r/2/R;

%set caculation parameters
T=1200;                                                                    %total time[s]
dt=0.1;                                                                    %time step size[s]
dp=0.005;                                                                  %pressure step size[kPa]
ep=0.01;                                                                   %convergence criterion of pressure[#]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PART Ⅱ
%calculation methods

%calculate the mass flow rate at t=0
if (p0>7377.3)&&(T0<304.1)                                                 %supercritical region
  ps=refpropm('P','T',T0,'Q',0,'CO2');
  v1=refpropm('I','T',T0,'Q',0,'CO2');
  v2=refpropm('I','P',1380,'Q',0,'CO2');
  pcc=ps*(1-0.284*v1/v2);
  W=as*(2*d0*((p0-pcc)*1000))^0.5;
else
  pr=p0/7377.3;
  dr=d0/354.36;
  eta=0.579+0.024*log(LD);
  W=as*(0.5463+0.0587*pr^2.07*dr^(-0.939))*(2*d0*p0*1000*(1-eta))^0.5;
end

%store initial results
data(i,1)=t; data(i,2)=p0/1000; data(i,3)=d0; data(i,4)=T0; data(i,5)=W; data(i,6)=1; data(i,7)=TW;

%dense-phase region
%critical point:(7377.3kPa,304.1K)
while p0>7377.3
     if T0<304.1
       ps=refpropm('P','T',T0,'Q',0,'CO2');
       v1=refpropm('I','T',T0,'Q',0,'CO2');
       v2=refpropm('I','P',1380,'Q',0,'CO2');
       pcc=ps*(1-0.284*v1/v2);
       W=as*(2*d0*((p0-pcc)*1000))^0.5;     
     else
       pr=p0/7377.3;
       dr=d0/354.36;
       eta=0.579+0.024*log(LD);
       W=as*(0.5463+0.0587*pr^2.07*dr^(-0.939))*(2*d0*p0*1000*(1-eta))^0.5;
     end
     u=W/(d0*as);                                                          %fluid velocity at the nozzle[m/s]
     ub=u*(as/As);                                                         %fluid velocity in the nozzle[m/s]
     cp=refpropm('C','P',p0,'D',d0,'CO2');                                 %fluid specific heat capacity[J/(kg K)]
     DV=refpropm('V','P',p0,'D',d0,'CO2');                                 %dynamic viscosity[Pa s]
     TC=refpropm('L','P',p0,'D',d0,'CO2');                                 %fluid thermal conductivity[W/(m K)]
     Pr=refpropm('^','P',p0,'D',d0,'CO2');                                 %the Prandtl number[#]                                              
     Re=d0*ub*(2*R*0.001)/DV;                                              %the Reynolds number[#]
     Nu=Re^0.35*Pr^1.9*((d0/denw)^-1.6)*(cp/cpw)^-3.4;                     %the Nussel number[#]
     H=Nu*TC/(2*R*0.001);                                                  %heat transfer coefficient[W/(m^2 K]
     QW=(V/As)*2*pi*R*H*(TW-T0);                                           %heat transfer[J]
     if T0>TW
       QW=0;
     end
     d0=d0+dt*(-W/V);                                                      %new fluid density[kg/m^3]
     TW=TW-QW/(denw*cpw*d);                                                %new wall temperature[K]
     h0=h0+QW/(d0*V);                                                      %new enthalpy[J/kg], he=h0   
     p0=refpropm('P','D',d0,'H',h0,'CO2');                                 %new pressure[kPa]
     T0=refpropm('T','D',d0,'H',h0,'CO2');                                 %new fluid temperature[K]
     i=i+1; t=t+dt,p0/1000                                                 %partial results monitoring
     if (t>T)||(p0<=101)                                                   %criterion for cessation
       break
     else
       %store calculation results
       data(i,1)=t; data(i,2)=p0/1000; data(i,3)=d0; data(i,4)=T0; data(i,5)=W; data(i,6)=1; data(i,7)=TW;
     end
end

%liquid
hl=refpropm('H','T',T0,'Q',0,'CO2');                                       %enthalpy of saturated liquid at pressure p0[J/kg]
while h0<hl
     %W=as*0.61*(2*d0*(1-0.7)*p0*1000)^0.5;                                 %L/D>20>7,therefore eta=0.7
     ps=refpropm('P','T',T0,'Q',0,'CO2');
       v1=refpropm('I','T',T0,'Q',0,'CO2');
       v2=refpropm('I','P',1380,'Q',0,'CO2');
       pcc=ps*(1-0.284*v1/v2);
       W=as*(2*d0*((p0-pcc)*1000))^0.5;     
     u=W/(d0*as);
     ub=u*(as/As);
     cp=refpropm('C','P',p0,'D',d0,'CO2');
     DV=refpropm('V','P',p0,'D',d0,'CO2');
     TC=refpropm('L','P',p0,'D',d0,'CO2');
     Pr=refpropm('^','P',p0,'D',d0,'CO2');
     Re=d0*ub*(2*R*0.001)/DV;
     Nu=Re^0.35*Pr^1.9*((d0/denw)^-1.6)*(cp/cpw)^-3.4;
     H=Nu*TC/(2*R*0.001);
     QW=(V/As)*2*pi*R*H*(TW-T0);
     if T0>TW
       QW=0;
     end
     d0=d0+dt*(-W/V);
     TW=TW-QW/(denw*cpw*d);
     h0=h0+QW/(d0*V);
     p0=refpropm('P','D',d0,'H',h0,'CO2');
     T0=refpropm('T','D',d0,'H',h0,'CO2');
     hl=refpropm('H','T',T0,'Q',0,'CO2');
     i=i+1; t=t+dt,p0/1000
     if (t>T)||(p0<=101)  
       break
     else
       data(i,1)=t; data(i,2)=p0/1000; data(i,3)=d0; data(i,4)=T0; data(i,5)=W; data(i,6)=1; data(i,7)=TW;
     end
end

%two-phase region
x1=1;                                                                      %initial relative level height[#]
x2=0;                                                                      %initial mass gas content[#]
hg=refpropm('H','P',p0,'Q',1,'CO2');                                       %enthalpy of saturated gas at pressure p0[J/kg]
d_guess=100000;                                                            %initial iteration pressure[kPa]
while (hl<=h0)&&(h0<=hg)  
     hl=refpropm('H','P',p0,'Q',0,'CO2'); hg=refpropm('H','P',p0,'Q',1,'CO2');
     dl=refpropm('D','P',p0,'Q',0,'CO2'); dg=refpropm('D','P',p0,'Q',1,'CO2');
     cpl=refpropm('C','P',p0,'Q',0,'CO2'); cpg=refpropm('C','P',p0,'Q',1,'CO2');
     Prl=refpropm('^','P',p0,'q',0,'CO2'); Prg=refpropm('^','P',p0,'q',1,'CO2');
     TCl=refpropm('L','P',p0,'Q',0,'CO2'); TCg=refpropm('L','P',p0,'Q',1,'CO2');
     %two-phase physical property parameters
     cp=(1-x2)*cpl+x2*cpg;
     Pr=(1-x2)*Prl+x2*Prg;
     TC=(1-x2)*TCl+x2*TCg;
     if x1>x_up                                                            %the liquid level is above the nozzle
       W=as*0.61*(2*dl*(1-0.7)*p0*1000)^0.5;                               %LD>20>7
       he=hl;
       u=W/(dl*as);
     elseif (x_down<x1)&&(x1<x_up)                                         %initial mass gas content[#] at the nozzle
       if x1>x0
         area_g=(r/2/R)^2*acos((x1-x0)/(r/2/R))-(x1-x0)*((r/2/R)^2-(x1-x0)^2)^0.5;
         area_l=pi*(r/2/R)^2-area_g;
         ratio_m=dg*area_g/(dg*area_g+dl*area_l);
       elseif x1<x0
         area_l=(r/2/R)^2*acos((x0-x1)/(r/2/R))-(x0-x1)*((r/2/R)^2-(x0-x1)^2)^0.5;
         area_g=pi*(r/2/R)^2-area_l;
         ratio_m=dg*area_g/(dg*area_g+dl*area_l);
       else
         ratio_m=dg/(dg+dl);
       end
       ratio_v=area_g/(area_g+area_l);
       d_lg=dl*(1-ratio_v)+dg*ratio_v;
       W=ratio_m*as*0.81*(2*dg*(1-0.7)*p0*1000)^0.5+(1-ratio_m)*as*0.61*(2*dl*(1-0.7)*p0*1000)^0.5;
       he=ratio_m*hg+(1-ratio_m)*hl;
       u=W/(d_lg*as);
     else                                                                  %the liquid level is under the nozzle
       W=as*0.81*(2*dg*(1-0.7)*p0*1000)^0.5;
       he=hg;
       u=W/(dg*as);
     end 
     ub=u*(as/As);
     Re=d0*ub*(2*R*0.001)/DV;
     Nu=Re^0.35*Pr^1.9*((d0/denw)^-1.6)*(cp/cpw)^-3.4;
     H=Nu*TC/(2*R*0.001);
     QW=(V/As)*2*pi*R*H*(TW-T0);
     if T0>TW
       QW=0;
     end
     d0=d0+dt*(-W/V);
     TW=TW-QW/(denw*cpw*d);
     h0=h0+(QW+(h0-he)*W*dt)/(d0*V);
     %pressure iterations
     while abs(1-d_guess/d0)>ep
          p0=p0-dp;
          dl=refpropm('D','P',p0,'Q',0,'CO2'); dg=refpropm('D','P',p0,'Q',1,'CO2');
          hl=refpropm('H','P',p0,'Q',0,'CO2'); hg=refpropm('H','P',p0,'Q',1,'CO2');
          x2=(h0-hl)/(hg-hl);
          vl=1/dl;vg=1/dg;
          v_guess=x2*vg+(1-x2)*vl;
          d_guess=1/v_guess;
     end
     T0=refpropm('T','P',p0,'Q',0,'CO2');
     x1=1-x2*d0/dg;                                                        %new relative height of liquid level[#]
     i=i+1; t=t+dt,p0/1000
     if (t>T)||(p0<=101)
       break
     else
       data(i,1)=t; data(i,2)=p0/1000; data(i,3)=d0; data(i,4)=T0; data(i,5)=W; data(i,6)=x1; data(i,7)=TW;
     end
end

%gas
while h0>hg
     W=as*0.81*(2*d0*(1-0.7)*p0*1000)^0.5;
     u=W/(d0*as);
     ub=u*(as/As);
     cp=refpropm('C','P',p0,'D',d0,'CO2');
     DV=refpropm('V','P',p0,'D',d0,'CO2');
     TC=refpropm('L','P',p0,'D',d0,'CO2');
     Pr=refpropm('^','P',p0,'D',d0,'CO2');
     Re=d0*ub*(2*R*0.001)/DV;
     Nu=Re^0.35*Pr^1.9*((d0/denw)^-1.6)*(cp/cpw)^-3.4;
     H=Nu*TC/(2*R*0.001);
     QW=(V/As)*2*pi*R*H*(TW-T0);
     if T0>TW
       QW=0;
     end
     d0=d0+dt*(-W/V);
     TW=TW-QW/(denw*cpw*d);
     h0=h0+QW/(d0*V);
     p0=refpropm('P','D',d0,'H',h0,'CO2');
     T0=refpropm('T','D',d0,'H',h0,'CO2');
     hg=refpropm('H','P',p0,'Q',1,'CO2');
     i=i+1; t=t+dt,p0/1000
     if (t>T)||(p0<=101)
       break
     else
       data(i,1)=t; data(i,2)=p0/1000; data(i,3)=d0; data(i,4)=T0; data(i,5)=W; data(i,6)=0; data(i,7)=TW;
     end
end