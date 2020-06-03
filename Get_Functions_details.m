
% This function containts full information and implementations of the benchmark 
% functions in Table 1, Table 2, and other test functins from the literature 

% lb is the lower bound: lb=[lb_1,lb_2,...,lb_d]
% up is the uppper bound: ub=[ub_1,ub_2,...,ub_d]
% dim is the number of variables (dimension of the problem)

function [lb,ub,dim,fobj] = Get_Functions_details(F)


switch F
    case 1
        fobj = @F1;
        lb=[0 0 0];
        ub=[1 10 20];
        dim=3;
        
    case 2
        fobj = @F2;
        lb=[0 0 0];
        ub=[1 10 20];
        dim=3;
        
    case 3
        fobj = @F3;
         lb=[0 0 0];
        ub=[1 10 20];
        dim=3;
        
    case 4
        fobj = @F4;
         lb=[0 0 0];
        ub=[1 10 20];
        dim=3;
        
    case 5
        fobj = @F5;
        lb=[0 0 0];
        ub=[1 10 20];
        dim=3;
        
    case 6
        fobj = @F6;
        lb=[0 0 0];
        ub=[1 10 20];
        dim=3;
        
    case 7
        fobj = @F7;
       lb=[0 0 0];
        ub=[1 10 20];
        dim=3;
        
    case 8
        fobj = @F8;
        lb=[0 0 0];
        ub=[1 10 20];
        dim=3;
        
                
end

end

% F1

function fobj = F1(filter_data)
format long g;
global VTHD_percent ITHD_percent HP Efficiency PF DPF IS VL PLOSS VLH ISH IS1 VL1;
VSH_f = 2400;
ILH_f = 987.83;
RL1 = 1.7421;
XL1 = 1.696;
RS1 = 0.02163;
XS1 = 0.2163;
h = [5 7 11 13];
VSH_percent = [0 0 0 0];
ILH_percent = [0.4 0.06 0.02 0.01];

%filter_data = [122.14/1000 2.7279 11];
XL = filter_data(1);
XC = filter_data(2);
K = filter_data(3);

VSH = VSH_f*VSH_percent;
RLH = RL1*ones(1,4);
ILH = ILH_f*ILH_percent;
XLH = XL1*h;
RSH = RS1*ones(1,4);
XSH = XS1*h;

RTLH = RLH.*RSH-XLH.*XSH;
XTLH = RLH.*XSH+RSH.*XLH;

A = VSH.*RLH - ILH.*XLH.*(h*XL-XC./h);
B = VSH.*(XLH+h*XL-XC./h) + ILH.*RLH.*(h*XL-XC./h);
C = RTLH + K*RLH - (XLH+XSH).*(h*XL-XC./h);
D = XTLH + K*XLH + (RSH+RLH).*(h*XL-XC./h);
E = VSH.*(K*RLH-XLH.*(h*XL-XC./h)) + ILH.*XTLH.*(h*XL-XC./h);
F = VSH.*(RLH.*(h*XL-XC./h)+K*XLH) - ILH.*RTLH.*(h*XL-XC./h);

ISH_comp = complex(A,B)./complex(C,D);
VLH_comp = complex(E,F)./complex(C,D);

ISH = abs(ISH_comp);
VLH = abs(VLH_comp);

% Find fundamental IS1 & VL1

ZLH1 = complex(RL1,XL1);
ZF1 = complex(0,XL-XC);
ZS1 = complex(RS1,XS1);

ZT1 = ZS1 + ZLH1*ZF1/(ZLH1+ZF1);

IS1_comp = VSH_f/ZT1;
IS1 = abs(IS1_comp); % Fundamental IS1
IS = sqrt(IS1^2+sum(ISH.^2));
VL1_comp = IS1_comp*ZLH1*ZF1/(ZLH1+ZF1);
VL1 = abs(VL1_comp); % Fundamental VL1
VL = sqrt(VL1^2+sum(VLH.^2));

GL1 = real(1/complex(RL1,XL1)); % Conductance for fundamental
GLH = real(1./complex(RLH,XLH)); % Conductance for each harmonic

PL1 = GL1*VL1^2;

DPF = PL1/(VL1*IS1)*100; % Displacement %PF

PL = PL1 + sum(GLH.*VLH.^2); % Load active power per phase in watt

PF = PL/(sqrt(IS1^2+sum(ISH.^2))*sqrt(VL1^2+sum(VLH.^2)))*100; % Compensated Load %PF

PLOSS = IS1^2*RS1 + sum((ISH.^2).*RSH); % Transmission loss
Efficiency = PL/(PL+PLOSS)*100;

VTHD = sqrt(sum(VLH.^2))/VL1; % VTHD p.u. value
ITHD = sqrt(sum(ISH.^2))/IS1; % ITHD p.u. value

VTHD_percent = VTHD*100;
ITHD_percent = ITHD*100;

HP = sqrt(VTHD^2 + ITHD^2)*100; % Harmonic pollution

if all(VLH./VL1*100<3) && (VTHD_percent<5) && (ITHD_percent<5) && (PF>=90) %&& (PF<=95)
fobj = -abs(VTHD_percent-5)-abs(ITHD_percent-5)+abs(PF-95);

else
fobj = 100;
end
end

% F2

function fobj = F2(filter_data)
format long g;
global VTHD_percent ITHD_percent HP Efficiency PF DPF IS VL PLOSS VLH ISH IS1 VL1;
VSH_f = 2400;
ILH_f = 987.83;
RL1 = 1.7421;
XL1 = 1.696;
RS1 = 0.02163;
XS1 = 0.2163;
h = [5 7 11 13];
VSH_percent = [0.02 0.015 0.01 0.005];
ILH_percent = [0.4 0.06 0.02 0.01];
XL = filter_data(1);
XC = filter_data(2);
K = filter_data(3);

VSH = VSH_f*VSH_percent;
RLH = RL1*ones(1,4);
ILH = ILH_f*ILH_percent;
XLH = XL1*h;
RSH = RS1*ones(1,4);
XSH = XS1*h;

RTLH = RLH.*RSH-XLH.*XSH;
XTLH = RLH.*XSH+RSH.*XLH;

A = VSH.*RLH - ILH.*XLH.*(h*XL-XC./h);
B = VSH.*(XLH+h*XL-XC./h) + ILH.*RLH.*(h*XL-XC./h);
C = RTLH + K*RLH - (XLH+XSH).*(h*XL-XC./h);
D = XTLH + K*XLH + (RSH+RLH).*(h*XL-XC./h);
E = VSH.*(K*RLH-XLH.*(h*XL-XC./h)) + ILH.*XTLH.*(h*XL-XC./h);
F = VSH.*(RLH.*(h*XL-XC./h)+K*XLH) - ILH.*RTLH.*(h*XL-XC./h);

ISH_comp = complex(A,B)./complex(C,D);
VLH_comp = complex(E,F)./complex(C,D);

ISH = abs(ISH_comp);
VLH = abs(VLH_comp);

% Find fundamental IS1 & VL1

ZLH1 = complex(RL1,XL1);
ZF1 = complex(0,XL-XC);
ZS1 = complex(RS1,XS1);

ZT1 = ZS1 + ZLH1*ZF1/(ZLH1+ZF1);

IS1_comp = VSH_f/ZT1;
IS1 = abs(IS1_comp); % Fundamental IS1
IS = sqrt(IS1^2+sum(ISH.^2));
VL1_comp = IS1_comp*ZLH1*ZF1/(ZLH1+ZF1);
VL1 = abs(VL1_comp); % Fundamental VL1
VL = sqrt(VL1^2+sum(VLH.^2));

GL1 = real(1/complex(RL1,XL1)); % Conductance for fundamental
GLH = real(1./complex(RLH,XLH)); % Conductance for each harmonic

PL1 = GL1*VL1^2;

DPF = PL1/(VL1*IS1)*100; % Displacement %PF

PL = PL1 + sum(GLH.*VLH.^2); % Load active power per phase in watt

PF = PL/(sqrt(IS1^2+sum(ISH.^2))*sqrt(VL1^2+sum(VLH.^2)))*100; % Compensated Load %PF

PLOSS = IS1^2*RS1 + sum((ISH.^2).*RSH); % Transmission loss
Efficiency = PL/(PL+PLOSS)*100;

VTHD = sqrt(sum(VLH.^2))/VL1; % VTHD p.u. value
ITHD = sqrt(sum(ISH.^2))/IS1; % ITHD p.u. value

VTHD_percent = VTHD*100;
ITHD_percent = ITHD*100;

HP = sqrt(VTHD^2 + ITHD^2)*100; % Harmonic pollution

if all(VLH./VL1*100<3) && (VTHD_percent<5) && (ITHD_percent<5) && (PF>=90) %&& (PF<=95)
fobj = -abs(VTHD_percent-5)-abs(ITHD_percent-5)+abs(PF-95);

else
fobj = 100;
end
end

% F3

function fobj = F3(filter_data)
format long g;
global VTHD_percent ITHD_percent HP Efficiency PF DPF IS VL PLOSS VLH ISH IS1 VL1;
VSH_f = 2400;
ILH_f = 987.83;
RL1 = 1.7421;
XL1 = 1.696;
RS1 = 0.02163;
XS1 = 0.2163;
h = [5 7 11 13];
VSH_percent = [0.04 0.03 0.02 0.01];
ILH_percent = [0.4 0.06 0.02 0.01];

%filter_data = [0.99*10^-9 2.6171 8.07];
XL = filter_data(1);
XC = filter_data(2);
K = filter_data(3);

VSH = VSH_f*VSH_percent;
RLH = RL1*ones(1,4);
ILH = ILH_f*ILH_percent;
XLH = XL1*h;
RSH = RS1*ones(1,4);
XSH = XS1*h;

RTLH = RLH.*RSH-XLH.*XSH;
XTLH = RLH.*XSH+RSH.*XLH;

A = VSH.*RLH - ILH.*XLH.*(h*XL-XC./h);
B = VSH.*(XLH+h*XL-XC./h) + ILH.*RLH.*(h*XL-XC./h);
C = RTLH + K*RLH - (XLH+XSH).*(h*XL-XC./h);
D = XTLH + K*XLH + (RSH+RLH).*(h*XL-XC./h);
E = VSH.*(K*RLH-XLH.*(h*XL-XC./h)) + ILH.*XTLH.*(h*XL-XC./h);
F = VSH.*(RLH.*(h*XL-XC./h)+K*XLH) - ILH.*RTLH.*(h*XL-XC./h);

ISH_comp = complex(A,B)./complex(C,D);
VLH_comp = complex(E,F)./complex(C,D);

ISH = abs(ISH_comp);
VLH = abs(VLH_comp);

% Find fundamental IS1 & VL1

ZLH1 = complex(RL1,XL1);
ZF1 = complex(0,XL-XC);
ZS1 = complex(RS1,XS1);

ZT1 = ZS1 + ZLH1*ZF1/(ZLH1+ZF1);

IS1_comp = VSH_f/ZT1;
IS1 = abs(IS1_comp); % Fundamental IS1
IS = sqrt(IS1^2+sum(ISH.^2));
VL1_comp = IS1_comp*ZLH1*ZF1/(ZLH1+ZF1);
VL1 = abs(VL1_comp); % Fundamental VL1
VL = sqrt(VL1^2+sum(VLH.^2));

GL1 = real(1/complex(RL1,XL1)); % Conductance for fundamental
GLH = real(1./complex(RLH,XLH)); % Conductance for each harmonic

PL1 = GL1*VL1^2;

DPF = PL1/(VL1*IS1)*100; % Displacement %PF

PL = PL1 + sum(GLH.*VLH.^2); % Load active power per phase in watt

PF = PL/(sqrt(IS1^2+sum(ISH.^2))*sqrt(VL1^2+sum(VLH.^2)))*100; % Compensated Load %PF

PLOSS = IS1^2*RS1 + sum((ISH.^2).*RSH); % Transmission loss
Efficiency = PL/(PL+PLOSS)*100;

VTHD = sqrt(sum(VLH.^2))/VL1; % VTHD p.u. value
ITHD = sqrt(sum(ISH.^2))/IS1; % ITHD p.u. value

VTHD_percent = VTHD*100;
ITHD_percent = ITHD*100;

HP = sqrt(VTHD^2 + ITHD^2)*100; % Harmonic pollution

if all(VLH./VL1*100<3) && (VTHD_percent<5) && (ITHD_percent<5) && (PF>=90) %&& (PF<=95)
fobj = -abs(VTHD_percent-5)-abs(ITHD_percent-5)+abs(PF-95);

else
fobj = 100;
end
end

% F4

function fobj = F4(filter_data)
format long g;
global VTHD_percent ITHD_percent HP Efficiency PF DPF IS VL PLOSS VLH ISH IS1 VL1;
VSH_f = 2400;
ILH_f = 987.83;
RL1 = 1.7421;
XL1 = 1.696;
RS1 = 0.02163;
XS1 = 0.2163;
h = [5 7 11 13];
% VSH_percent = [0.038 0.03 0.015 0.01]; %1st assumption
% VSH_percent = [0.04 0.03 0.03 0.012]; %2nd assumption
VSH_percent = [0.04 0.03 0.03 0.012]; %3rd assumption
ILH_percent = [0.4 0.06 0.02 0.01];
%ILH_percent = [0.3 0.05 0.01 0.005]; %3rd assumption

% filter_data = [1.40636*10^-14 5.22399 17.1476];
XL = filter_data(1);
XC = filter_data(2);
K = filter_data(3);

VSH = VSH_f*VSH_percent;
RLH = RL1*ones(1,4);
ILH = ILH_f*ILH_percent;
XLH = XL1*h;
RSH = RS1*ones(1,4);
XSH = XS1*h;

RTLH = RLH.*RSH-XLH.*XSH;
XTLH = RLH.*XSH+RSH.*XLH;

A = VSH.*RLH - ILH.*XLH.*(h*XL-XC./h);
B = VSH.*(XLH+h*XL-XC./h) + ILH.*RLH.*(h*XL-XC./h);
C = RTLH + K*RLH - (XLH+XSH).*(h*XL-XC./h);
D = XTLH + K*XLH + (RSH+RLH).*(h*XL-XC./h);
E = VSH.*(K*RLH-XLH.*(h*XL-XC./h)) + ILH.*XTLH.*(h*XL-XC./h);
F = VSH.*(RLH.*(h*XL-XC./h)+K*XLH) - ILH.*RTLH.*(h*XL-XC./h);

ISH_comp = complex(A,B)./complex(C,D);
VLH_comp = complex(E,F)./complex(C,D);

ISH = abs(ISH_comp);
VLH = abs(VLH_comp);

% Find fundamental IS1 & VL1

ZLH1 = complex(RL1,XL1);
ZF1 = complex(0,XL-XC);
ZS1 = complex(RS1,XS1);

ZT1 = ZS1 + ZLH1*ZF1/(ZLH1+ZF1);

IS1_comp = VSH_f/ZT1;
IS1 = abs(IS1_comp); % Fundamental IS1
IS = sqrt(IS1^2+sum(ISH.^2));
VL1_comp = IS1_comp*ZLH1*ZF1/(ZLH1+ZF1);
VL1 = abs(VL1_comp); % Fundamental VL1
VL = sqrt(VL1^2+sum(VLH.^2));

GL1 = real(1/complex(RL1,XL1)); % Conductance for fundamental
GLH = real(1./complex(RLH,XLH)); % Conductance for each harmonic

PL1 = GL1*VL1^2;

DPF = PL1/(VL1*IS1)*100; % Displacement %PF

PL = PL1 + sum(GLH.*VLH.^2); % Load active power per phase in watt

PF = PL/(sqrt(IS1^2+sum(ISH.^2))*sqrt(VL1^2+sum(VLH.^2)))*100; % Compensated Load %PF

PLOSS = IS1^2*RS1 + sum((ISH.^2).*RSH); % Transmission loss
Efficiency = PL/(PL+PLOSS)*100;

VTHD = sqrt(sum(VLH.^2))/VL1; % VTHD p.u. value
ITHD = sqrt(sum(ISH.^2))/IS1; % ITHD p.u. value

VTHD_percent = VTHD*100;
ITHD_percent = ITHD*100;

HP = sqrt(VTHD^2 + ITHD^2)*100; % Harmonic pollution

if all(VLH./VL1*100<3) && (VTHD_percent<5) && (ITHD_percent<5) && (PF>=90) %&& (PF<=95)
fobj = -abs(VTHD_percent-5)-abs(ITHD_percent-5)+abs(PF-95);

else
fobj = 1;
end
end

% F5

function fobj = F5(filter_data)
format long g;
global VTHD_percent ITHD_percent HP Efficiency PF DPF IS VL PLOSS VLH ISH IS1 VL1;
VSH_f = 2400;
ILH_f = 987.83;
RL1 = 1.7421;
XL1 = 1.696;
RS1 = 0.02163;
XS1 = 0.2163;
h = [5 7 11 13];
VSH_percent = [0 0 0 0];
ILH_percent = [0.4 0.06 0.02 0.01];

%filter_data = [97.93/1000 2.7061 12];
XL = filter_data(1);
XC = filter_data(2);
K = filter_data(3);

VSH = VSH_f*VSH_percent;
RLH = RL1*ones(1,4);
ILH = ILH_f*ILH_percent;
XLH = XL1*h;
RSH = RS1*ones(1,4);
XSH = XS1*h;

RTLH = RLH.*RSH-XLH.*XSH;
XTLH = RLH.*XSH+RSH.*XLH;

A = VSH.*RLH - ILH.*XLH.*(h*XL-XC./h);
B = VSH.*(XLH+h*XL-XC./h) + ILH.*RLH.*(h*XL-XC./h);
C = RTLH + K*RLH - (XLH+XSH).*(h*XL-XC./h);
D = XTLH + K*(XLH+h*XL-XC./h) + (RSH+RLH).*(h*XL-XC./h);
%E = -VSH.*XLH.*(h*XL-XC./h) + ILH.*(XTLH+K*XLH).*(h*XL-XC./h);%wrong VLH'
E = VSH.*(K*RLH-XLH.*(h*XL-XC./h)) + ILH.*XTLH.*(h*XL-XC./h);
%F = VSH.*RLH.*(h*XL-XC./h) - ILH.*(RTLH+K*RLH).*(h*XL-XC./h);%wrong VLH'
F = VSH.*((K+RLH).*(h*XL-XC./h)+K*XLH) - ILH.*RTLH.*(h*XL-XC./h);

ISH_comp = complex(A,B)./complex(C,D);
VLH_comp = complex(E,F)./complex(C,D);

ISH = abs(ISH_comp);
VLH = abs(VLH_comp);

% Find fundamental IS1 & VL1

ZLH1 = complex(RL1,XL1);
ZF1 = complex(0,XL-XC);
ZS1 = complex(RS1,XS1);

ZT1 = ZS1 + ZLH1*ZF1/(ZLH1+ZF1);

IS1_comp = VSH_f/ZT1;
IS1 = abs(IS1_comp); % Fundamental IS1
IS = sqrt(IS1^2+sum(ISH.^2));
VL1_comp = IS1_comp*ZLH1*ZF1/(ZLH1+ZF1);
VL1 = abs(VL1_comp); % Fundamental VL1
VL = sqrt(VL1^2+sum(VLH.^2));

GL1 = real(1/complex(RL1,XL1)); % Conductance for fundamental
GLH = real(1./complex(RLH,XLH)); % Conductance for each harmonic

PL1 = GL1*VL1^2;

DPF = PL1/(VL1*IS1)*100; % Displacement %PF

PL = PL1 + sum(GLH.*VLH.^2); % Load active power per phase in watt

PF = PL/(sqrt(IS1^2+sum(ISH.^2))*sqrt(VL1^2+sum(VLH.^2)))*100; % Compensated Load %PF

PLOSS = IS1^2*RS1 + sum((ISH.^2).*RSH); % Transmission loss
Efficiency = PL/(PL+PLOSS)*100;

VTHD = sqrt(sum(VLH.^2))/VL1; % VTHD p.u. value
ITHD = sqrt(sum(ISH.^2))/IS1; % ITHD p.u. value

VTHD_percent = VTHD*100;
ITHD_percent = ITHD*100;

HP = sqrt(VTHD^2 + ITHD^2)*100; % Harmonic pollution

if all(VLH./VL1*100<3) && (VTHD_percent<5) && (ITHD_percent<5) && (PF>=90) %&& (PF<=95)
fobj = -abs(VTHD_percent-5)-abs(ITHD_percent-5)+abs(PF-95);

else
fobj = 100;
end
end

% F6

function fobj = F6(filter_data)
format long g;
global VTHD_percent ITHD_percent HP Efficiency PF DPF IS VL PLOSS VLH ISH IS1 VL1;
VSH_f = 2400;
ILH_f = 987.83;
RL1 = 1.7421;
XL1 = 1.696;
RS1 = 0.02163;
XS1 = 0.2163;
h = [5 7 11 13];
VSH_percent = [0.02 0.015 0.01 0.005];
ILH_percent = [0.4 0.06 0.02 0.01];

% filter_data = [97.93/1000 2.7061 12];
XL = filter_data(1);
XC = filter_data(2);
K = filter_data(3);

VSH = VSH_f*VSH_percent;
RLH = RL1*ones(1,4);
ILH = ILH_f*ILH_percent;
XLH = XL1*h;
RSH = RS1*ones(1,4);
XSH = XS1*h;

RTLH = RLH.*RSH-XLH.*XSH;
XTLH = RLH.*XSH+RSH.*XLH;

A = VSH.*RLH - ILH.*XLH.*(h*XL-XC./h);
B = VSH.*(XLH+h*XL-XC./h) + ILH.*RLH.*(h*XL-XC./h);
C = RTLH + K*RLH - (XLH+XSH).*(h*XL-XC./h);
D = XTLH + K*(XLH+h*XL-XC./h) + (RSH+RLH).*(h*XL-XC./h);
%E = -VSH.*XLH.*(h*XL-XC./h) + ILH.*(XTLH+K*XLH).*(h*XL-XC./h);%wrong VLH'
E = VSH.*(K*RLH-XLH.*(h*XL-XC./h)) + ILH.*XTLH.*(h*XL-XC./h);
%F = VSH.*RLH.*(h*XL-XC./h) - ILH.*(RTLH+K*RLH).*(h*XL-XC./h);%wrong VLH'
F = VSH.*((K+RLH).*(h*XL-XC./h)+K*XLH) - ILH.*RTLH.*(h*XL-XC./h);

ISH_comp = complex(A,B)./complex(C,D);
VLH_comp = complex(E,F)./complex(C,D);

ISH = abs(ISH_comp);
VLH = abs(VLH_comp);

% Find fundamental IS1 & VL1

ZLH1 = complex(RL1,XL1);
ZF1 = complex(0,XL-XC);
ZS1 = complex(RS1,XS1);

ZT1 = ZS1 + ZLH1*ZF1/(ZLH1+ZF1);

IS1_comp = VSH_f/ZT1;
IS1 = abs(IS1_comp); % Fundamental IS1
IS = sqrt(IS1^2+sum(ISH.^2));
VL1_comp = IS1_comp*ZLH1*ZF1/(ZLH1+ZF1);
VL1 = abs(VL1_comp); % Fundamental VL1
VL = sqrt(VL1^2+sum(VLH.^2));

GL1 = real(1/complex(RL1,XL1)); % Conductance for fundamental
GLH = real(1./complex(RLH,XLH)); % Conductance for each harmonic

PL1 = GL1*VL1^2;

DPF = PL1/(VL1*IS1)*100; % Displacement %PF

PL = PL1 + sum(GLH.*VLH.^2); % Load active power per phase in watt

PF = PL/(sqrt(IS1^2+sum(ISH.^2))*sqrt(VL1^2+sum(VLH.^2)))*100; % Compensated Load %PF

PLOSS = IS1^2*RS1 + sum((ISH.^2).*RSH); % Transmission loss
Efficiency = PL/(PL+PLOSS)*100;

VTHD = sqrt(sum(VLH.^2))/VL1; % VTHD p.u. value
ITHD = sqrt(sum(ISH.^2))/IS1; % ITHD p.u. value

VTHD_percent = VTHD*100;
ITHD_percent = ITHD*100;

HP = sqrt(VTHD^2 + ITHD^2)*100; % Harmonic pollution

if all(VLH./VL1*100<3) && (VTHD_percent<5) && (ITHD_percent<5) && (PF>=90) %&& (PF<=95)
fobj = -abs(VTHD_percent-5)-abs(ITHD_percent-5)+abs(PF-95);

else
fobj = 1;
end
end

% F7

function fobj = F7(filter_data)
format long g;
global VTHD_percent ITHD_percent HP Efficiency PF DPF IS VL PLOSS VLH ISH IS1 VL1;
VSH_f = 2400;
ILH_f = 987.83;
RL1 = 1.7421;
XL1 = 1.696;
RS1 = 0.02163;
XS1 = 0.2163;
h = [5 7 11 13];
VSH_percent = [0.04 0.03 0.02 0.01];
ILH_percent = [0.4 0.06 0.02 0.01];

% filter_data = [9.9*10^-10 2.6171 8.07];
XL = filter_data(1);
XC = filter_data(2);
K = filter_data(3);

VSH = VSH_f*VSH_percent;
RLH = RL1*ones(1,4);
ILH = ILH_f*ILH_percent;
XLH = XL1*h;
RSH = RS1*ones(1,4);
XSH = XS1*h;

RTLH = RLH.*RSH-XLH.*XSH;
XTLH = RLH.*XSH+RSH.*XLH;

A = VSH.*RLH - ILH.*XLH.*(h*XL-XC./h);
B = VSH.*(XLH+h*XL-XC./h) + ILH.*RLH.*(h*XL-XC./h);
C = RTLH + K*RLH - (XLH+XSH).*(h*XL-XC./h);
D = XTLH + K*(XLH+h*XL-XC./h) + (RSH+RLH).*(h*XL-XC./h);
%E = -VSH.*XLH.*(h*XL-XC./h) + ILH.*(XTLH+K*XLH).*(h*XL-XC./h);%wrong VLH'
E = VSH.*(K*RLH-XLH.*(h*XL-XC./h)) + ILH.*XTLH.*(h*XL-XC./h);
%F = VSH.*RLH.*(h*XL-XC./h) - ILH.*(RTLH+K*RLH).*(h*XL-XC./h);%wrong VLH'
F = VSH.*((K+RLH).*(h*XL-XC./h)+K*XLH) - ILH.*RTLH.*(h*XL-XC./h);

ISH_comp = complex(A,B)./complex(C,D);
VLH_comp = complex(E,F)./complex(C,D);

ISH = abs(ISH_comp);
VLH = abs(VLH_comp);

% Find fundamental IS1 & VL1

ZLH1 = complex(RL1,XL1);
ZF1 = complex(0,XL-XC);
ZS1 = complex(RS1,XS1);

ZT1 = ZS1 + ZLH1*ZF1/(ZLH1+ZF1);

IS1_comp = VSH_f/ZT1;
IS1 = abs(IS1_comp); % Fundamental IS1
IS = sqrt(IS1^2+sum(ISH.^2));
VL1_comp = IS1_comp*ZLH1*ZF1/(ZLH1+ZF1);
VL1 = abs(VL1_comp); % Fundamental VL1
VL = sqrt(VL1^2+sum(VLH.^2));

GL1 = real(1/complex(RL1,XL1)); % Conductance for fundamental
GLH = real(1./complex(RLH,XLH)); % Conductance for each harmonic

PL1 = GL1*VL1^2;

DPF = PL1/(VL1*IS1)*100; % Displacement %PF

PL = PL1 + sum(GLH.*VLH.^2); % Load active power per phase in watt

PF = PL/(sqrt(IS1^2+sum(ISH.^2))*sqrt(VL1^2+sum(VLH.^2)))*100; % Compensated Load %PF

PLOSS = IS1^2*RS1 + sum((ISH.^2).*RSH); % Transmission loss
Efficiency = PL/(PL+PLOSS)*100;

VTHD = sqrt(sum(VLH.^2))/VL1; % VTHD p.u. value
ITHD = sqrt(sum(ISH.^2))/IS1; % ITHD p.u. value

VTHD_percent = VTHD*100;
ITHD_percent = ITHD*100;

HP = sqrt(VTHD^2 + ITHD^2)*100; % Harmonic pollution

if all(VLH./VL1*100<3) && (VTHD_percent<5) && (ITHD_percent<5) && (PF>=90) %&& (PF<=95)
fobj = -abs(VTHD_percent-5)-abs(ITHD_percent-5)+abs(PF-95);

else
fobj = 100;
end
end

% F8

function fobj = F8(filter_data)
format long g;
global VTHD_percent ITHD_percent HP Efficiency PF DPF IS VL PLOSS VLH ISH IS1 VL1;
VSH_f = 2400;
ILH_f = 987.83;
RL1 = 1.7421;
% RL1 = 0.5421;
XL1 = 1.696;
%XL1 = 5.696;
RS1 = 0.02163;
XS1 = 0.2163;
h = [5 7 11 13];
% VSH_percent = [0.038 0.03 0.015 0.01]; %1st assumption
% VSH_percent = [0.04 0.03 0.03 0.012]; %2nd assumption
VSH_percent = [0.04 0.03 0.03 0.012]; %3rd assumption
ILH_percent = [0.4 0.06 0.02 0.01];
% ILH_percent = [0.3 0.05 0.01 0.005]; %3rd assumption

%filter_data = [97.93/1000 2.7061 12];
XL = filter_data(1);
XC = filter_data(2);
K = filter_data(3);

VSH = VSH_f*VSH_percent;
RLH = RL1*ones(1,4);
ILH = ILH_f*ILH_percent;
XLH = XL1*h;
RSH = RS1*ones(1,4);
XSH = XS1*h;

RTLH = RLH.*RSH-XLH.*XSH;
XTLH = RLH.*XSH+RSH.*XLH;

A = VSH.*RLH - ILH.*XLH.*(h*XL-XC./h);
B = VSH.*(XLH+h*XL-XC./h) + ILH.*RLH.*(h*XL-XC./h);
C = RTLH + K*RLH - (XLH+XSH).*(h*XL-XC./h);
D = XTLH + K*(XLH+h*XL-XC./h) + (RSH+RLH).*(h*XL-XC./h);
%E = -VSH.*XLH.*(h*XL-XC./h) + ILH.*(XTLH+K*XLH).*(h*XL-XC./h);%wrong VLH'
E = VSH.*(K*RLH-XLH.*(h*XL-XC./h)) + ILH.*XTLH.*(h*XL-XC./h);
%F = VSH.*RLH.*(h*XL-XC./h) - ILH.*(RTLH+K*RLH).*(h*XL-XC./h);%wrong VLH'
F = VSH.*((K+RLH).*(h*XL-XC./h)+K*XLH) - ILH.*RTLH.*(h*XL-XC./h);

ISH_comp = complex(A,B)./complex(C,D);
VLH_comp = complex(E,F)./complex(C,D);

ISH = abs(ISH_comp);
VLH = abs(VLH_comp);

% Find fundamental IS1 & VL1

ZLH1 = complex(RL1,XL1);
ZF1 = complex(0,XL-XC);
ZS1 = complex(RS1,XS1);

ZT1 = ZS1 + ZLH1*ZF1/(ZLH1+ZF1);

IS1_comp = VSH_f/ZT1;
IS1 = abs(IS1_comp); % Fundamental IS1
IS = sqrt(IS1^2+sum(ISH.^2));
VL1_comp = IS1_comp*ZLH1*ZF1/(ZLH1+ZF1);
VL1 = abs(VL1_comp); % Fundamental VL1
VL = sqrt(VL1^2+sum(VLH.^2));

GL1 = real(1/complex(RL1,XL1)); % Conductance for fundamental
GLH = real(1./complex(RLH,XLH)); % Conductance for each harmonic

PL1 = GL1*VL1^2;

DPF = PL1/(VL1*IS1)*100; % Displacement %PF

PL = PL1 + sum(GLH.*VLH.^2); % Load active power per phase in watt

PF = PL/(sqrt(IS1^2+sum(ISH.^2))*sqrt(VL1^2+sum(VLH.^2)))*100; % Compensated Load %PF

PLOSS = IS1^2*RS1 + sum((ISH.^2).*RSH); % Transmission loss
Efficiency = PL/(PL+PLOSS)*100;

VTHD = sqrt(sum(VLH.^2))/VL1; % VTHD p.u. value
ITHD = sqrt(sum(ISH.^2))/IS1; % ITHD p.u. value

VTHD_percent = VTHD*100;
ITHD_percent = ITHD*100;

HP = sqrt(VTHD^2 + ITHD^2)*100; % Harmonic pollution

if all(VLH./VL1*100<3) && (VTHD_percent<5) && (ITHD_percent<5) && (PF>=90) %&& (PF<=95)
fobj = -abs(VTHD_percent-5)-abs(ITHD_percent-5)+abs(PF-95);

else
fobj = 1;
end

end

