%Expanded Comprehensive Disinfection and Water Quality Model - Batch Version (CDWQ-E2 BATCH)

%Originally written by: John Woolschlager
%Northwestern University, Evanston, USA
%Expanded and modified by: Precious Biyela
%Arizona State University, Tempe, Arizona, USA
%Expanded and modified by: Liam Culligan
%University of the Witwatersrand, Johannesburg, South Africa
%August 2014

clear all
clc
format short e

%Set run time and time step

days = 30;
timeInc = 0.05;
timeSteps = (1/timeInc);

%Set global conditions for distribution system
%Maximum of 1 vector value

logHArr = [-7.5]; 
CtCO3Arr = [1.0E-3, 3.0E-3, 5.0E-3];
TempArr = [20];

numLogH = numel(logHArr);
numCtCO3 = numel(CtCO3Arr);
numTemp = numel(TempArr);

%Display prompt for graphs plots

promptMessage = sprintf('Do you want to plot the graphs?');
promptButton = questdlg(promptMessage, 'Plot Graphs?', 'Yes', 'No', 'Yes');
drawnow;

%Determine number of iterations based on number of varying global
%conditions

numIterations = max([numLogH,numCtCO3,numTemp]);

%Determine which global condition is varying

if numIterations > 1
    if (numLogH == numIterations && numCtCO3 == 1 && numTemp == 1)
        varying = 'logH';
    elseif (numCtCO3 == numIterations && numLogH == 1 && numTemp == 1)
        varying = 'CtCO3';
    elseif (numTemp == numIterations && numLogH == 1 && numCtCO3 == 1)
        varying = 'Temp';
    else
        error('Only one varying parameter is allowed.');
    end
else
    varying = '';
end

iterationString = '';

tV = 1:(days*24*timeSteps);

%Loop through each iteration

for iteration = 1:numIterations
    
    %Set conversion factors
    
    moleN = 14E6;
    MC = 3.1211E-8;
    MN = 6.2422E-9;
    
    %Set biological parameters (per hour)
    
    bh = 0.1/24.0;
    bn1 = 0.05/24.0;
    bn2 = 0.05/24.0;
    fd = 0.8;
    kbaph = 0.1/24.0;
    kbapn1 = 0.1/24.0;
    kbapn2 = 0.1/24.0;
    kuaph = 0.2;
    kuapn1 = (1.54E6)*MN;
    kuapn2 = (4.2E5)*MN;
    kepsh = 7.5E-03;
    kepsn1 = 6.6E-03;
    kepsn2 = 1.8E-03;
    Kbom1 = 15000.0;
    Kbom2 = 120000.0;
    khydEPS = 0.17/24.0;
    KNH3 = 2.14E-6; 
    KNCl = 1000/moleN;
    KNO2 = 5.36E-5;
    KNO3 = 1000/moleN;
    KNH2Cl = 2000.0/moleN;
    KNHOHCl = 1000.0/moleN;
    KIDO = 0.05*1000;
    Kbap = 30000.0;
    Kuap = 20000.0;
    KO2 = 1000*0.20;
    qh = 10.0/24;
    qbom1 = 10.0/24;
    qbom2 = 10.0/24;
    qn1 = 1.6/(24*moleN);
    qn2 = 7.0/(24*moleN);
    qbap = 2.0/24;
    quap = 13.0/24;
    qNH2Cl = 1.6/(24.0*moleN);
    qNHOHCl = 1.6/(24.0*moleN);
    qNO2 = 7.0/(24.0*moleN);
    and1 = 0.6; 
    and2 = 0.6; 
    Yh = 0.6;
    Yn1 = 6.16E6;
    Yn2 = 1.68E6;
    Yp = 0.6;
    
    %Set chemical conditions
    
    kH = 2.5E7;
    kH2CO3 = 4.0E4;
    kHCO3 = 14.4E2;
    kox4a = 2.3E-2; 
    kox4b = 5.0; 
 
    kd1x1 = 0.0; 
    kd1x2 = 0.0;
    kd1x3 = 0.0;
    kd1x4 = 0.0;
    kd1x5 = 0.0;
    kd1x6 = 0.0;
    kd2x1 = 0.0;
    kd2x2 = 0.0;
    kd2x3 = 0.0;
    kd2x4 = 0.0;
    kd2x5 = 0.0;
    kd2x6 = 0.0;
    kd3x1 = 600.0; 
    kd3x2 = 300.0;
    kd3x3 = 300.0; 
    kd3x4 = 150.0;
    kd3x5 = 300.0; 
    kd3x6 = 150.0;
    
    kox1x1 = 0.0; 
    kox1x2 = 0.0; 
    kox2x1 = 0.0; 
    kox2x2 = 0.0; 
    kox3x1 = 30.0;
    kox3x2 = 3.0; 
     
    kox5a = 1.36E7; 
    kox5b = 2.17E2;
    kox5c = 5.5E5;
    
    kAD3 = 1.0E6; 
    kAD2 = 2.16E8; 
    kAD4 = 6.0E5; 
    kAD5 = 5.0E-3; 
    
    %Set physical and temperature parameters
    
    Tbio = 1.05; 
    Tbio2 = 1.00; 
    Tchem = 1.05; 
    
    %Set legend titles for graphs
    
    if strcmp(varying,'logH')
        logH = logHArr(iteration);
        iterationString{iteration} = strcat('pH=', num2str(abs(logH)));
    else
        logH = logHArr(1);
    end
    if strcmp(varying,'CtCO3')
        CtCO3 = CtCO3Arr(iteration);
        iterationString{iteration} = strcat('CtCO_3=', num2str(CtCO3*1000),'mM');
    else
        CtCO3 = CtCO3Arr(1);
    end
    if strcmp(varying,'Temp')
        Temp = TempArr(iteration);
        iterationString{iteration} = strcat('Temp=', num2str(Temp),'^oC');
    else
        Temp = TempArr(1);
    end
    
    %Determine oxygen saturation value
    
    kla = 1/24; 
    O2s = (14.59 - 0.3955*Temp + 0.0072*Temp^2 - 0.0000619*Temp^3)*1000; %Oxygen saturation - ugO2/l
    
    %Set initial conditions
    
    BOM1 = 900.0; %ugCOD/l --- Divide by 2670 to convert to mgC/l
    BOM2 = 4100.0; %ugCOD/l --- Divide by 2670 to convert to mgC/l
    CtNH3 = 2.28E-5;%mol N/l --- Multiply by 17543.9 to convert to mgN/l
    NO2 = 2.28E-6; %mol N/l --- Multiply by 17543.9 to convert to mgN/l
    NO3 = 3.80E-5; %mol N/l --- Multiply by 17543.9 to convert to mgN/l
    Xhs = 0.4; %ugCOD/l --- Multiply by 2404 to convert to cells/l
    Xn1s = 0.04; %ugCOD/l --- Multiply by 2404 to convert to cells/l
    Xn2s = 0.04;%ugCOD/l --- Multiply by 2404 to convert to cells/l
    UAP = 0.0; %ugCOD/l --- Divide by 2670 to convert to mgC/l
    BAP = 0.0; %ugCOD/l --- Divide by 2670 to convert to mgC/l
    EPSs = 0.0; %ugCOD/l --- Divide by 2670 to convert to mgC/l
    Xis = 0.0; %ugCOD/l --- Multiply by 2404 to convert to cells/l
    CO2 = 10.0; %ugCOD/l --- Divide by 2670 to convert to mgC/l
    CtOCl = 0.0; %mol Cl/l
    NH2Cl = 4.4E-5; %mol N/l --- Multiply by 70909.1 to convert to mgCl2/l
    NHCl2 = 0; %mol N/l --- Multiply by 70909.1 to convert to mgCl2/l
    NHOHCl = 0.0; %mol N/l --- Multiply by 70909.1 to convert to mgCl2/l
    Cl = 0.0; %mol Cl/l
    N2 = 0.0; %mol N/l
    O2 = O2s; %ug/l
    
    %Set initial acid-base equilibrium and temperature adjustment factors
    
    H = 10^logH;
    R = 8.314;
    KeNH3 = 7.029E-1*(exp(-52210/(R*(Temp + 273.15))));
    KeHCO3 = 1.119E-5*(exp(-7700/(R*(Temp + 273.15))));
    KeCO3 = 2.044E-8*(exp(-14900/(R*(Temp + 273.15))));
    KeOCl = (8.273E-6)*exp(-13800/(R*(Temp+273.15)));
    kCle = 0.249*exp(-7290/(Temp+273.15));
    
    H2CO3 = CtCO3*((H*H)/(H*H + H*KeHCO3 + KeHCO3*KeCO3));
    HCO3 = CtCO3*((H*KeHCO3)/(H*H + H*KeHCO3 + KeHCO3*KeCO3));
    CO3 = CtCO3*((KeHCO3*KeCO3)/(H*H + H*KeHCO3 + KeHCO3*KeCO3));
    OH = (1.0D-14)/H;
    
    TFb = Tbio^(Temp-20.0); 
    TFb2 = Tbio2^(Temp-20.0);
    TFc = Tchem^(Temp-20.0);
    
    %Initialise time-specific individual loss and production term variables
    
    currT = 0;
    
    ad1(iteration,1) = 0;
    ad2(iteration,1) = 0;
    ad3(iteration,1) = 0;
    ad5(iteration,1) = 0;
    ox5(iteration,1) = 0;
    comet(iteration,1) = 0;
    organicox(iteration,1) = 0;
    
    bom1a(iteration,1) = 0;
    bom1b(iteration,1) = 0;
    bom1c(iteration,1) = 0;
    bom1d(iteration,1) = 0;
    
    bom2a(iteration,1) = 0;
    bom2b(iteration,1) = 0;
    bom2c(iteration,1) = 0;
    bom2d(iteration,1) = 0;
    bom2e(iteration,1) = 0;
    
    %Loop through day, hour & time steps
    
    for day = 1:days
        for hour = 1:24
            for t = 1:timeSteps
                
                currT = currT + 1;
                
                %Set acid-base equilibrium
                
                NH4 = (CtNH3*H)/(H+KeNH3);
                if NH4  <= 0.0
                    NH4 = 1E-18;
                end
                
                NH3 = (CtNH3*KeNH3)/(H+KeNH3);
                if NH3  <= 0.0
                    NH3 = 1E-18;
                end
                
                HOCl = (CtOCl*H)/(H+KeOCl);
                if HOCl < 0
                    HOCl = 0.0;
                end
                
                OCl = (CtOCl*KeOCl)/(H+KeOCl);
                if OCl < 0
                    OCl = 0.0;
                end
                
                %Adjust flexible kinetic constants
                
                kAD1 = kH*H+kH2CO3*H2CO3; 
                kox5 = (kox5a*H*(1+kox5b*NO2))/(kox5c*NH3 ...
                    + (1+kox5b*NO2));
                
                kox4 = (kox4a + kox4b*NO2)/OH;
                
                %Substrate utilisation rates
                
                UTLbom1 = qh*(BOM1/(Kbom1+BOM1))*(O2/(KO2+O2));
                UTLbom2 = qh*(BOM2/(Kbom2+BOM2))*(O2/(KO2+O2));
                UTLuap = quap*(UAP/(Kuap+UAP))*(O2/(KO2+O2));
                UTLbap = qbap*(BAP/(Kbap+BAP))*(O2/(KO2+O2));
                UTLNH3 = qn1*(NH3/(KNH3*(1+(NH2Cl/KNH2Cl)) + NH3))*(O2/(KO2+O2));
                if CO2 <= 0
                    UTLNH3 = 0;
                end
                UTLNH2Cl = qNH2Cl*(NH2Cl/(KNH2Cl*(1+(NH3/KNH3)) + NH2Cl))*(O2/(KO2+O2)); 
                UTLNHOHCl = qNHOHCl*(NHOHCl/(KNHOHCl + NHOHCl))*(O2/(KO2+O2));
                UTLNO2 = qn2*(NO2/(KNO2+NO2))*(O2/(KO2+O2));
                if CO2 <= 0
                    UTLNO2 = 0;
                end
                UTLbom = UTLbom1 + UTLbom2;
                UTLsmp = UTLuap + UTLbap;
                and1bom1 = qbom1*and1*(BOM1/(Kbom1 + BOM1))*(KIDO/(KIDO + O2))*(NO2/(KNO2 + NO2)); 
                and1bom2 = qbom2*and1*(BOM2/(Kbom2 + BOM2))*(KIDO/(KIDO + O2))*(NO2/(KNO2 + NO2)); 
                and1bom = and1bom1 + and1bom2;
                and2bom1 = qbom1*and2*(BOM1/(Kbom1 + BOM1))*(KIDO/(KIDO + O2))*(NO3/(KNO3 + NO3)); 
                and2bom2 = qbom1*and2*(BOM2/(Kbom2 + BOM2))*(KIDO/(KIDO + O2))*(NO3/(KNO3 + NO3)); 
                and2bom = and2bom1 + and2bom2;
                andbom = and1bom +and2bom;
                and1uap = quap*and1*(UAP/(Kuap + UAP))*(KIDO/(KIDO + O2))*(NO2/(KNO2 + NO2)); 
                and1bap = qbap*and1*(BAP/(Kbap + BAP))*(KIDO/(KIDO + O2))*(NO2/(KNO2 + NO2)); 
                and1smp = and1uap + and1bap;
                and2uap = quap*and2*(UAP/(Kuap + UAP))*(KIDO/(KIDO + O2))*(NO3/(KNO3 + NO3)); 
                and2bap = qbap*and2*(BAP/(Kbap + BAP))*(KIDO/(KIDO + O2))*(NO3/(KNO3 + NO3)); 
                and2smp = and2uap + and2bap;
                andsmp = and1smp +and2smp;
                
                %Mass-balance equations
                
                %Suspended heterotrophs
                
                XhsNew = Xhs + timeInc*(TFb*Yh*(1-kuaph-kepsh)*UTLbom*Xhs...
                    + TFb*Yh*(1-kuaph-kepsh)*and1bom*Xhs + TFb*Yh*(1-kuaph-kepsh)*and2bom*Xhs...
                    + TFb*Yp*UTLsmp*Xhs + TFb*Yp*and1smp*Xhs + TFb*Yp*and2smp*Xhs...
                    - bh*((O2/(KO2+O2))+(NO2/(KNO2+NO2))*(KIDO/(KIDO*O2))+(NO3/(KNO3+NO3))*(KIDO/(KIDO*O2)))*Xhs ...
                    - (kd1x1*HOCl + kd2x1*OCl + kd3x1*NH2Cl)*Xhs);
                
                if XhsNew  < 0.0
                    XhsNew = 0.0;
                end
                
                %Suspended ammonia oxidisers
                
                Xn1sNew = Xn1s + timeInc*(TFb2*Yn1*(1-kuapn1-kepsn1)*UTLNH3*Xn1s ...
                    - bn1*(O2/(KO2+O2))*Xn1s ...
                    - (kd1x3*HOCl + kd2x3*OCl + kd3x3*NH2Cl)*Xn1s);
                
                
                if Xn1sNew  < 0.0
                    Xn1sNew = 0.0;
                end
                
                %Suspended nitrite oxidisers
                
                Xn2sNew = Xn2s + timeInc*(TFb2*Yn2*(1-kuapn2-kepsn2)*UTLNO2*Xn2s...
                    - bn2*(O2/(KO2+O2))*Xn2s - (kd1x5*HOCl + kd2x5*OCl + kd3x5*NH2Cl)*Xn2s);
                
                if Xn2sNew  < 0.0
                    Xn2sNew = 0.0;
                end
                
                %Suspended inert biomass
                
                XisNew = Xis + timeInc*((1-fd)*bh*((O2/(KO2+O2))+(NO2/(KNO2+NO2))*(KIDO/(KIDO*O2))+(NO3/(KNO3+NO3))*(KIDO/(KIDO*O2)))*Xhs ...
                    + (1-fd)*bn1*(O2/(KO2+O2))*Xn1s + (1-fd)*bn2*(O2/(KO2+O2))*Xn2s...
                    + (1-fd)*(kd1x1*HOCl + kd2x1*OCl + kd3x1*NH2Cl)*Xhs...
                    + (1-fd)*(kd1x3*HOCl + kd2x3*OCl + kd3x3*NH2Cl)*Xn1s...
                    + (1-fd)*(kd1x5*HOCl + kd2x5*OCl + kd3x5*NH2Cl)*Xn2s...
                    - TFc*(kox1x2*HOCl + kox2x2*OCl + kox3x2*NH2Cl)*Xis);
                
                if XisNew  < 0.0
                    XisNew = 0.0;
                end
                
                %Suspended extracellular polymeric substances
                
                EPSsNew = EPSs + timeInc*(TFb*kepsh*UTLbom*Xhs ...
                    + TFb2*(kepsn1/MN)*UTLNH3*Xn1s + TFb2*(kepsn2/MN)*UTLNO2*Xn2s...
                    + TFb*and1bom*kepsh*Xhs ...
                    + TFb*and2bom*kepsh*Xhs ...
                    - khydEPS*EPSs - TFc*(kox3x2*NH2Cl*EPSs - kox1x2*HOCl*EPSs ...
                    - kox2x2*OCl*EPSs));
                
                if EPSsNew  < 0.0
                    EPSsNew = 0.0;
                end
                
                %Biodegradable organic matter
                
                BOM1New = BOM1 + timeInc*(-TFb*UTLbom1*Xhs...
                    - TFb*and1bom1*Xhs - TFb*and2bom1*Xhs...
                    - TFc*(kox1x1*HOCl + kox2x1*OCl + kox3x1*NH2Cl)*BOM1);
                
                %Individual BOM1 loss terms
                
                bom1a(iteration) = bom1a(iteration) - timeInc*(TFb*UTLbom1*Xhs);
                bom1b(iteration) = bom1b(iteration) - timeInc*(TFb*and1bom1*Xhs);
                bom1c(iteration) = bom1c(iteration) - timeInc*(TFb*and2bom1*Xhs);
                bom1d(iteration) = bom1d(iteration) - timeInc*(TFc*(kox1x1*HOCl + kox2x1*OCl + kox3x1*NH2Cl)*BOM1);
                
                
                if BOM1New  < 0.0
                    BOM1New = 0.0;
                end
                
                BOM2New = BOM2 + timeInc*(-TFb*UTLbom2*Xhs - TFb*and1bom2*Xhs - TFb*and2bom2*Xhs...
                    + fd*(kd1x1*HOCl + kd2x1*OCl + kd3x1*NH2Cl)*Xhs...
                    + fd*(kd1x3*HOCl + kd2x3*OCl + kd3x3*NH2Cl)*Xn1s...
                    + fd*(kd1x5*HOCl + kd2x5*OCl + kd3x5*NH2Cl)*Xn2s...
                    - TFc*(kox1x2*HOCl + kox2x2*OCl + kox3x2*NH2Cl)*BOM2);
                
                %Individual BOM2 loss and production terms
                
                bom2a(iteration) = bom2a(iteration) - timeInc*(TFb*UTLbom2*Xhs);
                bom2b(iteration) = bom2b(iteration) - timeInc*(TFb*and1bom2*Xhs);
                bom2c(iteration) = bom2c(iteration) - timeInc*(TFb*and2bom2*Xhs);
                bom2d(iteration) = bom2d(iteration) - timeInc*(TFc*(kox1x2*HOCl + kox2x2*OCl + kox3x2*NH2Cl)*BOM2);
                bom2e(iteration) = bom2e(iteration) + timeInc*(fd*(kd1x1*HOCl + kd2x1*OCl + kd3x1*NH2Cl)*Xhs ...
                    + fd*(kd1x3*HOCl + kd2x3*OCl + kd3x3*NH2Cl)*Xn1s ...
                    + fd*(kd1x5*HOCl + kd2x5*OCl + kd3x5*NH2Cl)*Xn2s);
                
                if BOM2New  < 0.0
                    BOM2New = 0.0;
                end
                
                %Soluble microbial products
                
                UAPNew = UAP + timeInc*(TFb*kuaph*UTLbom*Xhs ...
                    + TFb2*(kuapn1/MN)*UTLNH3*Xn1s...
                    + TFb2*(kuapn2/MN)*UTLNO2*Xn2s - TFb*UTLuap*Xhs...
                    + TFb*kuaph*and1bom1*Xhs + TFb*kuaph*and2bom1*Xhs...
                    + TFb*kuaph*and1bom2*Xhs + TFb*kuaph*and2bom2*Xhs...
                    - TFb*and1uap*Xhs - TFb*and2uap*Xhs...
                    - TFc*(kox1x2*HOCl + kox2x2*OCl + kox3x2*NH2Cl)*UAP);
                
                if UAPNew  < 0.0
                    UAPNew = 0.0;
                end
                
                BAPNew = BAP + timeInc*(khydEPS*EPSs...
                    - TFb*UTLbap*Xhs - TFb*and1bap*Xhs - TFb*and2bap*Xhs...
                    - TFc*(kox1x2*HOCl + kox2x2*OCl + kox3x2*NH2Cl)*BAP);
                
                if BAPNew  < 0.0
                    BAPNew = 0.0;
                end
                
                %Total ammonia
                
                CtNH3New = CtNH3 + timeInc*(-TFb2*UTLNH3*Xn1s... 
                    - TFb2*MN*Yn1*(1-kuapn1-kepsn1)*UTLNH3*Xn1s... 
                    - TFb2*MN*Yn2*(1-kuapn2-kepsn2)*UTLNO2*Xn2s... 
                    - TFb2*MN*(kuapn1/MN)*UTLNH3*Xn1s ... 
                    - TFb2*MN*(kuapn2/MN)*UTLNO2*Xn2s ... 
                    - TFb2*MN*(kepsn1/MN)*UTLNH3*Xn1s ... 
                    - TFb2*MN*(kepsn2/MN)*UTLNO2*Xn2s ... 
                    + MN*fd*(bh*((O2/(KO2+O2))+(NO2/(KNO2+NO2))*(KIDO/(KIDO*O2))+(NO3/(KNO3+NO3))*(KIDO/(KIDO*O2)))*Xhs+bn1*(O2/(KO2+O2))*Xn1s+bn2*(O2/(KO2+O2))*Xn2s)...
                    + TFc*kAD1*NH2Cl*NH2Cl - TFc*kAD2*NHCl2*NH3*H...
                    + TFc*kox5*NH2Cl*NO2...
                    + TFc*kAD3*(NH2Cl*kCle/NH3)*NH2Cl ...
                    + TFc*11*MN*(kox3x1*BOM1+kox3x2*BOM2...
                    + kox3x2*EPSs + kox3x2*BAP+kox3x2*UAP+kox3x2*Xis)*NH2Cl...
                    + TFb*MN*(1-Yp)*UTLbap*Xhs ...
                    + TFb*MN*(1-Yp)*and1bap*Xhs ...
                    + TFb*MN*(1-Yp)*and2bap*Xhs ...
                    + TFb*MN*(1-Yp)*UTLuap*Xhs ...
                    + TFb*MN*(1-Yp)*and1uap*Xhs ...
                    + TFb*MN*(1-Yp)*and2uap*Xhs ...
                    + TFb*MN*(1-Yh*(1-kuaph-kepsh)-(kuaph+kepsh))*UTLbom*Xhs...
                    + TFb*MN*(1-Yh*(1-kuaph-kepsh)-(kuaph+kepsh))*and1bom*Xhs...
                    + TFb*MN*(1-Yh*(1-kuaph-kepsh)-(kuaph+kepsh))*and2bom*Xhs);
                
                if CtNH3New  < 0.0
                    CtNH3New = 1E-18;
                end
                
                %Nitrite 
                
                NO2New = NO2 + timeInc*(-TFb2*UTLNO2*Xn2s...
                    - TFb*MN*(1-kuaph-kepsh-Yh*(1-kuaph-kepsh))*and1bom*Xhs - TFb*MN*(1-Yp)*and1smp*Xhs... 
                    + TFb*MN*(1-kuaph-kepsh-Yh*(1-kuaph-kepsh))*and2bom*Xhs + TFb*MN*(1-Yp)*and2smp*Xhs ... 
                    -MN*fd*bh*(NO2/(KNO2+NO2))*(KIDO/(KIDO*O2))*Xhs ...
                    +MN*fd*bh*(NO3/(KNO3+NO3))*(KIDO/(KIDO*O2))*Xhs ...
                    + TFb2*UTLNH3*Xn1s ... 
                    + TFb2*UTLNHOHCl*Xn1s ...
                    - TFc*kox5*NH2Cl*NO2);
                
                if NO2New  < 0.0
                    NO2New = 0.0;
                end
                
                %Nitrate 
                
                NO3New = NO3 + timeInc*(TFb2*UTLNO2*Xn2s... 
                    - TFb*MN*(1-kuaph-kepsh-Yh*(1-kuaph-kepsh))*and2bom*Xhs - TFb*MN*(1-Yp)*and2smp*Xhs... 
                    -MN*fd*bh*(NO3/(KNO3+NO3))*(KIDO/(KIDO*O2))*Xhs ...
                    + TFc*kox5*NH2Cl*NO2);
                
                if NO3New  < 0.0
                    NO3New = 0.0;
                end
                
                %Nitrogen gas 
                
                N2New = N2 + timeInc*(TFb*MN*0.5*(1-kuaph-kepsh-Yh*(1-kuaph-kepsh))*and1bom*Xhs ...
                    + TFb*MN*0.5*(1-Yp)*and1smp*Xhs ...
                    + TFc*kAD5*NHOHCl ...
                    +MN*0.5*fd*bh*(NO2/(KNO2+NO2))*(KIDO/(KIDO*O2))*Xhs);
                   
                
                if N2New  < 0.0
                    N2New = 0.0;
                end
                
                %Free chlorine (HOCl + OCl)
                
                CtOClNew = CtOCl + timeInc*(-TFc*2.0*MC*(kox1x1*BOM1 + kox1x2*BOM2 ...
                    + kox1x2*BAP + kox1x2*EPSs + kox1x2*UAP + kox1x2*Xis)*HOCl ...
                    - TFc*2.0*MC*(kox2x1*BOM1 + kox2x2*BOM2 + kox2x2*BAP ...
                    + kox2x2*EPSs + kox2x2*UAP + kox2x2*Xis)*OCl - TFc*kox4*HOCl*NO2);
                
                if CtOClNew  < 0.0
                    CtOClNew = 0;
                end
                
                %Monochloramine 
                
                NH2ClNew = NH2Cl + timeInc*(-TFc*2.0*kAD3*(NH2Cl*kCle/NH3)*NH2Cl ...
                    - TFc*2.0*kAD1*NH2Cl*NH2Cl + TFc*2.0*kAD2*NHCl2*NH3*H...
                    - TFc*kAD5*NHOHCl - TFc*kox5*NH2Cl*NO2 - TFb2*UTLNH2Cl*Xn1s...
                    - TFc*2.0*MC*(kox3x1*BOM1 ...
                    + kox3x2*BOM2 +kox3x2*BAP...
                    + kox3x2*UAP+kox3x2*Xis + kox3x2*EPSs)*NH2Cl);
                
                %Individual monochloramine loss terms
                
                ad1(iteration) = ad1(iteration) - timeInc*(2.0*TFc*kAD1*NH2Cl*NH2Cl);
                ad2(iteration) =  ad2(iteration) + timeInc*(2.0*TFc*kAD2*NHCl2*NH3*H);
                ad3(iteration) = ad3(iteration) - timeInc*(2.0*TFc*kAD3*(NH2Cl*kCle/NH3)*NH2Cl);
                ad5(iteration) = ad5(iteration) - timeInc*(TFc*kAD5*NHOHCl);
                ox5(iteration) = ox5(iteration) - timeInc*(TFc*kox5*NH2Cl*NO2);
                comet(iteration) = comet(iteration) - timeInc*(TFb2*UTLNH2Cl*Xn1s);
                organicox(iteration) = organicox(iteration) - timeInc*(TFc*2.0*MC*(kox3x1*BOM1 ...
                    + kox3x2*BOM2 +kox3x2*BAP...
                    + kox3x2*UAP+kox3x2*Xis + kox3x2*EPSs)*NH2Cl);
                
                if NH2ClNew  < 0.0
                    NH2ClNew = 0.0;
                end
                
                if NH2ClNew < 1E-12 %If NH2Cl approaches 0, Reaction AD5 cannot occur. NHOHCl + NH2Cl -> N2 + H20 + 2HCl
                    kAD5 = 0;
                else
                    kAD5 = 5.0E-3;
                end
                
                %Dichloramine 
                            
                NHCl2New = NHCl2 + timeInc*(TFc*kAD1*NH2Cl*NH2Cl...
                    - TFc*kAD2*NHCl2*NH3*H + TFc*kAD3*(NH2Cl*kCle/NH3)*NH2Cl ...
                    - TFc*kAD4*NHCl2*OH);
                
                if NHCl2New  < 0.0
                    NHCl2New = 0.0;
                end
                
                
                %Chlorohydroxylamine 
                
                NHOHClNew = NHOHCl + timeInc*(TFb2*UTLNH2Cl*Xn1s...
                    - TFb2*UTLNHOHCl*Xn1s + TFc*kAD4*NHCl2*OH - TFc*kAD5*NHOHCl);
                
                if NHOHClNew  < 0.0
                    NHOHClNew = 0.0;
                end
                
                %Chloride 
                
                ClNew = Cl + timeInc*(TFc*kAD4*NHCl2*OH + TFc*2.0*kAD5*NHOHCl...
                    + TFb2*UTLNHOHCl*Xn1s + TFc*kox4*HOCl*NO2 + TFc*kox5*NH2Cl*NO2...
                    + TFc*2.0*MC*(kox1x1*BOM1 + kox1x2*BOM2...
                    + kox1x2*BAP + kox1x2*UAP + kox1x2*Xis + kox1x2*EPSs)*HOCl...
                    + TFc*2.0*MC*(kox2x1*BOM1 + kox2x2*BOM2 + kox2x2*BAP...
                    + kox2x2*UAP + kox2x2*Xis + kox2x2*EPSs)*OCl...
                    + TFc*2.0*MC*(kox3x1*BOM1 + kox3x2*BOM2 + kox3x2*BAP...
                    + kox3x2*UAP + kox3x2*Xis + kox3x2*EPSs)*NH2Cl);
                
                if ClNew < 0.0
                    ClNew = 0.0;
                end
                
                %Carbon dioxide
                
                CO2New = CO2 + ...
                    timeInc*(TFb*(1-Yh*(1-kuaph-kepsh)-(kuaph+kepsh))*UTLbom*Xhs...
                    +TFb*(1-Yh*(1-kuaph-kepsh)-(kuaph+kepsh))*and1bom*Xhs...
                    +TFb*(1-Yh*(1-kuaph-kepsh)-(kuaph+kepsh))*and2bom*Xhs...
                    + TFb*(1-Yp)*and1smp*Xhs + TFb*(1-Yp)*and2smp*Xhs...
                    + TFb*(1-Yp)*UTLsmp*Xhs - TFb2*(Yn1*(1-kuapn1-kepsn1)+((kuapn1/MN)+(kepsn1/MN)))*UTLNH3*Xn1s...
                    -  TFb2*(Yn2*(1-kuapn2-kepsn2)+((kuapn2/MN)+(kepsn2/MN)))*UTLNO2*Xn2s...
                    + fd*(bh*((O2/(KO2+O2))+(NO2/(KNO2+NO2))*(KIDO/(KIDO*O2))+(NO3/(KNO3+NO3))*(KIDO/(KIDO*O2)))*Xhs+bn1*(O2/(KO2+O2))*Xn1s+bn2*(O2/(KO2+O2))*Xn2s) + TFc*(kox3x1*BOM1...
                    + kox3x2*BOM2+ kox3x2*BAP + kox3x2*EPSs + kox3x2*UAP+kox3x2*Xis)*NH2Cl...
                    + TFc*(kox1x1*BOM1 + kox1x2*BOM2 + kox1x2*BAP + kox1x2*EPSs + kox1x2*UAP + kox1x2*Xis)*HOCl...
                    + TFc*(kox2x1*BOM1 + kox2x2*BOM2 + kox2x2*BAP + kox2x2*EPSs + kox2x2*UAP + kox2x2*Xis)*OCl);
                
                if CO2New < 0.0
                    CO2New = 0.0;
                end
                
                %Oxygen
                
                O2New = O2 + timeInc*((-TFb*(1-kuaph-kepsh-Yh*(1-kuaph-kepsh))*UTLbom*Xhs...
                    - TFb*(1-Yp)*UTLsmp*Xhs - TFb2*(1-kuapn1-kepsn1-Yn1*(1-kuapn1-kepsn1))*UTLNH3*Xn1s...
                    - TFb2*(1-kuapn2-kepsn2-Yn2*(1-kuapn2-kepsn2))*UTLNO2*Xn2s...
                    - fd*(bh*(O2/(KO2+O2))*Xhs+bn1*(O2/(KO2+O2))*Xn1s+bn2*(O2/(KO2+O2))*Xn2s))+kla*(O2s-O2));
                
                if O2New < 0.0
                    O2New = 0.0;
                elseif O2New > O2s
                    O2New = O2s;
                end
                
                %Stability factor AOB
                
                SFxn1New = TFb2*Yn1*(1 - kuapn1 - kepsn1)*UTLNH3 - bn1*(O2/(KO2+O2)) - kd3x3*NH2Cl;
                
                %Stability factor NOB
                
                SFxn2New = TFb2*Yn2*(1 - kuapn2 - kepsn2)*UTLNO2 - bn2*(O2/(KO2+O2)) - kd3x5*NH2Cl;
                
                %Stability factor Xhs
                
                SFxhNew = TFb*Yh*(1-kuaph-kepsh)*UTLbom ...
                    + TFb*Yh*(1-kuaph-kepsh)*and1bom + TFb*Yh*(1-kuaph-kepsh)*and2bom ...
                    + TFb*Yp*UTLsmp + TFb*Yp*and1smp + TFb*Yp*and2smp...
                    - bh*((O2/(KO2+O2))+(NO2/(KNO2+NO2))*(KIDO/(KIDO*O2))+(NO3/(KNO3+NO3))*(KIDO/(KIDO*O2))) - kd3x1*NH2Cl;
                
                %Update concentration variables for next iteration
                
                Xhs = XhsNew;
                Xn1s = Xn1sNew;
                Xn2s = Xn2sNew;
                Xis = XisNew;
                EPSs = EPSsNew;
                BOM1 = BOM1New;
                BOM2 = BOM2New;
                UAP = UAPNew;
                BAP = BAPNew;
                CtNH3 = CtNH3New;
                NH2Cl = NH2ClNew;
                NHCl2 = NHCl2New;
                NHOHCl = NHOHClNew;
                Cl = ClNew;
                NO2 = NO2New;
                NO3 = NO3New;
                N2 = N2New;
                CtOCl = CtOClNew;
                CO2 = CO2New;
                O2 = O2New;
                SFxn1 = SFxn1New;
                SFxn2 = SFxn2New;
                SFxh = SFxhNew;
                
                %Save initial variable values to output table
                
                if day == 1 && hour == 1 && t == 1
                    output(1, :) = [0 HOCl OCl NH2Cl NHCl2 NHOHCl Cl BOM1 BOM2 Xhs NH4 NH3 CtNH3 NO2 NO3 N2 Xn1s Xn2s EPSs Xis CO2 UAP BAP O2 OH H2CO3 HCO3 kox5 kAD1 CtOCl UTLbom UTLsmp UTLNH3 UTLNO2 kAD2 H SFxn1 SFxn2 SFxh];
                end
                
            end
        end
        
        %Save daily variable values to output table
        
        output(day+1, :) = [day HOCl OCl NH2Cl NHCl2 NHOHCl Cl BOM1 BOM2 Xhs NH4 NH3 CtNH3 NO2 NO3 N2 Xn1s Xn2s EPSs Xis CO2 UAP BAP O2 OH H2CO3 HCO3 kox5 kAD1 CtOCl UTLbom UTLsmp UTLNH3 UTLNO2 kAD2 H SFxn1 SFxn2 SFxh];
    end
    
    %Assign new variable names to values in output table
    
    daynumber(:,iteration) = output (:,1);
    hocl(:,iteration) = output (:,2);
    ocl(:,iteration) = output (:,3);
    monochloramine(:,iteration) = output (:,4);
    dichloramine(:,iteration) = output (:,5);
    nhohcl(:,iteration) = output (:,6);
    chloride(:,iteration) = output (:,7);
    bom1(:,iteration) = output (:,8);
    bom2(:,iteration) = output (:,9);
    heterotrophs(:,iteration) = output (:,10);
    ammonium(:,iteration) = output (:,11);
    ammonia(:,iteration) = output (:,12);
    totalnh3(:,iteration) = output (:,13);
    nitrite(:,iteration) = output (:,14);
    nitrate(:,iteration) = output (:,15);
    nitrogen(:,iteration) = output (:,16);
    aob(:,iteration) = output (:,17);
    nob(:,iteration) = output (:,18);
    eps(:,iteration) = output (:,19);
    inertbacteria(:,iteration) = output (:,20);
    carbondioxide(:,iteration) = output (:,21);
    uap(:,iteration) = output (:,22);
    bap(:,iteration) = output (:,23);
    oxygen(:,iteration) = output (:,24);
    oh(:,iteration) = output (:,25);
    h2co3(:,iteration) = output (:,26);
    hco3(:,iteration) = output (:,27);
    Kox5(:,iteration) = output (:,28);
    kad1(:,iteration) = output (:,29);
    ctocl(:,iteration) = output (:,30);
    utlbom(:,iteration) = output (:,31);
    utlsmp(:,iteration) = output (:,32);
    utlnh3(:,iteration) = output (:,33);
    utlno2(:,iteration) = output (:,34);
    kad2(:,iteration) = output (:,35);
    h(:,iteration) = output (:,36);
    sfxn1(:,iteration) = output (:,37);
    sfxn2(:,iteration) = output (:,38);
    sfxh(:,iteration) = output (:,39);
    
    %Determine total mass of Cl, N and C
    
    totalCl(:,iteration) = hocl(:,iteration) + ocl(:,iteration) + monochloramine(:,iteration) + 2*dichloramine(:,iteration) + nhohcl(:,iteration) + chloride(:,iteration);
    totalN(:,iteration) = monochloramine(:,iteration) + dichloramine(:,iteration) + nhohcl(:,iteration) + totalnh3(:,iteration) + nitrite(:,iteration) + nitrate(:,iteration) + 2*nitrogen(:,iteration) + MN*(bom1(:,iteration) + bom2(:,iteration) + heterotrophs(:,iteration) + aob(:,iteration) + nob(:,iteration) + eps(:,iteration) + inertbacteria(:,iteration) + uap(:,iteration) + bap(:,iteration));
    totalC(:,iteration) = MC*(bom1(:,iteration) + bom2(:,iteration) + heterotrophs(:,iteration) + aob(:,iteration) + nob(:,iteration) + eps(:,iteration) + inertbacteria(:,iteration) + uap(:,iteration) + bap(:,iteration) + carbondioxide(:,iteration));
    
    %Determine total change in mass of Cl, N and C
    
    totalChangeCl = (totalCl(days,iteration) - totalCl(1,iteration))/totalCl(days,iteration)*100; 
    totalChangeN = (totalN(days,iteration) - totalN(1,iteration))/totalN(days,iteration)*100; 
    totalChangeC = (totalC(days,iteration) - totalC(1,iteration))/totalC(days,iteration)*100; 

end

%If user has opted to plot graphs

if strcmpi(promptButton, 'Yes')
    
    lines = cellstr(char('-','--',':','-.','-','--',':','-.','-','--',':','-.'));
    
    % Graph unit conversions
    
    monochloramineP = monochloramine.*70909.1;
    ammoniaP = ammonia.*14000;
    bom1P = bom1./2670;
    bom2P = bom2./2670;
    heterotrophsP = heterotrophs.*2404;
    aobP = aob.*2404;
    nobP = nob.*2404;
    nitriteP = nitrite.*14000;
    nitrateP = nitrate.*14000;
    uapP = uap./2670;
    bapP = bap./2670;
    epsP = eps./2670;
    
    %Plot output graphs
    
    figure (1);
    plot (daynumber, monochloramineP);
    grid;
    xlabel ('Day');
    ylabel ('Monochloramine (mg/l as Cl_2)');
    if not(strcmp(iterationString,''))
        legend(iterationString);
    end
    
    figure (2);
    plot (daynumber, bom1P);
    grid;
    xlabel ('Day');
    ylabel ('BOM_1 (mg/l as C)');
    if not(strcmp(iterationString,''))
        legend(iterationString);
    end
    
    figure (3);
    plot (daynumber, bom2P);
    grid;
    xlabel ('Day');
    ylabel ('BOM_2 (mg/l as C)');
    if not(strcmp(iterationString,''))
        legend(iterationString);
    end
    
    figure (4);
    plot (daynumber, ammoniaP);
    grid;
    xlabel ('Day');
    ylabel ('Ammonia (mg/l as N)');
    if not(strcmp(iterationString,''))
        legend(iterationString);
    end
    
    figure (5);
    plot (daynumber, nitriteP);
    grid;
    xlabel ('Day');
    ylabel ('Nitrite (mg/l as N)');
    if not(strcmp(iterationString,''))
        legend(iterationString);
    end
    
    figure (6);
    plot (daynumber, nitrateP);
    grid;
    xlabel ('Day');
    ylabel ('Nitrate (mg/l as N)');
    if not(strcmp(iterationString,''))
        legend(iterationString);
    end
    
    figure(7);
    [f7,line1,line2] = plotyy(daynumber,heterotrophsP,daynumber,sfxh);
    grid on;
    set(f7,'XGrid','on','YGrid','off');
    xlabel ('Day');
    ylabel(f7(1),'Heterotrophs (cells/ml)')
    ylabel(f7(2),'Heterotroph Stability Factor (1/hour)');
    set(line1,'color','b');
    set(line2,'color',[0 .5 0]);
    if not(strcmp(iterationString,''))
        h = legend(iterationString);
        leg_line = findobj(h,'type','Line');
        for i = 1:length(leg_line)
            set(leg_line(i),'Color','black')
        end
        set(f7,{'ycolor'},{'b';[0 .5 0]});
        for i = 1:(length(leg_line)/2)
            set([line1(i);line2(i)],'LineStyle',lines{i});
        end
    end
    
    figure(8);
    [f8,line1,line2] = plotyy(daynumber,aobP,daynumber,sfxn1);
    grid on;
    set(f8,'XGrid','on','YGrid','off');
    xlabel ('Day');
    ylabel(f8(1),'AOB (cells/ml)')
    ylabel(f8(2),'AOB Stability Factor (1/hour)');
    set(line1,'color','b');
    set(line2,'color',[0 .5 0]);
    if not(strcmp(iterationString,''))
        h = legend(iterationString);
        leg_line = findobj(h,'type','Line');
        for i = 1:length(leg_line)
            set(leg_line(i),'Color','black')
        end
        set(f8,{'ycolor'},{'b';[0 .5 0]});
        for i = 1:(length(leg_line)/2)
            set([line1(i);line2(i)],'LineStyle',lines{i});
        end
    end
    
    figure(9);
    [f9,line1,line2] = plotyy(daynumber,nobP,daynumber,sfxn2);
    grid on;
    set(f9,'XGrid','on','YGrid','off');
    xlabel ('Day');
    ylabel(f9(1),'NOB (cells/ml)')
    ylabel(f9(2),'NOB Stability Factor (1/hour)');
    set(line1,'color','b');
    set(line2,'color',[0 .5 0]);
    if not(strcmp(iterationString,''))
        h = legend(iterationString);
        leg_line = findobj(h,'type','Line');
        for i = 1:length(leg_line)
            set(leg_line(i),'Color','black')
        end
        set(f9,{'ycolor'},{'b';[0 .5 0]});
        for i = 1:(length(leg_line)/2)
            set([line1(i);line2(i)],'LineStyle',lines{i});
        end
    end
    
    figure (10);
    plot (daynumber, uapP);
    grid;
    xlabel ('day');
    ylabel ('UAP (mg/l as C)');
    if not(strcmp(iterationString,''))
        legend(iterationString);
    end
    
    figure (11);
    plot (daynumber, bapP);
    grid;
    xlabel ('day');
    ylabel ('BAP (mg/l as C)');
    if not(strcmp(iterationString,''))
        legend(iterationString);
    end
    
    figure (12);
    plot (daynumber, epsP);
    grid;
    xlabel ('day');
    ylabel ('EPS (mg/l as C)');
    if not(strcmp(iterationString,''))
        legend(iterationString);
    end
    
end