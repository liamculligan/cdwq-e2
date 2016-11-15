%Expanded Comprehensive Disinfection and Water Quality Model (CDWQ-E2)

%Originally written by: John Woolschlager
%Northwestern University, Evanston, USA
%Expanded and modified by: Precious Biyela
%Arizona State University, Tempe, Arizona, USA
%Expanded and modified by: Liam Culligan
%University of the Witwatersrand, Johannesburg, South Africa
%November 2014

clear, clc
format short e

%Set global variables such that advection calculations can be performed in
%a separate function

global elementReservoir
global elementReservoirC1
global elementReservoirC2
global elementReservoirC3
global timeInc
global Qin
global elementPipes
global hour
global elementVolumes
global cElementStartNodes
global elementV
global Qout
global reservoirCompartments

%Set run time

days = 60;

%Create prompt for user to select input file

[file, pathname] = uigetfile('*.xlsx', 'Select a distribution system input file');
filename = strcat(pathname,file);

%Import Pipe Properties sheet from input file

properties = xlsread(filename,'Pipe Properties');
pipes = properties(:,1);
startNodes = properties(:,2);
endNodes = properties(:,3);
diameters = properties(:,4);
lengths = properties(:,5);
C = properties(:,6);
reservoir = properties(:,7);
[~, materials, ~] = xlsread(filename,'Pipe Properties','H:H');
%Remove pipe materials header
materials(1) = [];

totalPipes = numel(pipes);

%Import Pipe Flows sheet from input file

flows = xlsread(filename,'Flows');
maxFlows = max(abs(flows),[],2); %Find max value within each pipe

%Calculate pipe areas and maximum velocities

PI = 3.1416;
areas = (PI.*(diameters/2).*(diameters/2));
maxVelocities = maxFlows./areas;


%Find maximum time increment based on maximum velocities

maxTimeInc = floor(min(lengths./maxVelocities)*100)/100;

%Prompt user to set time increment

timeIncMessage = sprintf('The recommended time step is %d. Use this time step?', maxTimeInc);
timeIncQuery = questdlg(timeIncMessage, 'Use recommended time step?', 'Yes', 'No', 'Yes');
drawnow;

if strcmpi(timeIncQuery, 'Yes')
    timeInc = maxTimeInc;
else
    timeInc = input('Enter your selected time step: ');
end

%Set time steps based on selected time increment

timeSteps = (1/timeInc);

%Calculate total number of elements within each pipe

minElementLengths = maxVelocities.*timeInc;
numElements = floor(lengths./minElementLengths);

numElements(numElements==0) = 1;

%Prompt user to select reservoir compartment model

reservoirMessage = sprintf('Please select a reservoir compartment model to use:');
reservoirModel = questdlg(reservoirMessage, 'Compartment Model', 'One-Compartment', 'Two-Compartment', 'Three-Compartment', 'One-Compartment');
drawnow;

%Set reservoirs to selected number of compartments

if strcmpi(reservoirModel, 'One-Compartment')
    numElements(reservoir==1) = 1;
    reservoirCompartments = 1;
elseif strcmpi(reservoirModel, 'Two-Compartment')
    numElements(reservoir==1) = 2;
    reservoirCompartments = 2;
else
    numElements(reservoir==1) = 3;
    reservoirCompartments = 3;
end

%Set conversion factors

moleN = 14E6;
MC = 3.1211E-8;
MN = 6.2422E-9;

%Set equilibrium parameters and biological parameters (per hour)

bh = 0.1/24.0;
bn1 = 0.05/24.0;
bn2 = 0.05/24.0;
bxi = (((bh + bn1 + bn2)/3)/24.0);
fd = 0.8;
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
KNH4 = 7.14E-5;
KNCl = 1000/moleN;
KNO2 = 5.36E-5;
KNO3 = 1000/moleN;
KNH2Cl = 2000.0/moleN;
KNHOHCl = 1000.0/moleN;
KIDO = 0.05*1000;
Kbap = 30000.0;
Kuap = 20000.0;
KO2 = 1000*0.20;
KOA = 0.40;
KNH = 1.00;
KNO = 0.50;
qh = 10/24;
qn1 = 1.6/(24*moleN);
qn2 = 7.0/(24*moleN);
qnCl = 1.6/(24*moleN);
qbap = 2.0/24;
quap = 13.0/24;
qNH4 = 1.6/(24.0*moleN);
qNH2Cl = 1.6/(24.0*moleN);
qNHOHCl = 1.6/(24.0*moleN);
qNO2 = 7.0/(24.0*moleN);
qbom1 = 10/24;
qbom2 = 10/24;
and1 = 0.6;
and2 = 0.6;
Yh = 0.6;
Yn1 = 6.16E6;
Yn2 = 1.68E6;
Yp = 0.6;

%Set chemical conditions (per hour)

kH = 2.5E7;
kH2CO3 = 4.0E4;
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
KAD2 = 2.16E8;
KAD3 = 1.0E6;
KAD4 = 6.0E5;

%Set physical and temperature parameters

kdet = 0.013;
kads = 10.0;
AS = 0.58;
dens = 999.0;
g = 9.81;
R = 8.314;
BFden = 4.16E11;
kpro1 = 1;
kpro2 = 10;
Tbio = 1.05;
Tbio2 = 1.00;
Tchem = 1.05;

%Import Initial Node Concentrations sheet from input file

nConcentrations = xlsread(filename,'Initial Node Concentrations');
iNodeNums = nConcentrations(:,1);
iCtOCl = nConcentrations(:,2);
iNH2Cl = nConcentrations(:,3);
iNHCl2 = nConcentrations(:,4);
iNHOHCl = nConcentrations(:,5);
iCl = nConcentrations(:,6);
iBOM1 = nConcentrations(:,7);
iBOM2 = nConcentrations(:,8);
iXis = nConcentrations(:,9);
iXhs = nConcentrations(:,10);
iCtNH3 = nConcentrations(:,11);
iNO2 = nConcentrations(:,12);
iNO3 = nConcentrations(:,13);
iN2 = nConcentrations(:,14);
iXn1s = nConcentrations(:,15);
iXn2s = nConcentrations(:,16);
iCO2 = nConcentrations(:,17);
iUAP = nConcentrations(:,18);
iBAP = nConcentrations(:,19);
iO2 = nConcentrations(:,20);
iEPSs = nConcentrations(:,21);
iTime = nConcentrations(:,22);

%Create nodes array
%Merge all imported node arrays into nodes array, and apply
%Unique command to remove duplicates and sort from smallest to largest

nodes = unique([startNodes; endNodes; iNodeNums]);

totalNodes = numel(nodes);

nCtOCl = zeros(totalNodes,1);
nNH2Cl = zeros(totalNodes,1);
nNHCl2 = zeros(totalNodes,1);
nNHOHCl = zeros(totalNodes,1);
nCl = zeros(totalNodes,1);
nBOM1 = zeros(totalNodes,1);
nBOM2 = zeros(totalNodes,1);
nXis = zeros(totalNodes,1);
nXhs = zeros(totalNodes,1);
nCtNH3 = zeros(totalNodes,1);
nNO2 = zeros(totalNodes,1);
nNO3 = zeros(totalNodes,1);
nN2 = zeros(totalNodes,1);
nXn1s = zeros(totalNodes,1);
nXn2s = zeros(totalNodes,1);
nCO2 = zeros(totalNodes,1);
nUAP = zeros(totalNodes,1);
nBAP = zeros(totalNodes,1);
nO2 = zeros(totalNodes,1);
nEPSs = zeros(totalNodes,1);
nTime = zeros(totalNodes,1);

%Assign imported concentrations to relevant nodes

for node=1:numel(iNodeNums)
    nodeNum = find(nodes == iNodeNums(node));
    nCtOCl(nodeNum,1) = iCtOCl(node);
    nNH2Cl(nodeNum,1) = iNH2Cl(node);
    nNHCl2(nodeNum,1) = iNHCl2(node);
    nNHOHCl(nodeNum,1) = iNHOHCl(node);
    nCl(nodeNum,1) = iCl(node);
    nBOM1(nodeNum,1) = iBOM1(node);
    nBOM2(nodeNum,1) = iBOM2(node);
    nXis(nodeNum,1) = iXis(node);
    nXhs(nodeNum,1) = iXhs(node);
    nCtNH3(nodeNum,1) = iCtNH3(node);
    nNO2(nodeNum,1) = iNO2(node);
    nNO3(nodeNum,1) = iNO3(node);
    nN2(nodeNum,1) = iN2(node);
    nXn1s(nodeNum,1) = iXn1s(node);
    nXn2s(nodeNum,1) = iXn2s(node);
    nCO2(nodeNum,1) = iCO2(node);
    nUAP(nodeNum,1) = iUAP(node);
    nBAP(nodeNum,1) = iBAP(node);
    nO2(nodeNum,1) = iO2(node);
    nEPSs(nodeNum,1) = iEPSs(node);
    nTime(nodeNum,1) = iTime(node);
end

%Import Node Coordinates sheet from input file

nCoordinates = xlsread(filename,'Node Coordinates');
nCoordinatesNodeNums = nCoordinates(:,1);
iX = nCoordinates(:,2);
iY = nCoordinates(:,3);

%Assign imported coordinates to relevant nodes

for node=1:numel(nCoordinatesNodeNums)
    nodeNum = find(nodes == nCoordinatesNodeNums(node));
    nX(nodeNum,1) = iX(node);
    nY(nodeNum,1) = iY(node);
end

%Import Treatment Plant Node Properties sheet from input file

tPlant = xlsread(filename,'Treatment Plant Conditions');
nodeNum = tPlant(:,1);
TPnode = find(nodes == nodeNum); %Find position of treatment plant node number within nodes array

%Set treatment plant node concentrations

aCtOCl(TPnode,1) = tPlant(:,2);
aNH2Cl(TPnode,1) = tPlant(:,3);
aNHCl2(TPnode,1) = tPlant(:,4);
aNHOHCl(TPnode,1) = tPlant(:,5);
aCl(TPnode,1) = tPlant(:,6);
aBOM1(TPnode,1) = tPlant(:,7);
aBOM2(TPnode,1) = tPlant(:,8);
aXis(TPnode,1) = tPlant(:,9);
aXhs(TPnode,1) = tPlant(:,10);
aCtNH3(TPnode,1) = tPlant(:,11);
aNO2(TPnode,1) = tPlant(:,12);
aNO3(TPnode,1) = tPlant(:,13);
aN2(TPnode,1) = tPlant(:,14);
aXn1s(TPnode,1) = tPlant(:,15);
aXn2s(TPnode,1) = tPlant(:,16);
aCO2(TPnode,1) = tPlant(:,17);
aUAP(TPnode,1) = tPlant(:,18);
aBAP(TPnode,1) = tPlant(:,19);
aO2(TPnode,1) = tPlant(:,20);
aEPSs(TPnode,1) = tPlant(:,21);
aTime(TPnode,1) = tPlant(:,22);
Temp = tPlant(:,23); %Global system temperature
logH = tPlant(:,24); %Global system pH

%Import Initial Biofilm Concentrations sheet from input file

biofilm = xlsread(filename,'Initial Biofilm Concentrations');
pXhf = biofilm(:,2);
pXn1f = biofilm(:,3);
pXn2f = biofilm(:,4);
pEPSf = biofilm(:,5);
pXif = biofilm(:,6);

%Determine oxygen saturation value

kla = 1/24;
O2s = (14.59 - 0.3955*Temp + 0.0072*(Temp^2) - 0.0000619*(Temp^3))*1000;

CtCO3 = 5E-3;

%Determine equilibrium constants

kCle = 0.249*exp(-7290/(Temp + 273.15));
KeNH3 = 0.7029*exp(-52210/(R*(Temp + 273.15)));
KeHCO3 = (1.119E-5)*exp(-7700/(R*(Temp + 273.15)));
KeCO3 = (2.044E-8)*exp(-14900/(R*(Temp + 273.15)));
KeOCl = (8.273E-6)*exp(-13800/(R*(Temp + 273.15)));

%Calculate element properties

pElementLengths = lengths./numElements;
pElementVolumes = areas.*pElementLengths;

%Calculate total number of pipe elements and populate array

totalElements = sum(numElements);
elements = (1:totalElements)';

%Preallocate element variables

elementPipes = zeros(totalElements,1);
elementStartNodes = zeros(totalElements,1);
elementEndNodes = zeros(totalElements,1);
elementLengths = zeros(totalElements,1);
elementVolumes = zeros(totalElements,1);
elementVolumesNew = zeros(totalElements,1);
elementVolumesRatio = zeros(totalElements,1);
elementReservoir = zeros(totalElements,1);
elementReservoirC1 = zeros(totalElements,1);
elementReservoirC2 = zeros(totalElements,1);
elementReservoirC3 = zeros(totalElements,1);
Qin = zeros(totalPipes,24);
Qout = zeros(totalPipes,24);
cElementStartNodes = zeros(totalElements,24);
cElementEndNodes = zeros(totalElements,24);

% Preallocate substance variables

Xhf = zeros(totalElements,1);
Xn1f = zeros(totalElements,1);
Xn2f = zeros(totalElements,1);
EPSf = zeros(totalElements,1);
Xif = zeros(totalElements,1);
CtOCl = zeros(totalElements,1);
HOCl = zeros(totalElements,1);
OCl = zeros(totalElements,1);
NH2Cl = zeros(totalElements,1);
NHCl2 = zeros(totalElements,1);
NHOHCl = zeros(totalElements,1);
Cl = zeros(totalElements,1);
BOM1 = zeros(totalElements,1);
BOM2 = zeros(totalElements,1);
Xis = zeros(totalElements,1);
Xhs = zeros(totalElements,1);
CtNH3 = zeros(totalElements,1);
NO2 = zeros(totalElements,1);
NO3 = zeros(totalElements,1);
N2 = zeros(totalElements,1);
Xn1s = zeros(totalElements,1);
Xn2s = zeros(totalElements,1);
CO2 = zeros(totalElements,1);
UAP = zeros(totalElements,1);
BAP = zeros(totalElements,1);
O2 = zeros(totalElements,1);
EPSs = zeros(totalElements,1);
Time = zeros(totalElements,1);
EPS = zeros(totalElements,1);
Xh = zeros(totalElements,1);
Xi = zeros(totalElements,1);
Xn1 = zeros(totalElements,1);
Xn2 = zeros(totalElements,1);
NH4 = zeros(totalElements,1);
NH3 = zeros(totalElements,1);
kcr1 = zeros(totalElements,1);
kcr2 = zeros(totalElements,1);
kcr3 = zeros(totalElements,1);
kSC1 = zeros(totalElements,1);
kSC2 = zeros(totalElements,1);
kSC3 = zeros(totalElements,1);
kox4 = zeros(totalElements,1);
kcor11 = zeros(totalElements,1);
kcor12 = zeros(totalElements,1);
kcor13 = zeros(totalElements,1);
kcor21 = zeros(totalElements,1);
kcor22 = zeros(totalElements,1);
kcor23 = zeros(totalElements,1);
kcor31 = zeros(totalElements,1);
kcor32 = zeros(totalElements,1);
kcor33 = zeros(totalElements,1);
UTLbom1 = zeros(totalElements,1);
UTLbom2 = zeros(totalElements,1);
UTLuap = zeros(totalElements,1);
UTLbap = zeros(totalElements,1);
UTLNH3 = zeros(totalElements,1);
UTLNH2Cl = zeros(totalElements,1);
UTLNHOHCl = zeros(totalElements,1);
UTLNO2 = zeros(totalElements,1);
UTLbom = zeros(totalElements,1);
UTLsmp = zeros(totalElements,1);
and1bom1 = zeros(totalElements,1);
and1bom2 = zeros(totalElements,1);
and1bom = zeros(totalElements,1);
and2bom1 = zeros(totalElements,1);
and2bom2 = zeros(totalElements,1);
and2bom = zeros(totalElements,1);
andbom = zeros(totalElements,1);
and1uap = zeros(totalElements,1);
and1bap = zeros(totalElements,1);
and1smp = zeros(totalElements,1);
and2uap = zeros(totalElements,1);
and2bap = zeros(totalElements,1);
and2smp = zeros(totalElements,1);
andsmp = zeros(totalElements,1);
Lf = zeros(totalElements,1);
Xfpro = zeros(totalElements,1);
hf = zeros(totalElements,1);
ss = zeros(totalElements,1);
Xf = zeros(totalElements,1);
exposed = zeros(totalElements,1);
rdet = zeros(totalElements,1);
Xhdet = zeros(totalElements,1);
Xidet = zeros(totalElements,1);
Xiads = zeros(totalElements,1);
Xn1det = zeros(totalElements,1);
Xn2det = zeros(totalElements,1);
EPSdet = zeros(totalElements,1);
Xhads = zeros(totalElements,1);
Xn1ads = zeros(totalElements,1);
Xn2ads = zeros(totalElements,1);
EPSads = zeros(totalElements,1);
SFxn1s = zeros(totalElements,1);
SFxn1f = zeros(totalElements,1);
SFxn2s = zeros(totalElements,1);
SFxn2f = zeros(totalElements,1);
SFxhs = zeros(totalElements,1);
SFxhf = zeros(totalElements,1);
KAD5 = zeros(totalElements,1);

% Preallocate final day output variables

CtOClO = zeros(totalElements,24);
NH2ClO = zeros(totalElements,24);
NHCl2O = zeros(totalElements,24);
NHOHClO = zeros(totalElements,24);
ClO = zeros(totalElements,24);
CtNH3O = zeros(totalElements,24);
NO2O = zeros(totalElements,24);
NO3O = zeros(totalElements,24);
N2O = zeros(totalElements,24);
CO2O = zeros(totalElements,24);
O2O = zeros(totalElements,24);
UAPO = zeros(totalElements,24);
BAPO = zeros(totalElements,24);
BOM1O = zeros(totalElements,24);
BOM2O = zeros(totalElements,24);
Xn1sO = zeros(totalElements,24);
Xn1fO = zeros(totalElements,24);
Xn2sO = zeros(totalElements,24);
Xn2fO = zeros(totalElements,24);
EPSsO = zeros(totalElements,24);
EPSfO = zeros(totalElements,24);
XhsO = zeros(totalElements,24);
XhfO = zeros(totalElements,24);
XisO = zeros(totalElements,24);
XifO = zeros(totalElements,24);
XbtO = zeros(totalElements,24);
XatO = zeros(totalElements,24);
NH3O = zeros(totalElements,24);
NH4O = zeros(totalElements,24);
TimeO = zeros(totalElements,24);

%Preallocate individual loss term variables

rUTLbom1 = zeros(totalElements,1);
rand1bom1 = zeros(totalElements,1);
rand2bom1 = zeros(totalElements,1);
rkoxbom1 = zeros(totalElements,1);

rUTLbom2 = zeros(totalElements,1);
rand1bom2 = zeros(totalElements,1);
rand2bom2 = zeros(totalElements,1);
rkdbom2 = zeros(totalElements,1);
rkoxbom2 = zeros(totalElements,1);

rAD3 = zeros(totalElements,1);
rAD1 = zeros(totalElements,1);
rAD2 = zeros(totalElements,1);
rAD5 = zeros(totalElements,1);
rox5 = zeros(totalElements,1);
rUTLNH2Cl = zeros(totalElements,1);
rkox = zeros(totalElements,1);
rkcr = zeros(totalElements,1);
rkSC = zeros(totalElements,1);

%Preallocate mass check variables

totalCO = zeros(days,1);
totalClO = zeros(days,1);
totalNO = zeros(days,1);

%Initialise compartment modelling variable

currentCompartment = 0;

%Create elements

element = 0;

for pipe = 1:totalPipes

    %Get start and end nodes for this pipe
    cStartNode = find(nodes == startNodes(pipe));
    cEndNode = find(nodes == endNodes(pipe));

    for pipeElement=1:numElements(pipe)

        %Create a new element
        element = element + 1;

        %Specify pipe number that new element belongs to
        elementPipes(element,1) = pipe;

        %If first element within pipe
        if pipeElement == 1
            %set start node for element as start node of pipe
            elementStartNodes(element,1) = cStartNode;
        else
            %Set start node for element as the same as the end node
            %of previous element
            elementStartNodes(element,1) = numel(nodes);
        end

        %If last element within pipe
        if pipeElement == numElements(pipe)
            %set end node for element as end node of pipe
            elementEndNodes(element,1) = cEndNode;
        else
            %Add new node to nodes array, and set node value as 1
            %greater than previous largest node value
            nodes(numel(nodes)+1) = max(nodes)+1;
            currentNode = numel(nodes);

            %Set end node for element as newly created node number
            elementEndNodes(element,1) = currentNode;

            %Linear interpolation between known node concentrations
            nCtOCl(currentNode,1) = nCtOCl(cStartNode) + (((nCtOCl(cEndNode) - nCtOCl(cStartNode))/numElements(pipe))*pipeElement);
            nNH2Cl(currentNode,1) = nNH2Cl(cStartNode) + (((nNH2Cl(cEndNode) - nNH2Cl(cStartNode))/numElements(pipe))*pipeElement);
            nNHCl2(currentNode,1) = nNHCl2(cStartNode) + (((nNHCl2(cEndNode) - nNHCl2(cStartNode))/numElements(pipe))*pipeElement);
            nNHOHCl(currentNode,1) = nNHOHCl(cStartNode) + (((nNHOHCl(cEndNode) - nNHOHCl(cStartNode))/numElements(pipe))*pipeElement);
            nCl(currentNode,1) = nCl(cStartNode) + (((nCl(cEndNode) - nCl(cStartNode))/numElements(pipe))*pipeElement);
            nBOM1(currentNode,1) = nBOM1(cStartNode) + (((nBOM1(cEndNode) - nBOM1(cStartNode))/numElements(pipe))*pipeElement);
            nBOM2(currentNode,1) = nBOM2(cStartNode) + (((nBOM2(cEndNode) - nBOM2(cStartNode))/numElements(pipe))*pipeElement);
            nXis(currentNode,1) = nXis(cStartNode) + (((nXis(cEndNode) - nXis(cStartNode))/numElements(pipe))*pipeElement);
            nXhs(currentNode,1) = nXhs(cStartNode) + (((nXhs(cEndNode) - nXhs(cStartNode))/numElements(pipe))*pipeElement);
            nCtNH3(currentNode,1) = nCtNH3(cStartNode) + (((nCtNH3(cEndNode) - nCtNH3(cStartNode))/numElements(pipe))*pipeElement);
            nNO2(currentNode,1) = nNO2(cStartNode) + (((nNO2(cEndNode) - nNO2(cStartNode))/numElements(pipe))*pipeElement);
            nNO3(currentNode,1) = nNO3(cStartNode) + (((nNO3(cEndNode) - nNO3(cStartNode))/numElements(pipe))*pipeElement);
            nN2(currentNode,1) = nN2(cStartNode) + (((nN2(cEndNode) - nN2(cStartNode))/numElements(pipe))*pipeElement);
            nXn1s(currentNode,1) = nXn1s(cStartNode) + (((nXn1s(cEndNode) - nXn1s(cStartNode))/numElements(pipe))*pipeElement);
            nXn2s(currentNode,1) = nXn2s(cStartNode) + (((nXn2s(cEndNode) - nXn2s(cStartNode))/numElements(pipe))*pipeElement);
            nCO2(currentNode,1) = nCO2(cStartNode) + (((nCO2(cEndNode) - nCO2(cStartNode))/numElements(pipe))*pipeElement);
            nUAP(currentNode,1) = nUAP(cStartNode) + (((nUAP(cEndNode) - nUAP(cStartNode))/numElements(pipe))*pipeElement);
            nBAP(currentNode,1) = nBAP(cStartNode) + (((nBAP(cEndNode) - nBAP(cStartNode))/numElements(pipe))*pipeElement);
            nEPSs(currentNode,1) = nEPSs(cStartNode) + (((nEPSs(cEndNode) - nEPSs(cStartNode))/numElements(pipe))*pipeElement);
            nO2(currentNode,1) = nO2(cStartNode) + (((nO2(cEndNode) - nO2(cStartNode))/numElements(pipe))*pipeElement);
            nTime(currentNode,1) = nTime(cStartNode) + (((nTime(cEndNode) - nTime(cStartNode))/numElements(pipe))*pipeElement);
            
            nX(currentNode,1) = nX(cStartNode) + (((nX(cEndNode) - nX(cStartNode))/numElements(pipe))*pipeElement);
            nY(currentNode,1) = nY(cStartNode) + (((nY(cEndNode) - nY(cStartNode))/numElements(pipe))*pipeElement);
        end
        
        %Create reservoir
        
        elementReservoir(element,1) = reservoir(pipe);
        
        if elementReservoir(element) == 1
            if reservoirCompartments == 1
                elementReservoirC1(element,1) = 1;
                elementLengths(element,1) = pElementLengths(pipe);
                elementVolumes(element,1) = pElementVolumes(pipe);
            elseif reservoirCompartments == 2
                if currentCompartment == 1
                    elementReservoirC2(element,1) = 1;
                    elementLengths(element,1) = 0.6*pElementLengths(pipe)*2;
                    elementVolumes(element,1) = 0.6*pElementVolumes(pipe)*2;
                    currentCompartment = 2;
                else
                    elementReservoirC1(element,1) = 1;
                    elementLengths(element,1) = 0.4*pElementLengths(pipe)*2;
                    elementVolumes(element,1) = 0.4*pElementVolumes(pipe)*2;
                    currentCompartment = 1;
                end
            else
                if currentCompartment == 1
                    elementReservoirC2(element,1) = 1;
                    elementLengths(element,1) = 0.6*pElementLengths(pipe)*3;
                    elementVolumes(element,1) = 0.6*pElementVolumes(pipe)*3;
                    currentCompartment = 2;
                elseif currentCompartment == 2
                    elementReservoirC3(element,1) = 1;
                    elementLengths(element,1) = 0.2*pElementLengths(pipe)*3;
                    elementVolumes(element,1) = 0.2*pElementVolumes(pipe)*3;
                    currentCompartment = 3;
                else
                    elementReservoirC1(element,1) = 1;
                    elementLengths(element,1) = 0.2*pElementLengths(pipe)*3;
                    elementVolumes(element,1) = 0.2*pElementVolumes(pipe)*3;
                    currentCompartment = 1;
                end
            end
        else
            elementLengths(element,1) = pElementLengths(pipe);
            elementVolumes(element,1) = pElementVolumes(pipe);
        end

        %Set biofilm concentrations
        
        Xhf(element,1) = pXhf(pipe);
        Xn1f(element,1) = pXn1f(pipe);
        Xn2f(element,1) = pXn2f(pipe);
        EPSf(element,1) = pEPSf(pipe);
        Xif(element,1) = pXif(pipe);

        %Set material variables
        
        if strcmp(materials(pipe),'D.I.')
            kcr3(element,1) = 1.7E-2; 
            kSC1(element,1) = 0;
        elseif strcmp(materials(pipe),'R.C.')
            kcr3(element,1) = 0;
            kSC1(element,1) = 1.30E6;
        elseif strcmp(materials(pipe),'C.L.')
            kcr3(element,1) = 0;
            kSC1(element,1) = 1.30E5;
        end
    end
end

elementV = 1:numel(elements);

nodeEndElements = {zeros(totalNodes,24)};
nodeStartElements = {zeros(totalNodes,24)};

totalNodes = numel(nodes);

%Import Booster Chlorimination sheet from input file

bNH2Cl = zeros(totalNodes,1);
booster = xlsread(filename,'Booster Chlorimination');

% If there is at least one booster node
if size(booster) ~= [0,0]
    bNodeNums = booster(:,1);
    bNH2Cli = booster(:,2);
    for node = 1:numel(bNodeNums)
        %Find node number in nodes array
        nodeNum = find(nodes == bNodeNums(node));
        bNH2Cl(nodeNum,1) = bNH2Cli(node);
    end
end

TFb = Tbio^(Temp-20.0);
TFb2 = Tbio2^(Temp-20.0);
TFc = Tchem^(Temp-20.0);
H = 10^logH;
OH = (1.0D-14)/H;
H2CO3 = CtCO3*((H*H)/(H*H + H*KeHCO3 + KeHCO3*KeCO3));
HCO3 = CtCO3*((H*KeHCO3)/(H*H + H*KeHCO3 + KeHCO3*KeCO3));
CO3 = CtCO3*((KeHCO3*KeCO3)/(H*H + H*KeHCO3 + KeHCO3*KeCO3));
KAD1 = kH.*H+kH2CO3.*H2CO3;

%Initialise mass check variables

addedC = 0;
addedCl = 0;
addedN = 0;
lostC = 0;
lostCl = 0;
lostN = 0;
externalC = 0;
externalCl = 0;
externalN = 0;

for day = 1:days

    for hour = 1:24

        % Display current day and hour
        disp(day);
        disp(hour);

        for timeStep = 1:timeSteps

            % At the start of each hour on the first day
            if day == 1 && timeStep == 1
                for pipe = 1:totalPipes
                    Qin(pipe,hour) = flows(pipe,(hour*2)-1);
                    Qout(pipe,hour) = flows(pipe,(hour*2));
                end

                % Compute flow directions for every hour
                for element = 1:numel(elements)
                    if Qin(elementPipes(element),hour) < 0 || Qout(elementPipes(element),hour) < 0
                        % Switch start and end nodes
                        cElementStartNodes(element,hour) = elementEndNodes(element);
                        cElementEndNodes(element,hour) = elementStartNodes(element);
                    else
                        % Keep start and end nodes the same
                        cElementStartNodes(element,hour) = elementStartNodes(element);
                        cElementEndNodes(element,hour) = elementEndNodes(element);
                    end
                end

                % Find elements connected to every node for every hour
                for node = 1:numel(nodes)

                    %Find elements which have this node as end node
                    %Returns vector specifying elements with this node as
                    %end node
                    nodeEndElements{node,hour} = find(cElementEndNodes(:,hour) == node);
                    nodeStartElements{node,hour} = find(cElementStartNodes(:,hour) == node);

                end

                % At the start of the program, set element concentrations to
                % end node concentrations
                if hour == 1

                    CtOCl = (nCtOCl(cElementStartNodes(elementV,hour))+nCtOCl(cElementEndNodes(elementV,hour)))/2;
                    NH2Cl = (nNH2Cl(cElementStartNodes(elementV,hour))+nNH2Cl(cElementEndNodes(elementV,hour)))/2;
                    NHCl2 = (nNHCl2(cElementStartNodes(elementV,hour))+nNHCl2(cElementEndNodes(elementV,hour)))/2;
                    NHOHCl = (nNHOHCl(cElementStartNodes(elementV,hour))+nNHOHCl(cElementEndNodes(elementV,hour)))/2;
                    Cl = (nCl(cElementStartNodes(elementV,hour))+nCl(cElementEndNodes(elementV,hour)))/2;
                    BOM1 = (nBOM1(cElementStartNodes(elementV,hour))+nBOM1(cElementEndNodes(elementV,hour)))/2;
                    BOM2 = (nBOM2(cElementStartNodes(elementV,hour))+nBOM2(cElementEndNodes(elementV,hour)))/2;
                    Xis = (nXis(cElementStartNodes(elementV,hour))+nXis(cElementEndNodes(elementV,hour)))/2;
                    Xhs = (nXhs(cElementStartNodes(elementV,hour))+nXhs(cElementEndNodes(elementV,hour)))/2;
                    CtNH3 = (nCtNH3(cElementStartNodes(elementV,hour))+nCtNH3(cElementEndNodes(elementV,hour)))/2;
                    NO2 = (nNO2(cElementStartNodes(elementV,hour))+nNO2(cElementEndNodes(elementV,hour)))/2;
                    NO3 = (nNO3(cElementStartNodes(elementV,hour))+nNO3(cElementEndNodes(elementV,hour)))/2;
                    N2 = (nN2(cElementStartNodes(elementV,hour))+nN2(cElementEndNodes(elementV,hour)))/2;
                    Xn1s = (nXn1s(cElementStartNodes(elementV,hour))+nXn1s(cElementEndNodes(elementV,hour)))/2;
                    Xn2s = (nXn2s(cElementStartNodes(elementV,hour))+nXn2s(cElementEndNodes(elementV,hour)))/2;
                    CO2 = (nCO2(cElementStartNodes(elementV,hour))+nCO2(cElementEndNodes(elementV,hour)))/2;
                    UAP = (nUAP(cElementStartNodes(elementV,hour))+nUAP(cElementEndNodes(elementV,hour)))/2;
                    BAP = (nBAP(cElementStartNodes(elementV,hour))+nBAP(cElementEndNodes(elementV,hour)))/2;
                    O2 = (nO2(cElementStartNodes(elementV,hour))+nO2(cElementEndNodes(elementV,hour)))/2;
                    EPSs = (nEPSs(cElementStartNodes(elementV,hour))+nEPSs(cElementEndNodes(elementV,hour)))/2;
                    Time = (nTime(cElementStartNodes(elementV,hour))+nTime(cElementEndNodes(elementV,hour)))/2;

                    EPS = EPSs + EPSf;
                    Xh = Xhs + Xhf;
                    Xn1 = Xn1s + Xn1f;
                    Xn2 = Xn2s + Xn2f;
                    Xi = Xis + Xif;

                    HOCl = (CtOCl.*H)/(H+KeOCl);
                    HOCl(HOCl<0) = 0;

                    OCl = (CtOCl.*KeOCl)/(H+KeOCl);
                    OCl(OCl<0) = 0;

                    totalC(1) = sum(MC.*(BOM1 + BOM2 + Xhs + Xhf + Xn1s + Xn1f + Xn2s + Xn2f + EPSs + EPSf + Xis + Xif + UAP + BAP + CO2).*elementVolumes.*1000);
                    totalCl(1) = sum((HOCl + OCl + NH2Cl + 2.*NHCl2 + NHOHCl + Cl).*elementVolumes.*1000); 
                    totalN(1) = sum((NH2Cl + NHCl2 + NHOHCl + CtNH3 + NO2 + NO3 + 2.*N2 + MN.*(BOM1 + BOM2 + Xhs + Xhf + Xn1s + Xn1f + Xn2s + Xn2f + EPSs + EPSf + Xis + Xif + UAP + BAP)).*elementVolumes.*1000);

                end
            end
            
            %Define latest node concentrations
            %Loop through all nodes
            for node = 1:numel(nodes)

                denominator = 0;
                nCtOCl(node,1) = 0;
                nNH2Cl(node,1) = 0;
                nNHCl2(node,1) = 0;
                nNHOHCl(node,1) = 0;
                nCl(node,1) = 0;
                nBOM1(node,1) = 0;
                nBOM2(node,1) = 0;
                nXis(node,1) = 0;
                nXhs(node,1) = 0;
                nCtNH3(node,1) = 0;
                nNO2(node,1) = 0;
                nNO3(node,1) = 0;
                nN2(node,1) = 0;
                nXn1s(node,1) = 0;
                nXn2s(node,1) = 0;
                nCO2(node,1) = 0;
                nUAP(node,1) = 0;
                nBAP(node,1) = 0;
                nEPSs(node,1) = 0;
                nO2(node,1) = 0;
                nTime(node,1) = 0;

                cQin = 0;

                if isempty(nodeEndElements{node,hour}) == 0

                    %For all elements linked to his node
                    cQin = abs(Qout(elementPipes(nodeEndElements{node,hour}),hour));
                    denominator = sum(cQin);

                    if (denominator > 0)

                        nCtOCl(node,1) = sum(CtOCl(nodeEndElements{node,hour}).*cQin)/denominator;
                        nNH2Cl(node,1) = sum(NH2Cl(nodeEndElements{node,hour}).*cQin)/denominator;
                        nNHCl2(node,1) = sum(NHCl2(nodeEndElements{node,hour}).*cQin)/denominator;
                        nNHOHCl(node,1) = sum(NHOHCl(nodeEndElements{node,hour}).*cQin)/denominator;
                        nCl(node,1) = sum(Cl(nodeEndElements{node,hour}).*cQin)/denominator;
                        nBOM1(node,1) = sum(BOM1(nodeEndElements{node,hour}).*cQin)/denominator;
                        nBOM2(node,1) = sum(BOM2(nodeEndElements{node,hour}).*cQin)/denominator;
                        nXis(node,1) = sum(Xis(nodeEndElements{node,hour}).*cQin)/denominator;
                        nXhs(node,1) = sum(Xhs(nodeEndElements{node,hour}).*cQin)/denominator;
                        nCtNH3(node,1) = sum(CtNH3(nodeEndElements{node,hour}).*cQin)/denominator;
                        nNO2(node,1) = sum(NO2(nodeEndElements{node,hour}).*cQin)/denominator;
                        nNO3(node,1) = sum(NO3(nodeEndElements{node,hour}).*cQin)/denominator;
                        nN2(node,1) = sum(N2(nodeEndElements{node,hour}).*cQin)/denominator;
                        nXn1s(node,1) = sum(Xn1s(nodeEndElements{node,hour}).*cQin)/denominator;
                        nXn2s(node,1) = sum(Xn2s(nodeEndElements{node,hour}).*cQin)/denominator;
                        nCO2(node,1) = sum(CO2(nodeEndElements{node,hour}).*cQin)/denominator;
                        nUAP(node,1) = sum(UAP(nodeEndElements{node,hour}).*cQin)/denominator;
                        nBAP(node,1) = sum(BAP(nodeEndElements{node,hour}).*cQin)/denominator;
                        nEPSs(node,1) = sum(EPSs(nodeEndElements{node,hour}).*cQin)/denominator;
                        nO2(node,1) = sum(O2(nodeEndElements{node,hour}).*cQin)/denominator;
                        nTime(node,1) = sum(Time(nodeEndElements{node,hour}).*cQin)/denominator;
                        
                    end
                end
            end

            % Set treatment plant node concentrations
            nCtOCl(TPnode,1) = aCtOCl(TPnode);
            nNH2Cl(TPnode,1) = aNH2Cl(TPnode);
            nNHCl2(TPnode,1) = aNHCl2(TPnode);
            nNHOHCl(TPnode,1) = aNHOHCl(TPnode);
            nCl(TPnode,1) = aCl(TPnode);
            nBOM1(TPnode,1) = aBOM1(TPnode);
            nBOM2(TPnode,1) = aBOM2(TPnode);
            nXis(TPnode,1) = aXis(TPnode);
            nXhs(TPnode,1) = aXhs(TPnode);
            nCtNH3(TPnode,1) = aCtNH3(TPnode);
            nNO2(TPnode,1) = aNO2(TPnode);
            nNO3(TPnode,1) = aNO3(TPnode);
            nN2(TPnode,1) = aN2(TPnode);
            nXn1s(TPnode,1) = aXn1s(TPnode);
            nXn2s(TPnode,1) = aXn2s(TPnode);
            nCO2(TPnode,1) = aCO2(TPnode);
            nUAP(TPnode,1) = aUAP(TPnode);
            nBAP(TPnode,1) = aBAP(TPnode);
            nO2(TPnode,1) = aO2(TPnode);
            nEPSs(TPnode,1) = aEPSs(TPnode);
            nTime(TPnode,1) = aTime(TPnode);

            % Find treatment plant elements for every hour
            if timeStep == 1
                treatmentPlantElements = find(cElementStartNodes(:,hour) == TPnode);
            end
            
            % Calculate total mass added to the system from the treatment plant
            if isempty(treatmentPlantElements) == 0
                tpQin = sum(abs(Qin(elementPipes(treatmentPlantElements),hour)));
                
                addedC = addedC + timeInc*(MC*(aBOM1(TPnode) + aBOM2(TPnode) + aXhs(TPnode) + aXn1s(TPnode) + aXn2s(TPnode) + aEPSs(TPnode) + aXis(TPnode) + aUAP(TPnode) + aBAP(TPnode) + aCO2(TPnode)))*tpQin*1000;
                addedCl = addedCl + timeInc*(aNH2Cl(TPnode) + 2*aNHCl2(TPnode) + aNHOHCl(TPnode) + aCl(TPnode))*tpQin*1000;
                addedN = addedN + timeInc*(aNH2Cl(TPnode) + aNHCl2(TPnode) + aNHOHCl(TPnode) + aCtNH3(TPnode) + aNO2(TPnode) + aNO3(TPnode) + 2*aN2(TPnode) + MN*(aBOM1(TPnode) + aBOM2(TPnode) + aXhs(TPnode) + aXn1s(TPnode) + aXn2s(TPnode) + aXis(TPnode) + aEPSs(TPnode) + aUAP(TPnode) + aBAP(TPnode)))*tpQin*1000;
            end

            for node = 1:numel(nodes)

                demandQin = sum(abs(Qout(elementPipes(nodeEndElements{node,hour}),hour)));
                demandQout = sum(abs(Qin(elementPipes(nodeStartElements{node,hour}),hour)));

                demandQ = demandQin - demandQout;
                
                % Calculate total mass lost from the system
                
                if demandQ > 0
                    lostC = lostC + timeInc*(MC*(nBOM1(node) + nBOM2(node) + nXhs(node) + nXn1s(node) + nXn2s(node) + nXis(node) + nEPSs(node) + nUAP(node) + nBAP(node) + nCO2(node)))*demandQ*1000;
                    lostCl = lostCl + timeInc*(nNH2Cl(node) + 2*nNHCl2(node) + nNHOHCl(node) + nCl(node))*demandQ*1000;
                    lostN = lostN + timeInc*(nNH2Cl(node) + nNHCl2(node) + nNHOHCl(node) + nCtNH3(node) + nNO2(node) + nNO3(node) + 2*nN2(node) + MN*(nBOM1(node) + nBOM2(node) + nXhs(node) + nXn1s(node) + nXn2s(node) + nXis(node) + nEPSs(node) + nUAP(node) + nBAP(node)))*demandQ*1000;
                end
                
                % Calculate total mass added to the system

                if demandQ < 0 && node ~= TPnode
                    externalC = externalC - timeInc*(MC*(nBOM1(node) + nBOM2(node) + nXhs(node) + nXn1s(node) + nXn2s(node) + nXis(node) + nEPSs(node) + nUAP(node) + nBAP(node) + nCO2(node)))*demandQ*1000;
                    externalCl = externalCl - timeInc*(nNH2Cl(node) + 2*nNHCl2(node) + nNHOHCl(node) + nCl(node))*demandQ.*1000;
                    externalN = externalN - timeInc*(nNH2Cl(node) + nNHCl2(node) + nNHOHCl(node) + nCtNH3(node) + nNO2(node) + nNO3(node) + 2*nN2(node) + MN*(nBOM1(node) + nBOM2(node) + nXhs(node) + nXn1s(node) + nXn2s(node) + nXis(node) + nEPSs(node) + nUAP(node) + nBAP(node)))*demandQ*1000;
                end
                
                % Calculate total mass added to the system from booster
                % nodes
                
                if bNH2Cl(node) > 0
                    nNH2Cl(node,1) = nNH2Cl(node) + bNH2Cl(node);
                    bQout = sum(abs(Qout(elementPipes(nodeStartElements{node,hour}),hour)));
                    addedCl = addedCl + timeInc*(bNH2Cl(node))*bQout*1000;
                    addedN = addedN + timeInc*(bNH2Cl(node))*bQout*1000;
                end
                
            end

            NH4 = (CtNH3.*H)./(H+KeNH3);
            NH4(NH4<1E-14) = 1E-14;

            NH3 = (CtNH3.*KeNH3)./(H+KeNH3);
            NH3(NH3<1E-14) = 1E-14;

            HOCl = (CtOCl.*H)/(H+KeOCl);
            HOCl(HOCl<0) = 0;

            OCl = (CtOCl.*KeOCl)/(H+KeOCl);
            OCl(OCl<0) = 0;

            KAD5(NH2Cl<=1E-12,1) = 0;
            KAD5(NH2Cl>1E-12,1) = 5.0E-3;


            %Adjust flexible kinetic constants

            kox5 = (kox5a.*H.*(1+kox5b.*NO2))./(kox5c.*NH3 + (1+kox5b.*NO2));
            kox4 = (kox4a + kox4b.*NO2)./OH;


            %Substrate utilisation rate terms

            UTLbom1 = qh.*(BOM1./(Kbom1 + BOM1)).*(O2./(KO2+O2));

            UTLbom2 = qh.*(BOM2./(Kbom2 + BOM2)).*(O2./(KO2+O2));

            UTLbom = UTLbom1 + UTLbom2;

            UTLuap = quap.*(UAP./(Kuap + UAP)).*(O2./(KO2+O2));

            UTLbap = qbap.*(BAP./(Kbap + BAP)).*(O2./(KO2+O2));

            UTLsmp = UTLuap + UTLbap;

            UTLNH3 = qn1.*(NH3./(KNH3.*(1+(NH2Cl./KNH2Cl)) + NH3)).*(O2./(KO2+O2));
            UTLNH3(CO2<=0) = 0;

            UTLNO2 = qn2.*(NO2./(KNO2+ NO2)).*(O2./(KO2+O2));
            UTLNO2(CO2<=0) = 0;

            UTLNH2Cl = qNH2Cl.*(NH2Cl./(KNH2Cl.*(1+(NH3./KNH3)) + NH2Cl)).*(O2./(KO2+O2));

            UTLNHOHCl = qNHOHCl.*(NHOHCl./(KNHOHCl + NHOHCl)).*(O2./(KO2+O2));

            and1bom1 = qbom1.*and1.*(BOM1./(Kbom1 + BOM1)).*(KIDO./(KIDO + O2)).*(NO2./(KNO2 + NO2));

            and1bom2 = qbom2.*and1.*(BOM2./(Kbom2 + BOM2)).*(KIDO./(KIDO + O2)).*(NO2./(KNO2 + NO2));

            and1bom = and1bom1 + and1bom2;

            and2bom1 = qbom1.*and2.*(BOM1./(Kbom1 + BOM1)).*(KIDO./(KIDO + O2)).*(NO3./(KNO3 + NO3));

            and2bom2 = qbom2.*and2.*(BOM2./(Kbom2 + BOM2)).*(KIDO./(KIDO + O2)).*(NO3./(KNO3 + NO3));

            and2bom = and2bom1 + and2bom2;

            andbom = and1bom + and2bom;

            and1uap = quap.*and1.*(UAP./(Kuap + UAP)).*(KIDO./(KIDO + O2)).*(NO2./(KNO2 + NO2));

            and1bap = qbap.*and1.*(BAP./(Kbap + BAP)).*(KIDO./(KIDO + O2)).*(NO2./(KNO2 + NO2));

            and1smp = and1uap + and1bap;

            and2uap = quap.*and2.*(UAP./(Kuap + UAP)).*(KIDO./(KIDO + O2)).*(NO3./(KNO3 + NO3));

            and2bap = qbap.*and2.*(BAP./(Kbap + BAP)).*(KIDO./(KIDO + O2)).*(NO3./(KNO3 + NO3));

            and2smp = and2uap + and2bap;

            andsmp = and1smp + and2smp;


            %Biofilm terms

            Lf = kpro1.*(1./C(elementPipes)).^kpro2;

            Xfpro = Lf.*BFden.*(4./diameters(elementPipes).*1000);

            hf = (6.78.*((abs(Qin(elementPipes,hour))./(areas(elementPipes)*60*60)).^1.85).*elementLengths) ...
                ./((C(elementPipes).^1.85).*(diameters(elementPipes).^1.165));

            ss = (hf.*dens.*g.*diameters(elementPipes))./(4.*elementLengths);

            Xf = Xhf + Xn1f + Xn2f + Xif + EPSf;

            exposed(Xf>0 & Xf>Xfpro,1) = (Xf(Xf>0 & Xf>Xfpro) - Xfpro(Xf>0 & Xf>Xfpro))./Xf(Xf>0 & Xf>Xfpro);
            exposed(Xf<=0 | Xf<=Xfpro,1) = 0;

            rdet = kdet.*(ss.^AS);

            Xhdet = rdet.*Xhf.*exposed;
            Xhdet(Xhdet<=0) = 0;
            Xhdet(Xhdet>Xhf) = Xhf(Xhdet>Xhf);

            Xn1det = rdet.*Xn1f.*exposed;
            Xn1det(Xn1det<=0) = 0;
            Xn1det(Xn1det>Xn1f) = Xn1f(Xn1det>Xn1f);

            Xn2det = rdet.*Xn2f.*exposed;
            Xn2det(Xn2det<=0) = 0;
            Xn2det(Xn2det>Xn2f) = Xn2f(Xn2det>Xn2f);

            Xidet = rdet.*Xif.*exposed;
            Xidet(Xidet<=0) = 0;
            Xidet(Xidet>Xif) = Xif(Xidet>Xif);

            EPSdet = rdet.*EPSf.*exposed;
            EPSdet(EPSdet<=0) = 0;
            EPSdet(EPSdet>EPSf) = EPSf(EPSdet>EPSf);

            Xhads = kads.*Xhs.*(4./diameters(elementPipes).*1000);
            Xhads(Xhads>Xhs) = Xhs(Xhads>Xhs);

            Xn1ads = kads.*Xn1s.*(4./diameters(elementPipes).*1000);
            Xn1ads(Xn1ads>Xn1s) = Xn1s(Xn1ads>Xn1s);

            Xn2ads = kads.*Xn2s.*(4./diameters(elementPipes).*1000);
            Xn2ads(Xn2ads>Xn2s) = Xn2s(Xn2ads>Xn2s);

            Xiads = kads.*Xis.*(4./diameters(elementPipes).*1000);
            Xiads(Xiads>Xis) = Xis(Xiads>Xis);

            EPSads = kads.*EPSs.*(4./diameters(elementPipes).*1000);
            EPSads(EPSads>EPSs) = EPSs(EPSads>EPSs);


            %Mass-balance equations
            
            if reservoirCompartments == 1
                elementVolumesNew(elementReservoir==1 & elementReservoirC1==1,1) = elementVolumes(elementReservoir==1 & elementReservoirC1==1) + (timeInc.*abs(Qin(elementPipes(elementReservoir==1 & elementReservoirC1==1),hour))) - (timeInc.*abs(Qout(elementPipes(elementReservoir==1 & elementReservoirC1==1),hour)));
            elseif reservoirCompartments == 2
                elementVolumesNew(elementReservoir==1 & elementReservoirC1==1,1) = elementVolumes(elementReservoir==1 & elementReservoirC1==1);
                elementVolumesNew(elementReservoir==1 & elementReservoirC2==1,1) = elementVolumes(elementReservoir==1 & elementReservoirC2==1) + (timeInc.*abs(Qin(elementPipes(elementReservoir==1 & elementReservoirC2==1),hour))) - (timeInc.*abs(Qout(elementPipes(elementReservoir==1 & elementReservoirC2==1),hour)));
            elseif reservoirCompartments == 3
                elementVolumesNew(elementReservoir==1 & elementReservoirC1==1,1) = elementVolumes(elementReservoir==1 & elementReservoirC1==1);
                elementVolumesNew(elementReservoir==1 & elementReservoirC2==1,1) = elementVolumes(elementReservoir==1 & elementReservoirC2==1) + (timeInc.*abs(Qin(elementPipes(elementReservoir==1 & elementReservoirC2==1),hour))) - (timeInc.*abs(Qout(elementPipes(elementReservoir==1 & elementReservoirC2==1),hour)));
                elementVolumesNew(elementReservoir==1 & elementReservoirC3==1,1) = elementVolumes(elementReservoir==1 & elementReservoirC3==1);
            end
            elementVolumesRatio(elementReservoir==1) = elementVolumes(elementReservoir==1)./elementVolumesNew(elementReservoir==1);
            elementVolumesRatio(elementReservoir~=1) = 1;
            
            % Update reservoir volumes
            elementVolumes(elementReservoir==1) = elementVolumesNew(elementReservoir==1);
            if any(elementVolumes(elementReservoir==1) <= 0)
                error('Reservoir empty.')
            end
            elementLengths(elementReservoir==1) = elementVolumes(elementReservoir==1)./areas(elementPipes(elementReservoir==1));

            TimeNew = Time ...
                + timeInc.*((abs(Qin(elementPipes,hour)./elementVolumes)).*nTime(cElementStartNodes(elementV,hour)) ...
                - (abs(Qout(elementPipes,hour)./elementVolumes)).*Time + 1.00);

            XhsNew = advection(Xhs,nXhs);
            Xn1sNew = advection(Xn1s,nXn1s);
            Xn2sNew = advection(Xn2s,nXn2s);
            XisNew = advection(Xis,nXis);
            EPSsNew = advection(EPSs,nEPSs);
            BOM1New = advection(BOM1,nBOM1);
            BOM2New = advection(BOM2,nBOM2);
            UAPNew = advection(UAP,nUAP);
            BAPNew = advection(BAP,nBAP);
            CtNH3New = advection(CtNH3,nCtNH3);
            NO2New = advection(NO2,nNO2);
            NO3New = advection(NO3,nNO3);
            N2New = advection(N2,nN2);
            CtOClNew = advection(CtOCl,nCtOCl);
            NH2ClNew = advection(NH2Cl,nNH2Cl);
            NHCl2New = advection(NHCl2,nNHCl2);
            NHOHClNew = advection(NHOHCl,nNHOHCl);
            ClNew = advection(Cl,nCl);
            CO2New = advection(CO2,nCO2);
            O2New = advection(O2,nO2);
            
            
            UTLNH3loss(:,1) = timeInc.*(TFb2.*(Yn1.*(1-kuapn1-kepsn1)+((kuapn1./MN)+(kepsn1./MN))).*UTLNH3.*Xn1);
            UTLNO2loss(:,1) = timeInc.*(TFb2.*(Yn2.*(1-kuapn2-kepsn2)+((kuapn2./MN)+(kepsn2./MN))).*UTLNO2.*Xn2);
            
            totalUTLloss(:,1) = UTLNH3loss + UTLNO2loss;
            
            CO2NewT = CO2New ...
                + timeInc.*(TFb.*(1-Yh.*(1-kuaph-kepsh) - (kuaph+kepsh)).*UTLbom.*Xh ...
                + TFb.*(1-Yh.*(1-kuaph-kepsh) - (kuaph+kepsh)).*and1bom.*Xh ...
                + TFb.*(1-Yh.*(1-kuaph-kepsh) - (kuaph+kepsh)).*and2bom.*Xh ...
                + TFb.*(1-Yp).*and1smp.*Xh ...
                + TFb.*(1-Yp).*and2smp.*Xh ...
                + TFb.*(1-Yp).*UTLsmp.*Xh ...
                - TFb2.*(Yn1.*(1-kuapn1-kepsn1)+((kuapn1./MN)+(kepsn1./MN))).*UTLNH3.*Xn1 ...
                - TFb2.*(Yn2.*(1-kuapn2-kepsn2)+((kuapn2./MN)+(kepsn2./MN))).*UTLNO2.*Xn2 ...
                + fd.*(bh.*((O2./(KO2+O2))+(NO2./(KNO2+NO2)).*(KIDO./(KIDO.*O2))+(NO3./(KNO3+NO3)).*(KIDO./(KIDO.*O2))).*Xh + bn1.*(O2./(KO2+O2)).*Xn1 + bn2.*(O2./(KO2+O2)).*Xn2) ...
                + TFc.*(kox3x1.*BOM1 + kox3x2.*BOM2 + kox3x2.*BAP + kox3x2.*EPS + kox3x2.*UAP + kox3x2.*Xi).*NH2Cl ...
                + TFc.*(kox1x1.*BOM1 + kox1x2.*BOM2 + kox1x2.*BAP + kox1x2.*EPS + kox1x2.*UAP + kox1x2.*Xi).*HOCl ...
                + TFc.*(kox2x1.*BOM1 + kox2x2.*BOM2 + kox2x2.*BAP + kox2x2.*EPS + kox2x2.*UAP + kox2x2.*Xi).*OCl);
            
            addTotal(CO2NewT<0,1) = abs(CO2NewT(CO2NewT<0));
            addTotal(CO2NewT>=0,1) = 0;
            scaleRatio(totalUTLloss>0,1) = 1-(addTotal(totalUTLloss>0)./totalUTLloss(totalUTLloss>0));
            scaleRatio(totalUTLloss<=0,1) = 1;
            UTLNH3 = UTLNH3.*scaleRatio;
            UTLNO2 = UTLNO2.*scaleRatio;

            CO2New = CO2New ...
                + timeInc.*(TFb.*(1-Yh.*(1-kuaph-kepsh) - (kuaph+kepsh)).*UTLbom.*Xh ...
                + TFb.*(1-Yh.*(1-kuaph-kepsh) - (kuaph+kepsh)).*and1bom.*Xh ...
                + TFb.*(1-Yh.*(1-kuaph-kepsh) - (kuaph+kepsh)).*and2bom.*Xh ...
                + TFb.*(1-Yp).*and1smp.*Xh ...
                + TFb.*(1-Yp).*and2smp.*Xh ...
                + TFb.*(1-Yp).*UTLsmp.*Xh ...
                - TFb2.*(Yn1.*(1-kuapn1-kepsn1)+((kuapn1./MN)+(kepsn1./MN))).*UTLNH3.*Xn1 ...
                - TFb2.*(Yn2.*(1-kuapn2-kepsn2)+((kuapn2./MN)+(kepsn2./MN))).*UTLNO2.*Xn2 ...
                + fd.*(bh.*((O2./(KO2+O2))+(NO2./(KNO2+NO2)).*(KIDO./(KIDO.*O2))+(NO3./(KNO3+NO3)).*(KIDO./(KIDO.*O2))).*Xh + bn1.*(O2./(KO2+O2)).*Xn1 + bn2.*(O2./(KO2+O2)).*Xn2) ...
                + TFc.*(kox3x1.*BOM1 + kox3x2.*BOM2 + kox3x2.*BAP + kox3x2.*EPS + kox3x2.*UAP + kox3x2.*Xi).*NH2Cl ...
                + TFc.*(kox1x1.*BOM1 + kox1x2.*BOM2 + kox1x2.*BAP + kox1x2.*EPS + kox1x2.*UAP + kox1x2.*Xi).*HOCl ...
                + TFc.*(kox2x1.*BOM1 + kox2x2.*BOM2 + kox2x2.*BAP + kox2x2.*EPS + kox2x2.*UAP + kox2x2.*Xi).*OCl);
                    
            CO2New(CO2New<0) = 0;
            
            
            XhsNew = XhsNew ...
                + timeInc.*(TFb.*Yh.*(1-kuaph-kepsh).*UTLbom.*Xhs ...
                + TFb.*Yh.*(1-kuaph-kepsh).*and1bom.*Xhs ...
                + TFb.*Yh.*(1-kuaph-kepsh).*and2bom.*Xhs ...
                + TFb.*Yp.*UTLsmp.*Xhs ...
                + TFb.*Yp.*and1smp.*Xhs ...
                + TFb.*Yp.*and2smp.*Xhs ...
                - bh.*((O2./(KO2+O2))+(NO2./(KNO2+NO2)).*(KIDO./(KIDO.*O2))+(NO3./(KNO3+NO3)).*(KIDO./(KIDO.*O2))).*Xhs ...
                - (kd1x1.*HOCl + kd2x1.*OCl+ kd3x1.*NH2Cl).*Xhs ...
                + Xhdet ...
                - Xhads);
            XhsNew(XhsNew<0) = 0;

            XhfNew = (Xhf ...
                + timeInc.*(TFb.*Yh.*(1-kuaph-kepsh).*UTLbom.*Xhf ...
                + TFb.*Yh.*(1-kuaph-kepsh).*and1bom.*Xhf ...
                + TFb.*Yh.*(1-kuaph-kepsh).*and2bom.*Xhf ...
                + TFb.*Yp.*UTLsmp.*Xhf ...
                + TFb.*Yp.*and1smp.*Xhf ...
                + TFb.*Yp.*and2smp.*Xhf ...
                - bh.*((O2./(KO2+O2))+(NO2./(KNO2+NO2)).*(KIDO./(KIDO.*O2))+(NO3./(KNO3+NO3)).*(KIDO./(KIDO.*O2))).*Xhf ...
                - (kd1x2.*HOCl + kd2x2.*OCl + kd3x2.*NH2Cl).*Xhf ...
                - Xhdet ...
                + Xhads)).*elementVolumesRatio;
            XhfNew(XhfNew<0) = 0;

            XhNew = XhsNew + XhfNew;
            
            Xn1sNew = Xn1sNew ...
                + timeInc.*(TFb2.*Yn1.*(1 - kuapn1 - kepsn1).*UTLNH3.*Xn1s ...
                - bn1.*(O2./(KO2+O2)).*Xn1s ...
                - (kd1x3.*HOCl + kd2x3.*OCl + kd3x3.*NH2Cl).*Xn1s ...
                + Xn1det ...
                - Xn1ads);
            Xn1sNew(Xn1sNew<0) = 0;
            
            Xn1fNew = (Xn1f ...
                + timeInc.*(TFb2.*Yn1.*(1 - kuapn1 - kepsn1).*UTLNH3.*Xn1f ...
                - bn1.*(O2./(KO2+O2)).*Xn1f ...
                - (kd1x4.*HOCl + kd2x4.*OCl + kd3x4.*NH2Cl).*Xn1f ...
                - Xn1det ...
                + Xn1ads)).*elementVolumesRatio;
            Xn1fNew(Xn1fNew<0) = 0;

            Xn1New = Xn1sNew + Xn1fNew;
            
            Xn2sNew = Xn2sNew ...
                + timeInc.*(TFb2.*Yn2.*(1 - kuapn2 - kepsn2).*UTLNO2.*Xn2s ...
                - bn2.*(O2./(KO2+O2)).*Xn2s ...
                - (kd1x5.*HOCl + kd2x5.*OCl + kd3x5.*NH2Cl).*Xn2s ...
                + Xn2det ...
                - Xn2ads);
            Xn2sNew(Xn2sNew<0) = 0;

            Xn2fNew = (Xn2f ...
                + timeInc.*(TFb2.*Yn2.*(1 - kuapn2 - kepsn2).*UTLNO2.*Xn2f...
                - bn2.*(O2./(KO2+O2)).*Xn2f ...
                - (kd1x6.*HOCl + kd2x6.*OCl + kd3x6.*NH2Cl).*Xn2f ...
                - Xn2det ...
                + Xn2ads)).*elementVolumesRatio;
            Xn2fNew(Xn2fNew<0) = 0;

            Xn2New = Xn2sNew + Xn2fNew;

            XisNew = XisNew ...
                + timeInc.*((1-fd).*bh.*((O2./(KO2+O2))+(NO2./(KNO2+NO2)).*(KIDO./(KIDO.*O2))+(NO3./(KNO3+NO3)).*(KIDO./(KIDO.*O2))).*Xhs ...
                + (1-fd).*bn1.*(O2./(KO2+O2)).*Xn1s ...
                + (1-fd).*bn2.*(O2./(KO2+O2)).*Xn2s ...
                + (1-fd).*(kd1x1.*HOCl + kd2x1.*OCl + kd3x1.*NH2Cl).*Xhs ...
                + (1-fd).*(kd1x3.*HOCl + kd2x3.*OCl + kd3x3.*NH2Cl).*Xn1s ...
                + (1-fd).*(kd1x5.*HOCl + kd2x5.*OCl + kd3x5.*NH2Cl).*Xn2s...
                - TFc.*(kox1x2.*HOCl + kox2x2.*OCl + kox3x2.*NH2Cl).*Xis ...
                + Xidet ...
                - Xiads);
            XisNew(XisNew<0) = 0;

            XifNew = (Xif ...
                + timeInc.*((1-fd).*bh.*((O2./(KO2+O2))+(NO2./(KNO2+NO2)).*(KIDO./(KIDO.*O2))+(NO3./(KNO3+NO3)).*(KIDO./(KIDO.*O2))).*Xhf ...
                + (1-fd).*bn1.*(O2./(KO2+O2)).*Xn1f ...
                + (1-fd).*bn2.*(O2./(KO2+O2)).*Xn2f ...
                + (1-fd).*(kd1x2.*HOCl + kd2x2.*OCl + kd3x2.*NH2Cl).*Xhf ...
                + (1-fd).*(kd1x4.*HOCl + kd2x4.*OCl + kd3x4.*NH2Cl).*Xn1f ...
                + (1-fd).*(kd1x6.*HOCl + kd2x6.*OCl + kd3x6.*NH2Cl).*Xn2f ...
                - TFc.*(kox1x2.*HOCl + kox2x2.*OCl + kox3x2.*NH2Cl).*Xif ...
                - Xidet ...
                + Xiads)).*elementVolumesRatio;
            XifNew(XifNew<0) = 0;

            XiNew = XisNew + XifNew;

            EPSsNew = EPSsNew ...
                + timeInc.*(TFb.*kepsh.*UTLbom.*Xhs ...
                + TFb2.*(kepsn1./MN).*UTLNH3.*Xn1s ...
                + TFb2.*(kepsn2./MN).*UTLNO2.*Xn2s ...
                + TFb.*and1bom.*kepsh.*Xhs ...
                + TFb.*and2bom.*kepsh.*Xhs ...
                - khydEPS.*EPSs ...
                - TFc.*kox1x2.*HOCl.*EPSs ...
                - TFc.*kox2x2.*OCl.*EPSs ...
                - TFc.*kox3x2.*NH2Cl.*EPSs ...
                + EPSdet ...
                - EPSads);
            EPSsNew(EPSsNew<0) = 0;

            EPSfNew = (EPSf ...
                + timeInc.*(TFb.*kepsh.*UTLbom.*Xhf ...
                + TFb2.*(kepsn1./MN).*UTLNH3.*Xn1f ...
                + TFb2.*(kepsn2./MN).*UTLNO2.*Xn2f ...
                + TFb.*and1bom.*kepsh.*Xhf ...
                + TFb.*and2bom.*kepsh.*Xhf ...
                - khydEPS.*EPSf ...
                - TFc.*kox1x2.*HOCl.*EPSf ...
                - TFc.*kox2x2.*OCl.*EPSf ...
                - TFc.*kox3x2.*NH2Cl.*EPSf ...
                - EPSdet...
                + EPSads)).*elementVolumesRatio;
            EPSfNew(EPSfNew<0) = 0;

            EPSNew = EPSsNew + EPSfNew;

            BOM1New = BOM1New ...
                + timeInc.*(-TFb.*UTLbom1.*Xh ...
                - TFb.*and1bom1.*Xh ...
                - TFb.*and2bom1.*Xh ...
                - TFc.*(kox1x1.*HOCl + kox2x1.*OCl + kox3x1.*NH2Cl).*BOM1);
            BOM1New(BOM1New<0) = 0;
            
            %Individual BOM1 loss terms
            rUTLbom1 = rUTLbom1 - TFb.*UTLbom1.*Xh;
            rand1bom1 = rand1bom1 - TFb.*and1bom1.*Xh;
            rand2bom1 = rand2bom1 - TFb.*and2bom1.*Xh;
            rkoxbom1 = rkoxbom1 - TFc.*(kox1x1.*HOCl + kox2x1.*OCl + kox3x1.*NH2Cl).*BOM1;

            BOM2New = BOM2New ...
                + timeInc.*(-TFb.*UTLbom2.*Xh ...
                - TFb.*and1bom2.*Xh ...
                - TFb.*and2bom2.*Xh ...
                + fd.*(kd1x1.*HOCl + kd2x1.*OCl + kd3x1.*NH2Cl).*Xhs ...
                + fd.*(kd1x3.*HOCl + kd2x3.*OCl + kd3x3.*NH2Cl).*Xn1s ...
                + fd.*(kd1x5.*HOCl + kd2x5.*OCl + kd3x5.*NH2Cl).*Xn2s ...
                + fd.*(kd1x2.*HOCl + kd2x2.*OCl + kd3x2.*NH2Cl).*Xhf ...
                + fd.*(kd1x4.*HOCl + kd2x4.*OCl + kd3x4.*NH2Cl).*Xn1f ...
                + fd.*(kd1x6.*HOCl + kd2x6.*OCl + kd3x6.*NH2Cl).*Xn2f ...
                - TFc.*(kox1x2.*HOCl + kox2x2.*OCl + kox3x2.*NH2Cl).*BOM2);
            BOM2New(BOM2New<0) = 0;
            
            %Individual BOM2 production loss terms
            rUTLbom2 = rUTLbom2 - TFb.*UTLbom2.*Xh;
            rand1bom2 = rand1bom2 - TFb.*and1bom2.*Xh;
            rand2bom2 = rand2bom2 - TFb.*and2bom2.*Xh;
            rkdbom2 = rkdbom2 + fd.*(kd1x1.*HOCl + kd2x1.*OCl + kd3x1.*NH2Cl).*Xhs ...
                + fd.*(kd1x3.*HOCl + kd2x3.*OCl + kd3x3.*NH2Cl).*Xn1s ...
                + fd.*(kd1x5.*HOCl + kd2x5.*OCl + kd3x5.*NH2Cl).*Xn2s ...
                + fd.*(kd1x2.*HOCl + kd2x2.*OCl + kd3x2.*NH2Cl).*Xhf ...
                + fd.*(kd1x4.*HOCl + kd2x4.*OCl + kd3x4.*NH2Cl).*Xn1f ...
                + fd.*(kd1x6.*HOCl + kd2x6.*OCl + kd3x6.*NH2Cl).*Xn2f;
            rkoxbom2 = rkoxbom2 - TFc.*(kox1x2.*HOCl + kox2x2.*OCl + kox3x2.*NH2Cl).*BOM2;

            UAPNew = UAPNew ...
                + timeInc.*(TFb.*kuaph.*UTLbom.*Xh ...
                + TFb2.*(kuapn1./MN).*UTLNH3.*Xn1 ...
                + TFb2.*(kuapn2./MN).*UTLNO2.*Xn2 ...
                - TFb.*UTLuap.*Xh ...
                + TFb.*kuaph.*and1bom1.*Xh ...
                + TFb.*kuaph.*and2bom1.*Xh ...
                + TFb.*kuaph.*and1bom2.*Xh ...
                + TFb.*kuaph.*and2bom2.*Xh ...
                - TFb.*and1uap.*Xh ...
                - TFb.*and2uap.*Xh ...
                - TFc.*(kox1x2.*HOCl + kox2x2.*OCl + kox3x2.*NH2Cl).*UAP);
            UAPNew(UAPNew<0) = 0;

            BAPNew = BAPNew ...
                + timeInc.*(khydEPS.*EPS ...
                - TFb.*UTLbap.*Xh ...
                - TFb.*and1bap.*Xh ...
                - TFb.*and2bap.*Xh ...
                - TFc.*(kox1x2.*HOCl + kox2x2.*OCl + kox3x2.*NH2Cl).*BAP);
            BAPNew(BAPNew<0) = 0;

            CtNH3New = CtNH3New ...
                + timeInc.*(-TFb2.*UTLNH3.*Xn1 ...
                - TFb2.*MN.*Yn1.*(1-kuapn1-kepsn1).*UTLNH3.*Xn1 ...
                - TFb2.*MN.*Yn2.*(1-kuapn2-kepsn2).*UTLNO2.*Xn2 ...
                - TFb2.*MN.*(kuapn1./MN).*UTLNH3.*Xn1 ...
                - TFb2.*MN.*(kuapn2./MN).*UTLNO2.*Xn2 ...
                - TFb2.*MN.*(kepsn1./MN).*UTLNH3.*Xn1 ...
                - TFb2.*MN.*(kepsn2./MN).*UTLNO2.*Xn2 ...
                + TFc.*kSC1.*NH2Cl.*NH2Cl.*(4./(diameters(elementPipes).*1000)) ...
                + MN.*fd.*(bh.*((O2./(KO2+O2))+(NO2./(KNO2+NO2)).*(KIDO./(KIDO.*O2))+(NO3./(KNO3+NO3)).*(KIDO./(KIDO.*O2))).*Xh + bn1.*(O2./(KO2+O2)).*Xn1 + bn2.*(O2./(KO2+O2)).*Xn2) ...
                + TFc.*KAD1.*NH2Cl.*NH2Cl ...
                - TFc.*KAD2.*NHCl2.*NH3.*H ...
                + TFc.*kox5.*NH2Cl.*NO2 ...
                + TFc.*KAD3.*(NH2Cl.*kCle./NH3).*NH2Cl ...
                + TFc.*11.*MN.*(kox3x1.*BOM1 + kox3x2.*BOM2 + kox3x2.*EPS + kox3x2.*BAP + kox3x2.*UAP + kox3x2.*Xi).*NH2Cl ...
                + TFc.*kcr3.*NH2Cl ...
                + TFb.*MN.*(1-Yp).*UTLbap.*Xh ...
                + TFb.*MN.*(1-Yp).*and1bap.*Xh ...
                + TFb.*MN.*(1-Yp).*and2bap.*Xh ...
                + TFb.*MN.*(1-Yp).*UTLuap.*Xh ...
                + TFb.*MN.*(1-Yp).*and1uap.*Xh ...
                + TFb.*MN.*(1-Yp).*and2uap.*Xh ...
                + TFb.*MN.*(1-Yh.*(1-kuaph-kepsh) - (kuaph+kepsh)).*UTLbom.*Xh ...
                + TFb.*MN.*(1-Yh.*(1-kuaph-kepsh) - (kuaph+kepsh)).*and1bom.*Xh ...
                + TFb.*MN.*(1-Yh.*(1-kuaph-kepsh) - (kuaph+kepsh)).*and2bom.*Xh);
            CtNH3New(CtNH3New<1E-14) = 1E-14;

            NO2New = NO2New ...
                + timeInc.*(-TFb2.*UTLNO2.*Xn2 ...
                - TFb.*MN.*(1-kuaph-kepsh-Yh.*(1-kuaph-kepsh)).*and1bom.*Xh ...
                - TFb.*MN.*(1-Yp).*and1smp.*Xh ...
                + TFb.*MN.*(1-kuaph-kepsh-Yh.*(1-kuaph-kepsh)).*and2bom.*Xh ...
                + TFb.*MN.*(1-Yp).*and2smp.*Xh ...
                - MN.*fd.*bh.*(NO2./(KNO2+NO2)).*(KIDO./(KIDO.*O2)).*Xh ...
                + MN.*fd.*bh.*(NO3./(KNO3+NO3)).*(KIDO./(KIDO.*O2)).*Xh ...
                + TFb2.*UTLNH3.*Xn1 ...
                + TFb2.*UTLNHOHCl.*Xn1 ...
                - TFc.*kox5.*NH2Cl.*NO2);
            NO2New(NO2New<0) = 0;

            NO3New = NO3New ...
                + timeInc.*(TFb2.*UTLNO2.*Xn2 ...
                - TFb.*MN.*(1-kuaph-kepsh-Yh.*(1-kuaph-kepsh)).*and2bom.*Xh ...
                - TFb.*MN.*(1-Yp).*and2smp.*Xh ...
                - MN.*fd.*bh.*(NO3./(KNO3+NO3)).*(KIDO./(KIDO.*O2)).*Xh ...
                + TFc.*kox5.*NH2Cl.*NO2);
            NO3New(NO3New<0) = 0;

            N2New = N2New ...
                + timeInc.*(TFb.*MN.*0.5.*(1-kuaph-kepsh-Yh.*(1-kuaph-kepsh)).*and1bom.*Xh ...
                + TFb.*MN.*0.5.*(1-Yp).*and1smp.*Xh ...
                + TFc.*KAD5.*NHOHCl ...
                + MN.*0.5.*fd.*bh.*(NO3./(KNO3+NO3)).*(KIDO./(KIDO.*O2)).*Xh);
            N2New(N2New<0) = 0;

            CtOClNew = CtOClNew ...
                + timeInc.*(-TFc.*2.0.*MC.*(kox1x1.*BOM1 + kox1x2.*BOM2 + kox1x2.*BAP + kox1x2.*EPS + kox1x2.*UAP + kox1x2.*Xi).*HOCl ...
                - TFc.*2.0.*MC.*(kox2x1.*BOM1 + kox2x2.*BOM2 + kox2x2.*BAP + kox2x2.*EPS + kox2x2.*UAP + kox2x2.*Xi).*OCl ...
                - TFc.*kox4.*HOCl.*NO2);
            CtOClNew(CtOClNew<0) = 0;

            NH2ClNew = NH2ClNew ...
                + timeInc.*(-TFc.*2.0.*KAD3.*(NH2Cl.*kCle./NH3).*NH2Cl ...
                - TFc.*2.0.*KAD1.*NH2Cl.*NH2Cl ...
                + TFc.*2.0.*KAD2.*NHCl2.*NH3.*H ...
                - TFc.*KAD5.*NHOHCl ...
                - TFc.*kox5.*NH2Cl.*NO2 ...
                - TFb2.*UTLNH2Cl.*Xn1 ...
                - TFc.*2.0.*MC.*(kox3x1.*BOM1 + kox3x2.*BOM2 + kox3x2.*BAP + kox3x2.*UAP + kox3x2.*Xi + kox3x2.*EPS).*NH2Cl ...
                - TFc.*kcr3.*NH2Cl ...
                - TFc.*2.*kSC1.*NH2Cl.*NH2Cl.*(4./(diameters(elementPipes).*1000)));
            NH2ClNew(NH2ClNew<0) = 0;
            
            % Individual NH2Cl loss terms
            rAD3 = rAD3 + timeInc.*(-TFc.*2.0.*KAD3.*(NH2Cl.*kCle./NH3).*NH2Cl);
            rAD1 = rAD1 + timeInc.*(- TFc.*2.0.*KAD1.*NH2Cl.*NH2Cl);
            rAD2 = rAD2 + timeInc.*(TFc.*2.0.*KAD2.*NHCl2.*NH3.*H);
            rAD5 = rAD5 + timeInc.*(- TFc.*KAD5.*NHOHCl);
            rox5 = rox5 + timeInc.*(- TFc.*kox5.*NH2Cl.*NO2);
            rUTLNH2Cl = rUTLNH2Cl + timeInc.*(- TFb2.*UTLNH2Cl.*Xn1);
            rkox = rkox + timeInc.*(- TFc.*2.0.*MC.*(kox3x1.*BOM1 + kox3x2.*BOM2 + kox3x2.*BAP + kox3x2.*UAP + kox3x2.*Xi + kox3x2.*EPS).*NH2Cl);
            rkcr = rkcr + timeInc.*(- TFc.*kcr3.*NH2Cl);
            rkSC = rkSC + timeInc.*(- TFc.*2.*kSC1.*NH2Cl.*NH2Cl.*(4./(diameters(elementPipes).*1000)));

            NHCl2New = NHCl2New ...
                + timeInc.*(TFc.*KAD1.*NH2Cl.*NH2Cl ...
                - TFc.*KAD2.*NHCl2.*NH3.*H ...
                + TFc.*KAD3.*(NH2Cl.*kCle./NH3).*NH2Cl ...
                - TFc.*KAD4.*NHCl2.*OH ...
                + TFc.*kSC1.*NH2Cl.*NH2Cl.*(4./(diameters(elementPipes).*1000)));
            NHCl2New(NHCl2New<0) = 0;

            NHOHClNew = NHOHClNew ...
                + timeInc.*(TFb2.*UTLNH2Cl.*Xn1 ...
                - TFb2.*UTLNHOHCl.*Xn1 ...
                + TFc.*KAD4.*NHCl2.*OH ...
                - TFc.*KAD5.*NHOHCl);
            NHOHClNew(NHOHClNew<0) = 0;

            ClNew = ClNew ...
                + timeInc.*(TFc.*KAD4.*NHCl2.*OH ...
                + TFc.*2.0.*KAD5.*NHOHCl ...
                + TFb2.*UTLNHOHCl.*Xn1 ...
                + TFc.*kox4.*HOCl.*NO2 ...
                + TFc.*kox5.*NH2Cl.*NO2 ...
                + TFc.*2.0.*MC.*(kox1x1.*BOM1 + kox1x2.*BOM2 + kox1x2.*BAP + kox1x2.*UAP + kox1x2.*Xi + kox1x2.*EPS).*HOCl ...
                + TFc.*2.0.*MC.*(kox2x1.*BOM1 + kox2x2.*BOM2 + kox2x2.*BAP + kox2x2.*UAP + kox2x2.*Xi + kox2x2.*EPS).*OCl ...
                + TFc.*2.0.*MC.*(kox3x1.*BOM1 + kox3x2.*BOM2 + kox3x2.*BAP + kox3x2.*UAP + kox3x2.*Xi + kox3x2.*EPS).*NH2Cl ...
                + TFc.*kcr1.*(4./(diameters(elementPipes).*1000)) ...
                + TFc.*kcr2.*(4./(diameters(elementPipes).*1000)) ...
                + TFc.*kcr3.*NH2Cl);
            ClNew(ClNew<0) = 0;

            O2New = O2New ...
                + timeInc.*(-TFb.*(1-kuaph-kepsh-Yh.*(1-kuaph-kepsh)).*UTLbom.*Xh ...
                - TFb.*(1-Yp).*UTLsmp.*Xh ...
                - TFb2.*(1-kuapn1-kepsn1-Yn1.*(1-kuapn1-kepsn1)).*UTLNH3.*Xn1 ...
                - TFb2.*(1-kuapn2-kepsn2-Yn2.*(1-kuapn2-kepsn2)).*UTLNO2.*Xn2 ...
                - fd.*(bh.*(O2./(KO2+O2)).*Xh ...
                + bn1.*(O2./(KO2+O2)).*Xn1+bn2.*(O2./(KO2+O2)).*Xn2) ...
                + kla.*(O2s-O2));
            O2New(O2New<0) = 0;
            O2New(O2New>O2s) = O2s;

            SFxn1sNew = TFb2.*Yn1.*(1 - kuapn1- kepsn1).*UTLNH3 ...
                - bn1.*(O2./(KO2+O2)) ...
                - kd3x3.*NH2Cl;

            SFxn1fNew = TFb2.*Yn1.*(1 - kuapn1- kepsn1).*UTLNH3 ...
                - bn1.*(O2./(KO2+O2)) ...
                - kd3x4.*NH2Cl;

            SFxn2sNew = TFb2.*Yn2.*(1 - kuapn2 - kepsn2).*UTLNO2 ...
                - bn2.*(O2./(KO2+O2)) ...
                - kd3x5.*NH2Cl;

            SFxn2fNew = TFb2.*Yn2.*(1 - kuapn2 - kepsn2).*UTLNO2 ...
                - bn2.*(O2./(KO2+O2)) ...
                - kd3x6.*NH2Cl;

            SFxhsNew = TFb.*Yh.*(1-kuaph-kepsh).*UTLbom ...
                + TFb.*Yh.*(1-kuaph-kepsh).*and1bom ...
                + TFb.*Yh.*(1-kuaph-kepsh).*and2bom ...
                + TFb.*Yp.*UTLsmp ...
                + TFb.*Yp.*and1smp ...
                + TFb.*Yp.*and2smp ...
                - bh.*((O2./(KO2+O2))+(NO2./(KNO2+NO2)).*(KIDO./(KIDO.*O2))+(NO3./(KNO3+NO3)).*(KIDO./(KIDO.*O2))) ...
                - kd3x1.*NH2Cl;

            SFxhfNew = TFb.*Yh.*(1-kuaph-kepsh).*UTLbom ...
                + TFb.*Yh.*(1-kuaph-kepsh).*and1bom ...
                + TFb.*Yh.*(1-kuaph-kepsh).*and2bom ...
                + TFb.*Yp.*UTLsmp ...
                + TFb.*Yp.*and1smp ...
                + TFb.*Yp.*and2smp ...
                - bh.*((O2./(KO2+O2))+(NO2./(KNO2+NO2)).*(KIDO./(KIDO.*O2))+(NO3./(KNO3+NO3)).*(KIDO./(KIDO.*O2))) ...
                - kd3x2.*NH2Cl;

            Time = TimeNew;
            Xhs = XhsNew;
            Xhf = XhfNew;
            Xh = XhNew;
            Xn1s = Xn1sNew;
            Xn1f = Xn1fNew;
            Xn1 = Xn1New;
            Xn2s = Xn2sNew;
            Xn2f = Xn2fNew;
            Xn2 = Xn2New;
            Xis = XisNew;
            Xif = XifNew;
            Xi = XiNew;
            EPSs = EPSsNew;
            EPSf = EPSfNew;
            EPS = EPSNew;
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
            SFxn1s = SFxn1sNew;
            SFxn1f = SFxn1fNew;
            SFxn2s = SFxn2sNew;
            SFxn2f = SFxn2fNew;
            SFxhs = SFxhsNew;
            SFxhf = SFxhfNew;

            % Check if solution is diverging by looking at one of the 
            % chemical species
            if any(isnan(NH2Cl))
                error('The solution is diverging. Try using a smaller time step.')
            end

        end

        % Final day hourly outputs
        if day == days
            CtOClO(:,hour) = CtOCl;
            NH2ClO(:,hour) = NH2Cl;
            NHCl2O(:,hour) = NHCl2;
            NHOHClO(:,hour) = NHOHCl;
            ClO(:,hour) = Cl;
            CtNH3O(:,hour) = CtNH3;
            NO2O(:,hour) = NO2;
            NO3O(:,hour) = NO3;
            N2O(:,hour) = N2;
            CO2O(:,hour) = CO2;
            O2O(:,hour) = O2;
            UAPO(:,hour) = UAP;
            BAPO(:,hour) = BAP;
            BOM1O(:,hour) = BOM1;
            BOM2O(:,hour) = BOM2;
            Xn1sO(:,hour) = Xn1s;
            Xn1fO(:,hour) = Xn1f;
            Xn2sO(:,hour) = Xn2s;
            Xn2fO(:,hour) = Xn2f;
            EPSsO(:,hour) = EPSs;
            EPSfO(:,hour) = EPSf;
            XhsO(:,hour) = Xhs;
            XhfO(:,hour) = Xhf;
            XisO(:,hour) = Xis;
            XifO(:,hour) = Xif;
            XbtO(:,hour) = Xhs + Xhf + Xn1s + Xn1f + Xn2s + Xn2f + EPSs + EPSf + Xis + Xif + UAP + BAP;
            XatO(:,hour) = Xhs + Xhf + Xn1s + Xn1f + Xn2s + Xn2f;
            NH4O(:,hour) = NH4;
            NH3O(:,hour) = NH3;
            TimeO(:,hour) = Time;
        end
    end
     
end

%Calculate total mass at end of simulation
totalC(2) = sum(MC.*(BOM1 + BOM2 + Xhs + Xhf + Xn1s + Xn1f + Xn2s + Xn2f + Xis + Xif + EPSs + EPSf + UAP + BAP + CO2).*elementVolumes.*1000);
totalCl(2) = sum((HOCl + OCl + NH2Cl + 2.*NHCl2 + NHOHCl + Cl).*elementVolumes.*1000);
totalN(2) = sum((NH2Cl + NHCl2 + NHOHCl + CtNH3 + NO2 + NO3 + 2.*N2 + MN.*(BOM1 + BOM2 + Xhs + Xhf + Xn1s + Xn1f + Xn2s + Xn2f + EPSs + EPSf + Xis + Xif + UAP + BAP)).*elementVolumes.*1000);

%Compute mass checks
checkC = ((totalC(2) - totalC(1) + lostC - externalC - addedC)/totalC(2))*100;
checkCl = ((totalCl(2) - totalCl(1) + lostCl - externalCl - addedCl)/totalCl(2))*100;
checkN = ((totalN(2) - totalN(1) + lostN - externalN - addedN)/totalN(2))*100;