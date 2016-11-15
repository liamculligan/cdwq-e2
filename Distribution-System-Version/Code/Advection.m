function Cnew = advection(C,nC)
    
    global elementReservoir
    global elementReservoirC1
    global elementReservoirC2
    global elementReservoirC3
    global reservoirCompartments
    global timeInc
    global Qin
    global Qout
    global hour
    global elementPipes
    global elementVolumes
    global cElementStartNodes
    global elementV
    
    Cnew(elementReservoir~=1,1) = C(elementReservoir~=1) ...
        + timeInc.*(((abs(Qin(elementPipes(elementReservoir~=1),hour))./elementVolumes(elementReservoir~=1)).*(nC(cElementStartNodes(elementV(elementReservoir~=1),hour)))) ...
        - (abs(Qout(elementPipes(elementReservoir~=1),hour))./elementVolumes(elementReservoir~=1)).*C(elementReservoir~=1));
    if reservoirCompartments == 1
        Cnew(elementReservoir==1 & elementReservoirC1==1,1) = C(elementReservoir==1 & elementReservoirC1==1) ...
            + timeInc.*(((abs(Qin(elementPipes(elementReservoir==1 & elementReservoirC1==1),hour))./elementVolumes(elementReservoir==1 & elementReservoirC1==1)).*(nC(cElementStartNodes(elementV(elementReservoir==1 & elementReservoirC1==1),hour))-C(elementReservoir==1 & elementReservoirC1==1))));
    elseif reservoirCompartments == 2
        Cnew(elementReservoir==1 & elementReservoirC1==1,1) = C(elementReservoir==1 & elementReservoirC1==1) ...
            + timeInc.*(((abs(Qin(elementPipes(elementReservoir==1 & elementReservoirC1==1),hour))./elementVolumes(elementReservoir==1 & elementReservoirC1==1)).*(nC(cElementStartNodes(elementV(elementReservoir==1 & elementReservoirC1==1),hour))-C(elementReservoir==1 & elementReservoirC1==1))) ...
            - (abs(Qout(elementPipes(elementReservoir==1 & elementReservoirC1==1),hour))./elementVolumes(elementReservoir==1 & elementReservoirC1==1)).*(C(elementReservoir==1 & elementReservoirC2==1)-C(elementReservoir==1 & elementReservoirC1==1)));
        Cnew(elementReservoir==1 & elementReservoirC2==1,1) = C(elementReservoir==1 & elementReservoirC2==1) ...
            + timeInc.*(((abs(Qin(elementPipes(elementReservoir==1 & elementReservoirC2==1),hour))./elementVolumes(elementReservoir==1 & elementReservoirC2==1)).*(C(elementReservoir==1 & elementReservoirC1==1)-C(elementReservoir==1 & elementReservoirC2==1))));
    else
        Cnew(elementReservoir==1 & elementReservoirC1==1,1) = C(elementReservoir==1 & elementReservoirC1==1) ...
            + timeInc.*(((abs(Qin(elementPipes(elementReservoir==1 & elementReservoirC1==1),hour))./elementVolumes(elementReservoir==1 & elementReservoirC1==1)).*(nC(cElementStartNodes(elementV(elementReservoir==1 & elementReservoirC1==1),hour))-C(elementReservoir==1 & elementReservoirC1==1))) ...
            - (abs(Qout(elementPipes(elementReservoir==1 & elementReservoirC1==1),hour))./elementVolumes(elementReservoir==1 & elementReservoirC1==1)).*(C(elementReservoir==1 & elementReservoirC2==1)-C(elementReservoir==1 & elementReservoirC1==1)));
        Cnew(elementReservoir==1 & elementReservoirC2==1,1) = C(elementReservoir==1 & elementReservoirC2==1) ...
            + timeInc.*(((abs(Qin(elementPipes(elementReservoir==1 & elementReservoirC2==1),hour))./elementVolumes(elementReservoir==1 & elementReservoirC2==1)).*(C(elementReservoir==1 & elementReservoirC1==1)-C(elementReservoir==1 & elementReservoirC2==1)+C(elementReservoir==1 & elementReservoirC3==1)-C(elementReservoir==1 & elementReservoirC2==1))) ...
            - (abs(Qout(elementPipes(elementReservoir==1 & elementReservoirC1==1),hour))./elementVolumes(elementReservoir==1 & elementReservoirC1==1)).*(C(elementReservoir==1 & elementReservoirC3==1)-C(elementReservoir==1 & elementReservoirC2==1)));
        Cnew(elementReservoir==1 & elementReservoirC3==1,1) = C(elementReservoir==1 & elementReservoirC3==1) ...
            + timeInc.*(((abs(Qin(elementPipes(elementReservoir==1 & elementReservoirC3==1),hour))./elementVolumes(elementReservoir==1 & elementReservoirC3==1)).*(C(cElementStartNodes(elementV(elementReservoir==1 & elementReservoirC2==1),hour))-C(elementReservoir==1 & elementReservoirC3==1))) ...
            - (abs(Qout(elementPipes(elementReservoir==1 & elementReservoirC1==1),hour))./elementVolumes(elementReservoir==1 & elementReservoirC1==1)).*(C(elementReservoir==1 & elementReservoirC2==1)-C(elementReservoir==1 & elementReservoirC3==1)));
    end
end