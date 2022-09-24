function [fitness, Ebonded, Enonbonded, ED] = findETotal(KTijk,Tijk,Teq,KRij,Rij,Req,Vmijkl,Fiijkl,Gijkl,Aij,Bij,qi,qj,chro)
EB=findEB(KRij,Rij,Req);
EA=findEA(KTijk,Tijk,Teq);
ED=findED(Vmijkl,Fiijkl,Gijkl);
Ebonded=EB+EA+ED;
if Rij<=3
%     ELJ=findELJ(Aij,Bij,Rij);
    ELJ=0;
    ES=findES(qi,qj,Rij,chro);
    Enonbonded=ELJ+ES;
else
    Enonbonded=0;
end
Etotal=Ebonded+Enonbonded;
fitness=Etotal;
end

