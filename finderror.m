function [ fitness, Ebonded, Enonbonded, ED, sq]=finderror(chromosomes,Dim,nop,reff)
    sum=0;    
    for i=1:length(reff)
        sum=sum+(reff(i)-chromosomes(i))^2;
    end
    sq=sqrt(sum);
    for i=1:nop
        chro(i,1:Dim)=chromosomes(1,(i-1)*Dim+1:i*Dim);
    end
    cod_bond=csvread('cod_bond.csv');
    cod_angle=csvread('cod_angle.csv');
    cod_dihedral=csvread('cod_dihedral.csv');
    charge=csvread('charge.csv');
    KTijk=cod_angle(:,4);
    Teq=cod_angle(:,5);
    ppp=length(KTijk);
    for i=1:ppp
        Tijk(i)=findAngle(chro(cod_angle(i,1),:),chro(cod_angle(i,2),:),chro(cod_angle(i,3),:));
    end
    KRij=cod_bond(:,3);
    Req=cod_bond(:,4);
    for i=1:size(cod_bond,1)
        Rij(i)=findDist(chro(cod_bond(i,1),:),chro(cod_bond(i,2),:));
    end
    Vmijkl=cod_dihedral(:,5);
    Fiijkl=cod_dihedral(:,6);
    
    for i=1:size(cod_dihedral,1)
        Gijkl(i)=findDihedral(chro(cod_dihedral(i,1),:),chro(cod_dihedral(i,2),:),chro(cod_dihedral(i,3),:),chro(cod_dihedral(i,4),:));
    end
    Aij=0;
    Bij=0;
    qi=charge;
    qj=charge;
    [fitness, Ebonded, Enonbonded, ED]=findETotal(KTijk,Tijk,Teq,KRij,Rij,Req,Vmijkl,Fiijkl,Gijkl,Aij,Bij,qi,qj,chro);
end