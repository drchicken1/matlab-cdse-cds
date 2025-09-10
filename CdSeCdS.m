clear all; 
nmtoao=18.8973;     %Bohrs in a nm
evtoH=0.0367493;    %Hartrees in an eV
lengsim=10000;      %Use a 10^4 point grid
bg=1.75*evtoH;      %bandgap of CdSe in eV converted to H
radius1=2.0*nmtoao; %core CdSe radius
radius2=2.0*nmtoao; %thickness of the CdS shell
delr=(radius1+radius2)/(lengsim/4); %Define the grid spacing

% Here are the parameters (potential energy offsets, effective masses,
% static dielectric constants. Values from the literature.
VeCdSe=0.0*evtoH;   %CdSe is the reference state
VeCdS=0.24*evtoH;   %Potential barrier for the e- in CdS
Veligands=3.0*evtoH;%High potential barrier for the exterior
VhCdSe=0.0*evtoH;   %CdSe is the reference state for holes
VhCdS=0.42*evtoH;   %Potential barrier for the h+ in CdS
Vhligands=3.0*evtoH;%High potential barrier for the exterior
mecdse=0.13;        %mass of the electron in CdSe
mecds=0.21;         %electron mass in CdS shell
meligands=1.0;      %normal electron mass
mhcdse=0.50;        %hole CdSe mass
mhcds=0.70;         %mass of CdS hole
mhligands=1.0;      %normal electron mass
dieleccdse=6.2;     %static dielectric of CdSe
dieleccds=8.9;      %static dielec of CdS
dielectol=2.38;     %static dielectric of toluene
buffersize=0.05*nmtoao; %switch function for smoothing the mass/potential sufaces

% These are the electric fields from the S-state holes, S-state electrons and
% the P-state electrons.
% They are initialized to 0 because we have no electrostatics in the first
% loop. There are three hole state wavefunctions, and four electron fields
% due to the fact that the P state produced 2 (a quadrupolar field)
Uh=zeros([lengsim 3]); Ue=zeros([lengsim 2]); UeP=zeros([lengsim 2]);
save Uh Uh; save Ue Ue; save UeP UeP; 

% This prepares the potential, mass, and electrostatic surfaces
% VXe,h are charge carrier potential surfaces in Hartrees
% me,h are masses, and eps is the dielectric surface.
for j=1:lengsim
    rax(j)=j*delr; 
    cdserdf(j)=rax(j)^2/radius1^2*0.5*(1-erf((rax(j)-radius1)/sqrt(2)/buffersize));
    temp=0.5*(1-erf((rax(j)-radius1)/sqrt(2)/buffersize));
    cdsrdf(j)=0.5*(1-erf((rax(j)-radius1-radius2)/sqrt(2)/buffersize))*(1-temp);
    temp=0.5*(1-erf((rax(j)-radius1-radius2)/sqrt(2)/buffersize));
    ligandsrdf(j)=(1-temp);
    % Electron potential surface:
    VXe(j)=VeCdSe*cdserdf(j)+(VeCdS-VeCdSe)*cdsrdf(j)+(Veligands-VeCdSe)*ligandsrdf(j);
    % Mass surface:
    me(j)=mecdse+(mecds-mecdse)*cdsrdf(j)+(meligands-mecdse)*ligandsrdf(j);
    % Hole potential surface:
    VXh(j)=VhCdSe*cdserdf(j)+(VhCdS-VhCdSe)*cdsrdf(j)+(Vhligands-VhCdSe)*ligandsrdf(j);
    % Hole mass surface:
    mh(j)=mhcdse+(mhcds-mhcdse)*cdsrdf(j)+(mhligands-mhcdse)*ligandsrdf(j);
    % Dielectric constant:
    eps(j)=dieleccdse+(dieleccds-dieleccdse)*cdsrdf(j)+(dielectol-dieleccdse)*ligandsrdf(j);
end

%The following fix an issue related to the FDM matrix size
me(lengsim+1:lengsim+2)=1.0; 
mh(lengsim+1:lengsim+2)=1.0; 
eps(lengsim+1)=eps(lengsim);

% Here we calculate the e,h energies with no electric fields.
% The functions return energy, which has to be recalled later on
% to address the fact that the Coulomb energy appearing 2X in the
% calculations (one for the e-, the other for h+).
initialEe(1)=CdSeCdS_E1SH1S(VXe,me,Uh,lengsim,delr); 
initialEe(2)=CdSeCdS_E1PH1S(VXe,me,Uh,lengsim,delr);
initialEe(3)=CdSeCdS_E2SH1S(VXe,me,Uh,lengsim,delr);
initialEh(1)=CdSeCdS_H1SE1S(VXh,mh,Ue,lengsim,delr); 
initialEh(2)=CdSeCdS_H1SE1P(VXh,mh,UeP,lengsim,delr);
initialEh(3)=CdSeCdS_H1SE2S(VXh,mh,Ue,lengsim,delr);
prevE=[initialEe initialEh];

% The electric fields are calculated from the initial wavefunctions, which
% then necessitate recalculation of the wavefunctions. The process loops
% until energy convergence is achieved.
errorfun=1;
while max(abs(errorfun))>0.0001
    Poisson_h(eps,lengsim,delr); 
    Poisson_e(eps,lengsim,delr);
    Poisson_eP(eps,lengsim,delr);
    load Uh; load Ue; load UeP;
    E1SH1S=CdSeCdS_E1SH1S(VXe,me,Uh,lengsim,delr); 
    E1PH1S=CdSeCdS_E1PH1S(VXe,me,Uh,lengsim,delr); 
    E2SH1S=CdSeCdS_E2SH1S(VXe,me,Uh,lengsim,delr); 
    EH1SE1S=CdSeCdS_H1SE1S(VXh,mh,Ue,lengsim,delr); 
    EH1SE1P=CdSeCdS_H1SE1P(VXh,mh,UeP,lengsim,delr); 
    EH1SE2S=CdSeCdS_H1SE2S(VXh,mh,Ue,lengsim,delr);
    %This is the percent change in energy for all states
    errorfun=([E1SH1S E1PH1S E2SH1S EH1SE1S EH1SE1P EH1SE2S]-prevE)./prevE;
    prevE=[E1SH1S E1PH1S E2SH1S EH1SE1S EH1SE1P EH1SE2S];
end
load WF_E1SH1S; load WF_E2SH1S; load WF_E1PH1S; 
load WF_H1SE1S; load WF_H1SE1P; load WF_H1SE2S; 

%Calculate the wave function overlaps between electron and holes
ehoverlap(1:3)=0; 
for j=1:lengsim
    ehoverlap(1)=ehoverlap(1)+rax(j)^2*WF_E1SH1S(j)*WF_H1SE1S(j)*delr;
    ehoverlap(2)=ehoverlap(2)+rax(j)^2*WF_E1PH1S(j)*WF_H1SE1P(j)*delr;
    ehoverlap(3)=ehoverlap(3)+rax(j)^2*WF_E2SH1S(j)*WF_H1SE2S(j)*delr;
end

% Exciton energies equal to the bandgap plus confinement plus 
% electrostatic energies of the e- and h+, followed by a correction factor
% for the fact that both the e- and h+ each have the same quantity of 
% electric potential (each should just have half that energy).
EnergyE1SH1S=bg+E1SH1S+EH1SE1S+(initialEe(1)-E1SH1S+initialEh(1)-EH1SE1S)/2;
EnergyE1PH1S=bg+E1PH1S+EH1SE1P+(initialEe(2)-E1PH1S+initialEh(2)-EH1SE1P)/2;
EnergyE2SH1S=bg+E2SH1S+EH1SE2S+(initialEe(3)-E2SH1S+initialEh(3)-EH1SE2S)/2;
%Convert to eV
EnergyE1SH1S=EnergyE1SH1S/evtoH;
EnergyE1PH1S=EnergyE1PH1S/evtoH;
EnergyE2SH1S=EnergyE2SH1S/evtoH;

% Biexciton state energy calculations, aslo in eV
% Initialization of variables
load Ue; load Uh; load UeP; load UePQ;
E1s1s=0; E1sh1s=0; E1sh1p=0; E1sh2s=0; 
E2s2s=0; E2sh1s=0; E2sh1p=0; E2sh2s=0;
E1p1pzz=0; E1p1pzx=0; E1ph1s=0; E1ph1p=0; E1ph2s=0; 
E1s1p=0; E1s2s=0; E1p2s=0;
Eh1sh1s=0; Eh1sh1p=0; Eh1sh2s=0; Eh1ph1p=0; Eh1ph2s=0; Eh2sh2s=0;
% Coulomb interactions. There are many of these because there are a large
% number of permutations of interactions for a biexciton.
for j=1:lengsim
    E1s1s=E1s1s    -    rax(j)^2*WF_E1SH1S(j)*Ue(j,1)*WF_E1SH1S(j)*delr/evtoH;
    E1sh1s=E1sh1s  -    rax(j)^2*WF_E1SH1S(j)*Uh(j,1)*WF_E1SH1S(j)*delr/evtoH;
    E1sh1p=E1sh1p  -0.5*rax(j)^2*WF_E1SH1S(j)*Uh(j,2)*WF_E1SH1S(j)*delr/evtoH;
    E1sh1p=E1sh1p  +0.5*rax(j)^2*WF_H1SE1P(j)*Ue(j,1)*WF_H1SE1P(j)*delr/evtoH;
    E1sh2s=E1sh2s  -0.5*rax(j)^2*WF_E1SH1S(j)*Uh(j,3)*WF_E1SH1S(j)*delr/evtoH;
    E1sh2s=E1sh2s  +0.5*rax(j)^2*WF_H1SE2S(j)*Ue(j,1)*WF_H1SE2S(j)*delr/evtoH;

    E1p1pzz=E1p1pzz-    rax(j)^2*WF_E1PH1S(j)*UeP(j)*WF_E1PH1S(j)*delr/evtoH;
    E1p1pzz=E1p1pzz-2/5*rax(j)^2*WF_E1PH1S(j)*UePQ(j)*WF_E1PH1S(j)*delr/evtoH;
    E1p1pzx=E1p1pzx-    rax(j)^2*WF_E1PH1S(j)*UeP(j)*WF_E1PH1S(j)*delr/evtoH;
    E1p1pzx=E1p1pzx+1/5*rax(j)^2*WF_E1PH1S(j)*UePQ(j)*WF_E1PH1S(j)*delr/evtoH;
    E1ph1s=E1ph1s  -0.5*rax(j)^2*WF_E1PH1S(j)*Uh(j,1)*WF_E1PH1S(j)*delr/evtoH;
    E1ph1s=E1ph1s  +0.5*rax(j)^2*WF_H1SE1S(j)*UeP(j)*WF_H1SE1S(j)*delr/evtoH;
    E1ph1p=E1ph1p  -    rax(j)^2*WF_E1PH1S(j)*Uh(j,2)*WF_E1PH1S(j)*delr/evtoH;
    E1ph2s=E1ph2s  -0.5*rax(j)^2*WF_E1PH1S(j)*Uh(j,3)*WF_E1PH1S(j)*delr/evtoH;
    E1ph2s=E1ph2s  +0.5*rax(j)^2*WF_H1SE2S(j)*UeP(j)*WF_H1SE2S(j)*delr/evtoH;

    E2s2s=E2s2s    -    rax(j)^2*WF_E2SH1S(j)*Ue(j,2)*WF_E2SH1S(j)*delr/evtoH;
    E2sh1s=E2sh1s  -0.5*rax(j)^2*WF_E2SH1S(j)*Uh(j,1)*WF_E2SH1S(j)*delr/evtoH;
    E2sh1s=E2sh1s  +0.5*rax(j)^2*WF_H1SE1S(j)*Ue(j,2)*WF_H1SE1S(j)*delr/evtoH;
    E2sh1p=E2sh1p  -0.5*rax(j)^2*WF_E2SH1S(j)*Uh(j,2)*WF_E2SH1S(j)*delr/evtoH;
    E2sh1p=E2sh1p  +0.5*rax(j)^2*WF_H1SE1P(j)*Ue(j,2)*WF_H1SE1P(j)*delr/evtoH;
    E2sh2s=E2sh2s  -    rax(j)^2*WF_E2SH1S(j)*Uh(j,3)*WF_E2SH1S(j)*delr/evtoH;
    
    E1s1p=E1s1p-0.5*    rax(j)^2*WF_E1SH1S(j)*UeP(j)*WF_E1SH1S(j)*delr/evtoH;
    E1s1p=E1s1p-0.5*    rax(j)^2*WF_E1PH1S(j)*Ue(j,1)*WF_E1PH1S(j)*delr/evtoH;
    E1s2s=E1s2s-0.5*    rax(j)^2*WF_E1SH1S(j)*Ue(j,2)*WF_E1SH1S(j)*delr/evtoH;
    E1s2s=E1s2s-0.5*    rax(j)^2*WF_E2SH1S(j)*Ue(j,1)*WF_E2SH1S(j)*delr/evtoH;
    E1p2s=E1p2s-0.5*    rax(j)^2*WF_E1PH1S(j)*Ue(j,2)*WF_E1PH1S(j)*delr/evtoH;
    E1p2s=E1p2s-0.5*    rax(j)^2*WF_E2SH1S(j)*UeP(j) *WF_E2SH1S(j)*delr/evtoH;
    
    Eh1sh1s=Eh1sh1s    +rax(j)^2*WF_H1SE1S(j)*Uh(j,1)*WF_H1SE1S(j)*delr/evtoH;
    Eh1sh1p=Eh1sh1p+0.5*rax(j)^2*WF_H1SE1S(j)*Uh(j,2)*WF_H1SE1S(j)*delr/evtoH;
    Eh1sh1p=Eh1sh1p+0.5*rax(j)^2*WF_H1SE1P(j)*Uh(j,1)*WF_H1SE1P(j)*delr/evtoH;
    Eh1sh2s=Eh1sh2s+0.5*rax(j)^2*WF_H1SE1S(j)*Uh(j,3)*WF_H1SE1S(j)*delr/evtoH;
    Eh1sh2s=Eh1sh2s+0.5*rax(j)^2*WF_H1SE2S(j)*Uh(j,1)*WF_H1SE2S(j)*delr/evtoH;
    Eh1ph1p=Eh1ph1p    +rax(j)^2*WF_H1SE1P(j)*Uh(j,2)*WF_H1SE1P(j)*delr/evtoH;
    Eh1ph2s=Eh1ph2s+0.5*rax(j)^2*WF_H1SE1P(j)*Uh(j,3)*WF_H1SE1P(j)*delr/evtoH;
    Eh1ph2s=Eh1ph2s+0.5*rax(j)^2*WF_H1SE2S(j)*Uh(j,2)*WF_H1SE2S(j)*delr/evtoH;
    Eh2sh2s=Eh2sh2s    +rax(j)^2*WF_H1SE2S(j)*Uh(j,3)*WF_H1SE2S(j)*delr/evtoH;
end

%Energies of the states plus electrostatic contributions:
BE1S1S=2*EnergyE1SH1S+E1s1s+Eh1sh1s+2*E1sh1s;
BE1PE1Pzz=2*EnergyE1PH1S+E1p1pzz+Eh1ph1p+2*E1ph1p;
BE1PE1Pzx=2*EnergyE1PH1S+E1p1pzx+Eh1ph1p+2*E1ph1p;
BE2SE2S=2*EnergyE2SH1S+E2s2s+Eh2sh2s+2*E2sh2s;

BE1SE2S=EnergyE1SH1S+EnergyE2SH1S+E1s2s+Eh1sh2s+1/2*E1sh2s+E2sh1s;
BE1SE1P=EnergyE1SH1S+EnergyE1PH1S+E1s1p+Eh1sh1p+1/2*E1sh1p+E1ph1s;
BE1PE2S=EnergyE1PH1S+EnergyE2SH1S+E1p2s+Eh1ph2s+1/2*E1ph2s+E2sh1p;

%Convert the wavefunctions back to nm for plotting:
delr_nm=delr/nmtoao;
norm(1:4)=0;
for j=1:lengsim
    rax_nm(j)=delr_nm*j;
    norm(1)=norm(1)+delr_nm*rax_nm(j)^2*WF_E1SH1S(j)^2;
    norm(2)=norm(2)+delr_nm*rax_nm(j)^2*WF_E1PH1S(j)^2;
    norm(3)=norm(3)+delr_nm*rax_nm(j)^2*WF_E2SH1S(j)^2;
    norm(4)=norm(4)+delr_nm*rax_nm(j)^2*WF_H1SE1S(j)^2;
end
%These are the RDF's
for j=1:lengsim
    WF_E1SH1S_2(j)=rax_nm(j)^2*WF_E1SH1S(j)^2/norm(1); 
    WF_E1PH1S_2(j)=rax_nm(j)^2*WF_E1PH1S(j)^2/norm(2); 
    WF_E2SH1S_2(j)=rax_nm(j)^2*WF_E2SH1S(j)^2/norm(3); 
    WF_H1SE1S_2(j)=rax_nm(j)^2*WF_H1SE1S(j)^2/norm(4);
end
close all;
plot(rax_nm,WF_E1SH1S_2,'b')
hold on
plot(rax_nm,WF_H1SE1S_2,'y--')
plot(rax_nm,WF_E2SH1S_2,'r')
plot(rax_nm,WF_E1PH1S_2,'g')
axis([0 5 0 1.1])
set(gca,'TickDir','out');
dim = [.2 .5 .3 .3]; 
str = 'Blue: E1S|H1S. Green: E1P|H1S. Red: E2S|H1S';
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on'); 
xlabel('r(nm)');ylabel('r^2*Y(r)^2');
fprintf("Exciton energy 1Se|1Sh %.3f eV\n",EnergyE1SH1S);
fprintf("Exciton energy 1Pe|1Sh %.3f eV\n",EnergyE1PH1S);
fprintf("Exciton energy 2Se|1Sh %.3f eV\n",EnergyE2SH1S);
fprintf("Bixciton energy 1Se1Se|1Sh1Sh %.3f eV Binding energy %.3f eV\n",BE1S1S,E1s1s+Eh1sh1s+2*E1sh1s);
fprintf("Bixciton energy 1Pe(z)1Pe(z)|1Sh1Sh %.3f eV. Binding energy %.3f eV\n",BE1PE1Pzz,E1p1pzz+Eh1ph1p+2*E1ph1p);
fprintf("Bixciton energy 1Pe(z)1Pe(x,y)|1Sh1Sh %.3f eV. Binding energy %.3f eV\n",BE1PE1Pzx,E1p1pzx+Eh1ph1p+2*E1ph1p);
fprintf("Bixciton energy 2Se2Se|1Sh1Sh %.3f eV. Binding energy %.3f eV\n",BE2SE2S,E2s2s+Eh2sh2s+2*E2sh2s);
fprintf("Bixciton energy 2Se1Se|1Sh1Sh %.3f eV. Binding energy %.3f eV\n",BE1SE2S,E1s2s+Eh1sh2s+1/2*E1sh2s+E2sh1s);
fprintf("Bixciton energy 1Pe1Se|1Sh1Sh %.3f eV. Binding energy %.3f eV \n",BE1SE1P,E1s1p+Eh1sh1p+1/2*E1sh1p+E1ph1s);
fprintf("Bixciton energy 1Pe2Se|1Sh1Sh %.3f eV. Binding energy %.3f eV\n",BE1PE2S,E1p2s+Eh1ph2s+1/2*E1ph2s+E2sh1p);



function [E1S]=CdSeCdS_E1SH1S(VXe,me,Uh,lengsim,delr) 
    % The following uses a 3 point stencil to solve the SE for a 1S state.
    
    l=0; %To make it clear this is a 1S state, angular momentum is 0.
    delrsq=delr*delr;
    
    % Initialization of the Hamiltonian FDM matrix. h is the traditional form, 
    % h2 takes into account the corrections necessary due to the varialbe mass 
    % between core and shell.
    h=zeros([lengsim lengsim+2]);
    h2=zeros([lengsim lengsim+2]);
    
    % Implementing boundary conditions into FDM. Below sets the wavefunction
    % derivative equal 0 at the origin.
    r=delr;
    h(1,1)=(1/2/me(1))*(2/r)*1/delr; 
    h(1,2)=0;
    h(1,3)=(1/2/me(1))*-(2/r)*1/delr;
    % Given the need to implement a boundary condition, h2 must be 0.
    h2(1,2)=0; 
    h2(1,1)=0; 
    
    % The 3point stencil structure of the FDM requires a separate
    % initialization for the 2nd row compared to the remaining matrix.
    for j=2:2
        r=j*delr;
        h(j,j)=(-1/2/me(j))*(-5/2/delr/delr)-Uh(j,1)+VXe(j)+(1/2/me(2))*l*(l+1)/r^2; 
        h(j,j+1)=(-1/2/me(j))*(4/3/delrsq+(2/r)*2/3/delr);  
        h(j,j+2)=(-1/2/me(j))*(-1/12/delrsq-(2/r)*1/12/delr); 
        h(j,j-1)=(-1/2/me(j))*(4/3/delrsq-(2/r)*2/3/delr);
        % Variable mass term:
        h2(j,j)=0; 
        h2(j,j+1)=-(1/2)*(-1/me(j)/me(j))*((me(j+1)-me(j-1))/2/delr)*2/3/delr;  
        h2(j,j+2)=-(1/2)*(-1/me(j)/me(j))*((me(j+1)-me(j-1))/2/delr)*(-1/12/delr);
        h2(j,j-1)=-(1/2)*(-1/me(j)/me(j))*((me(j+1)-me(j-1))/2/delr)*(-2/3/delr);
    end
    % The remaining FDM matrix:
    for j=3:lengsim
        r=j*delr;
        h(j,j)=(-1/2/me(j))*(-5/2/delrsq)-Uh(j,1)+VXe(j)+(1/2/me(j))*l*(l+1)/r^2; 
        h(j,j+1)=(-1/2/me(j))*(4/3/delrsq+(2/r)*2/3/delr);  
        h(j,j+2)=(-1/2/me(j))*(-1/12/delrsq-(2/r)*1/12/delr);
        h(j,j-1)=(-1/2/me(j))*(4/3/delrsq-(2/r)*2/3/delr);
        h(j,j-2)=(-1/2/me(j))*(-1/12/delrsq+(2/r)*1/12/delr);
        h2(j,j)=0; 
        h2(j,j+1)=-(1/2)*(-1/me(j)/me(j))*(-1/12*me(j+2)+2/3*me(j+1)-2/3*me(j-1)+1/12*me(j-2))/delr*(2/3/delr);  
        h2(j,j+2)=-(1/2)*(-1/me(j)/me(j))*(-1/12*me(j+2)+2/3*me(j+1)-2/3*me(j-1)+1/12*me(j-2))/delr*(-1/12/delr);
        h2(j,j-1)=-(1/2)*(-1/me(j)/me(j))*(-1/12*me(j+2)+2/3*me(j+1)-2/3*me(j-1)+1/12*me(j-2))/delr*(-2/3/delr);
        h2(j,j-2)=-(1/2)*(-1/me(j)/me(j))*(-1/12*me(j+2)+2/3*me(j+1)-2/3*me(j-1)+1/12*me(j-2))/delr*(1/12/delr);
    end
    
    % Sums the "normal" and variable mass correctd Hamiltonian matrices:
    q=h2(1:lengsim,1:lengsim)+h(1:lengsim,1:lengsim);
    
    % Sparse matrix eigenvale and eigenvector calculator.
    [V, D] = eigs(sparse(q), 15, 0);
    
    % MATLAB doesn't necessary present the eigenvalues and vectors in order. 
    % This part of the code addresses that:
    [d_notserp, idx] = sort(diag(D), 'ascend');
    Eens = diag(d_notserp);
    Vs = V(:, idx);
    
    %Lowest energy eigenvalue and eigenvector:
    E1S=Eens(1,1);
    WF_E1SH1S=Vs(:,1)/Vs(1,1);%this is done to flip the WF over if it is negative
    
    %Normalization of the 1S state:
    norm1=0.0;
    for j=1:lengsim 
        norm1=norm1+(j*delr)^2*WF_E1SH1S(j)^2*delr;
    end 
    WF_E1SH1S=WF_E1SH1S/sqrt(norm1);
    save WF_E1SH1S WF_E1SH1S;
end


function [E1P]=CdSeCdS_E1PH1S(VXe,me,Uh,lengsim,delr) 
    % The following uses a 3 point stencil to solve the SE for a 1S state.
    
    l=1;    %This is 1P state
    delrsq=delr*delr;
    
    % Initialization of the Hamiltonian FDM matrix. h is the traditional form, 
    % h2 takes into account the corrections necessary due to the varialbe mass 
    % between core and shell.
    h=zeros([lengsim lengsim+2]);h2=zeros([lengsim lengsim+2]);
    
    % Unlike an S-state there is no need to initialze the derivative to 0 for a
    % p-symmetry state. Thus we begin with initializing the Hamiltonian.
    % The 3point stencil structure of the FDM requires a separate
    % initialization for the 2nd row compared to the remaining matrix.
    r=delr; 
    h(1,1)=(-1/2/me(1))*(-5/2/delr/delr)-Uh(1)+VXe(1)+(1/2/me(1))*l*(l+1)/r^2;
    h(1,2)=(-1/2/me(1))*(4/3/delrsq+(2/r)*2/3/delr);  
    h(1,3)=(-1/2/me(1))*(-1/12/delrsq-(2/r)*1/12/delr);
    h2(1,1)=0; 
    h2(1,2)=-(1/2)*(-1/me(1)/me(1))*((me(2)-me(1))/delr)*2/3/delr;  
    h2(1,3)=-(1/2)*(-1/me(1)/me(1))*((me(2)-me(1))/delr)*(-1/12/delr);
    for j=2:2
        r=j*delr;
        h(j,j)=(-1/2/me(j))*(-5/2/delr/delr)-Uh(j,2)+VXe(j)+(1/2/me(2))*l*(l+1)/r^2; 
        h(j,j+1)=(-1/2/me(j))*(4/3/delrsq+(2/r)*2/3/delr);  
        h(j,j+2)=(-1/2/me(j))*(-1/12/delrsq-(2/r)*1/12/delr); 
        h(j,j-1)=(-1/2/me(j))*(4/3/delrsq-(2/r)*2/3/delr);
        h2(j,j)=0; 
        h2(j,j+1)=-(1/2)*(-1/me(j)/me(j))*((me(j+1)-me(j-1))/2/delr)*2/3/delr;  
        h2(j,j+2)=-(1/2)*(-1/me(j)/me(j))*((me(j+1)-me(j-1))/2/delr)*(-1/12/delr);
        h2(j,j-1)=-(1/2)*(-1/me(j)/me(j))*((me(j+1)-me(j-1))/2/delr)*(-2/3/delr);
    end
    % The remaining FDM matrix:
    for j=3:lengsim
        r=j*delr;
        h(j,j)=(-1/2/me(j))*(-5/2/delrsq)-Uh(j,2)+VXe(j)+(1/2/me(j))*l*(l+1)/r^2; 
        h(j,j+1)=(-1/2/me(j))*(4/3/delrsq+(2/r)*2/3/delr);  
        h(j,j+2)=(-1/2/me(j))*(-1/12/delrsq-(2/r)*1/12/delr);
        h(j,j-1)=(-1/2/me(j))*(4/3/delrsq-(2/r)*2/3/delr);
        h(j,j-2)=(-1/2/me(j))*(-1/12/delrsq+(2/r)*1/12/delr);
        h2(j,j)=0; 
        h2(j,j+1)=-(1/2)*(-1/me(j)/me(j))*(-1/12*me(j+2)+2/3*me(j+1)-2/3*me(j-1)+1/12*me(j-2))/delr*(2/3/delr);  
        h2(j,j+2)=-(1/2)*(-1/me(j)/me(j))*(-1/12*me(j+2)+2/3*me(j+1)-2/3*me(j-1)+1/12*me(j-2))/delr*(-1/12/delr);
        h2(j,j-1)=-(1/2)*(-1/me(j)/me(j))*(-1/12*me(j+2)+2/3*me(j+1)-2/3*me(j-1)+1/12*me(j-2))/delr*(-2/3/delr);
        h2(j,j-2)=-(1/2)*(-1/me(j)/me(j))*(-1/12*me(j+2)+2/3*me(j+1)-2/3*me(j-1)+1/12*me(j-2))/delr*(1/12/delr);
    end
    
    % Sums the "normal" and variable mass correctd Hamiltonian matrices:
    q=h2(1:lengsim,1:lengsim)+h(1:lengsim,1:lengsim);
    
    % Sparse matrix eigenvale and eigenvector calculator.
    [V, D] = eigs(sparse(q), 15, 0);
    
    % MATLAB doesn't necessary present the eigenvalues and vectors in order. 
    % This part of the code addresses that:
    [d_notserp, idx] = sort(diag(D), 'ascend');
    Eensp = diag(d_notserp);
    Vsp = V(:, idx);
    
    %Lowest energy eigenvalue and eigenvector:
    E1P=Eensp(1,1);
    WF_E1PH1S=Vsp(:,1)/Vsp(1,1); %this is done to flip the WF over if it is negative
    norm1=0.0;
    for j=1:lengsim 
        norm1=norm1+(j*delr)^2*WF_E1PH1S(j)^2*delr;
    end 
    WF_E1PH1S=WF_E1PH1S/sqrt(norm1);
    save WF_E1PH1S WF_E1PH1S;
end


function [E2S]=CdSeCdS_E2SH1S(VXe,me,Uh,lengsim,delr) 
    % The following uses a 3 point stencil to solve the SE for a 2S state.
    
    l=0; %To make it clear this is a 1S state, angular momentum is 0.
    delrsq=delr*delr;
    
    % Initialization of the Hamiltonian FDM matrix. h is the traditional form, 
    % h2 takes into account the corrections necessary due to the varialbe mass 
    % between core and shell.
    h=zeros([lengsim lengsim+2]);
    h2=zeros([lengsim lengsim+2]);
    
    % Implementing boundary conditions into FDM. Below sets the wavefunction
    % derivative equal 0 at the origin.
    r=delr; 
    h(1,1)=(1/2/me(1))*(2/r)*1/delr; %This sets the derivative equal 0 at the origin
    h(1,2)=0;
    h(1,3)=(1/2/me(1))*-(2/r)*1/delr;
    % Given the need to implement a boundary condition, h2 must be 0.
    h2(1,2)=0; 
    h2(1,1)=0; 
    
    % The 3point stencil structure of the FDM requires a separate
    % initialization for the 2nd row compared to the remaining matrix.
    for j=2:2
        r=j*delr;
        h(j,j)=(-1/2/me(j))*(-5/2/delr/delr)-Uh(j,3)+VXe(j)+(1/2/me(2))*l*(l+1)/r^2; 
        h(j,j+1)=(-1/2/me(j))*(4/3/delrsq+(2/r)*2/3/delr);  
        h(j,j+2)=(-1/2/me(j))*(-1/12/delrsq-(2/r)*1/12/delr); 
        h(j,j-1)=(-1/2/me(j))*(4/3/delrsq-(2/r)*2/3/delr);
        h2(j,j)=0; 
        h2(j,j+1)=-(1/2)*(-1/me(j)/me(j))*((me(j+1)-me(j-1))/2/delr)*2/3/delr;  
        h2(j,j+2)=-(1/2)*(-1/me(j)/me(j))*((me(j+1)-me(j-1))/2/delr)*(-1/12/delr);
        h2(j,j-1)=-(1/2)*(-1/me(j)/me(j))*((me(j+1)-me(j-1))/2/delr)*(-2/3/delr);
    end
    %The remaining FDM matrix:
    for j=3:lengsim
        r=j*delr;
        h(j,j)=(-1/2/me(j))*(-5/2/delrsq)-Uh(j,3)+VXe(j)+(1/2/me(j))*l*(l+1)/r^2; 
        h(j,j+1)=(-1/2/me(j))*(4/3/delrsq+(2/r)*2/3/delr);  
        h(j,j+2)=(-1/2/me(j))*(-1/12/delrsq-(2/r)*1/12/delr);
        h(j,j-1)=(-1/2/me(j))*(4/3/delrsq-(2/r)*2/3/delr);
        h(j,j-2)=(-1/2/me(j))*(-1/12/delrsq+(2/r)*1/12/delr);
        h2(j,j)=0; 
        h2(j,j+1)=-(1/2)*(-1/me(j)/me(j))*(-1/12*me(j+2)+2/3*me(j+1)-2/3*me(j-1)+1/12*me(j-2))/delr*(2/3/delr);  
        h2(j,j+2)=-(1/2)*(-1/me(j)/me(j))*(-1/12*me(j+2)+2/3*me(j+1)-2/3*me(j-1)+1/12*me(j-2))/delr*(-1/12/delr);
        h2(j,j-1)=-(1/2)*(-1/me(j)/me(j))*(-1/12*me(j+2)+2/3*me(j+1)-2/3*me(j-1)+1/12*me(j-2))/delr*(-2/3/delr);
        h2(j,j-2)=-(1/2)*(-1/me(j)/me(j))*(-1/12*me(j+2)+2/3*me(j+1)-2/3*me(j-1)+1/12*me(j-2))/delr*(1/12/delr);
    end
    
    % Sums the "normal" and variable mass correctd Hamiltonian matrices:
    q=h2(1:lengsim,1:lengsim)+h(1:lengsim,1:lengsim);
    
    % Sparse matrix eigenvale and eigenvector calculator.
    [V, D] = eigs(sparse(q), 15, 0);
    
    %This part of the code makes sure the eigenvalues and vectors are in order
    [d_notserp, idx] = sort(diag(D), 'ascend');
    Eens = diag(d_notserp);
    Vs = V(:, idx);
    
    % 1st excited state eigenvalue and eigenvector:
    E2S=Eens(2,2);
    WF_E2SH1S=Vs(:,2)/Vs(1,2); %this is done to flip the WF over if it is negative
    
    %Normalization of the 2S state:
    norm1=0.0;
    for j=1:lengsim 
        norm1=norm1+(j*delr)^2*WF_E2SH1S(j)^2*delr;
    end 
    WF_E2SH1S=WF_E2SH1S/sqrt(norm1);
    save WF_E2SH1S WF_E2SH1S;
end


function [EH1S]=CdSeCdS_H1SE1S(VXh,mh,Ue,lengsim,delr) 
    % The following uses a 3 point stencil to solve the SE for a 1S hole state.
    %To make it clear this is a 1S state, angular momentum is 0.
    l=0; 
    delrsq=delr*delr;
    
    % Initialization of the Hamiltonian FDM matrix. h is the traditional form, 
    % h2 takes into account the corrections necessary due to the varialbe mass 
    % between core and shell.
    h=zeros([lengsim lengsim+2]);
    h2=zeros([lengsim lengsim+2]);
    
    % Implementing boundary conditions into FDM. Below sets the wavefunction
    % derivative equal 0 at the origin.
    r=delr; 
    h(1,1)=(1/2/mh(1))*(2/r)*1/delr; 
    h(1,2)=0;
    h(1,3)=(1/2/mh(1))*-(2/r)*1/delr;
    
    % Given the need to implement a boundary condition, h2 must be 0.
    h2(1,2)=0; 
    h2(1,1)=0; 
    
    % The 3point stencil structure of the FDM requires a separate
    % initialization for the 2nd row compared to the remaining matrix.
    for j=2:2
        r=j*delr;
        h(j,j)=(-1/2/mh(j))*(-5/2/delr/delr)+Ue(j,1)+VXh(j)+1/2/mh(2)*l*(l+1)/r^2; 
        h(j,j+1)=(-1/2/mh(j))*(4/3/delrsq+(2/r)*2/3/delr);  
        h(j,j+2)=(-1/2/mh(j))*(-1/12/delrsq-(2/r)*1/12/delr); 
        h(j,j-1)=(-1/2/mh(j))*(4/3/delrsq-(2/r)*2/3/delr);
        h2(j,j)=0; 
        h2(j,j+1)=-(1/2)*(-1/mh(j)/mh(j))*((mh(j+1)-mh(j-1))/2/delr)*2/3/delr;  
        h2(j,j+2)=-(1/2)*(-1/mh(j)/mh(j))*((mh(j+1)-mh(j-1))/2/delr)*(-1/12/delr);
        h2(j,j-1)=-(1/2)*(-1/mh(j)/mh(j))*((mh(j+1)-mh(j-1))/2/delr)*(-2/3/delr);
    end
    % The remaining FDM matrix:
    for j=3:lengsim
        r=j*delr;
        h(j,j)=(-1/2/mh(j))*(-5/2/delrsq)+Ue(j,1)+VXh(j)+1/2/mh(j)*l*(l+1)/r^2; 
        h(j,j+1)=(-1/2/mh(j))*(4/3/delrsq+(2/r)*2/3/delr);  
        h(j,j+2)=(-1/2/mh(j))*(-1/12/delrsq-(2/r)*1/12/delr);
        h(j,j-1)=(-1/2/mh(j))*(4/3/delrsq-(2/r)*2/3/delr);
        h(j,j-2)=(-1/2/mh(j))*(-1/12/delrsq+(2/r)*1/12/delr);
        h2(j,j)=0; 
        h2(j,j+1)=-(1/2)*(-1/mh(j)/mh(j))*(-1/12*mh(j+2)+2/3*mh(j+1)-2/3*mh(j-1)+1/12*mh(j-2))/delr*(2/3/delr);  
        h2(j,j+2)=-(1/2)*(-1/mh(j)/mh(j))*(-1/12*mh(j+2)+2/3*mh(j+1)-2/3*mh(j-1)+1/12*mh(j-2))/delr*(-1/12/delr);
        h2(j,j-1)=-(1/2)*(-1/mh(j)/mh(j))*(-1/12*mh(j+2)+2/3*mh(j+1)-2/3*mh(j-1)+1/12*mh(j-2))/delr*(-2/3/delr);
        h2(j,j-2)=-(1/2)*(-1/mh(j)/mh(j))*(-1/12*mh(j+2)+2/3*mh(j+1)-2/3*mh(j-1)+1/12*mh(j-2))/delr*(1/12/delr);
    end
    
    % Sums the "normal" and variable mass correctd Hamiltonian matrices:
    q=h2(1:lengsim,1:lengsim)+h(1:lengsim,1:lengsim);
    
    % Sparse matrix eigenvale and eigenvector calculator.
    [V, D] = eigs(sparse(q), 15, 0);
    
    % MATLAB doesn't necessary present the eigenvalues and vectors in order. 
    % This part of the code addresses that:
    [d_notserp, idx] = sort(diag(D), 'ascend');
    Eens = diag(d_notserp);
    Vs = V(:, idx);
    
    %Lowest energy eigenvalue and eigenvector:
    EH1S=Eens(1,1);
    WF_H1SE1S=Vs(:,1)/Vs(1,1); %this is done to flip the WF over if it is negative
    
    %Normalization of the 1S hole state:
    norm1=0.0;
    for j=1:lengsim 
        norm1=norm1+(j*delr)^2*WF_H1SE1S(j)^2*delr;
    end
    WF_H1SE1S=WF_H1SE1S/sqrt(norm1);
    save WF_H1SE1S WF_H1SE1S;
end


function [EH1P]=CdSeCdS_H1SE1P(VXh,mh,UeP,lengsim,delr) 
    % This is a 1S hole with the electrostatic potentials arising form the 
    % electron 1P state.
    l=0;   
    delrsq=delr*delr;
    
    % Initialization of the Hamiltonian FDM matrix. h is the traditional form, 
    % h2 takes into account the corrections necessary due to the varialbe mass 
    % between core and shell.
    h=zeros([lengsim lengsim+2]);
    h2=zeros([lengsim lengsim+2]);
    
    % Implementing boundary conditions into FDM. Below sets the wavefunction
    % derivative equal 0 at the origin.
    r=delr;
    h(1,1)=(1/2/mh(1))*(2/r)*1/delr; 
    h(1,2)=0;
    h(1,3)=(1/2/mh(1))*-(2/r)*1/delr;
    
    % Given the need to implement a boundary condition, h2 must be 0.
    h2(1,2)=0; 
    h2(1,1)=0; 
    
    % The 3point stencil structure of the FDM requires a separate
    % initialization for the 2nd row compared to the remaining matrix.
    for j=2:2
        r=j*delr;
        h(j,j)=(-1/2/mh(j))*(-5/2/delr/delr)+UeP(j)+VXh(j)+1/2/mh(2)*l*(l+1)/r^2; 
        h(j,j+1)=(-1/2/mh(j))*(4/3/delrsq+(2/r)*2/3/delr);  
        h(j,j+2)=(-1/2/mh(j))*(-1/12/delrsq-(2/r)*1/12/delr); 
        h(j,j-1)=(-1/2/mh(j))*(4/3/delrsq-(2/r)*2/3/delr);
        h2(j,j)=0; 
        h2(j,j+1)=-(1/2)*(-1/mh(j)/mh(j))*((mh(j+1)-mh(j-1))/2/delr)*2/3/delr;  
        h2(j,j+2)=-(1/2)*(-1/mh(j)/mh(j))*((mh(j+1)-mh(j-1))/2/delr)*(-1/12/delr);
        h2(j,j-1)=-(1/2)*(-1/mh(j)/mh(j))*((mh(j+1)-mh(j-1))/2/delr)*(-2/3/delr);
    end
    % The remaining FDM matrix:
    for j=3:lengsim
        r=j*delr;
        h(j,j)=(-1/2/mh(j))*(-5/2/delrsq)+UeP(j)+VXh(j)+1/2/mh(j)*l*(l+1)/r^2; 
        h(j,j+1)=(-1/2/mh(j))*(4/3/delrsq+(2/r)*2/3/delr);  
        h(j,j+2)=(-1/2/mh(j))*(-1/12/delrsq-(2/r)*1/12/delr);
        h(j,j-1)=(-1/2/mh(j))*(4/3/delrsq-(2/r)*2/3/delr);
        h(j,j-2)=(-1/2/mh(j))*(-1/12/delrsq+(2/r)*1/12/delr);
        h2(j,j)=0; 
        h2(j,j+1)=-(1/2)*(-1/mh(j)/mh(j))*(-1/12*mh(j+2)+2/3*mh(j+1)-2/3*mh(j-1)+1/12*mh(j-2))/delr*(2/3/delr);  
        h2(j,j+2)=-(1/2)*(-1/mh(j)/mh(j))*(-1/12*mh(j+2)+2/3*mh(j+1)-2/3*mh(j-1)+1/12*mh(j-2))/delr*(-1/12/delr);
        h2(j,j-1)=-(1/2)*(-1/mh(j)/mh(j))*(-1/12*mh(j+2)+2/3*mh(j+1)-2/3*mh(j-1)+1/12*mh(j-2))/delr*(-2/3/delr);
        h2(j,j-2)=-(1/2)*(-1/mh(j)/mh(j))*(-1/12*mh(j+2)+2/3*mh(j+1)-2/3*mh(j-1)+1/12*mh(j-2))/delr*(1/12/delr);
    end
    
    % Sums the "normal" and variable mass correctd Hamiltonian matrices:
    q=h2(1:lengsim,1:lengsim)+h(1:lengsim,1:lengsim);
    
    % Sparse matrix eigenvale and eigenvector calculator.
    [V, D] = eigs(sparse(q), 15, 0);
    
    % MATLAB doesn't necessary present the eigenvalues and vectors in order. 
    % This part of the code addresses that:
    [d_notserp, idx] = sort(diag(D), 'ascend');
    Eens = diag(d_notserp);
    Vs = V(:, idx);
    
    %Lowest energy eigenvalue and eigenvector:
    EH1P=Eens(1,1);
    WF_H1SE1P=Vs(:,1)/Vs(1,1); %this is done to flip the WF over if it is negative
    
    %Normalization of the 1S state:
    norm1=0.0;
    for j=1:lengsim  
        norm1=norm1+(j*delr)^2*WF_H1SE1P(j)^2*delr;
    end
    WF_H1SE1P=WF_H1SE1P/sqrt(norm1);
    save WF_H1SE1P WF_H1SE1P;
end

function [EH2S]=CdSeCdS_H1SE2S(VXh,mh,Ue,lengsim,delr)
    % The following uses a 3 point stencil to solve the SE for a 2S hole state.
    %To make it clear this is a 2S state, angular momentum is 0.
    l=0; 
    delrsq=delr*delr;
    
    % Initialization of the Hamiltonian FDM matrix. h is the traditional form, 
    % h2 takes into account the corrections necessary due to the varialbe mass 
    % between core and shell.
    h=zeros([lengsim lengsim+2]);
    h2=zeros([lengsim lengsim+2]);
    
    % Implementing boundary conditions into FDM. Below sets the wavefunction
    % derivative equal 0 at the origin.
    r=delr; 
    h(1,1)=(1/2/mh(1))*(2/r)*1/delr; 
    h(1,2)=0;
    h(1,3)=(1/2/mh(1))*-(2/r)*1/delr;
    
    % Given the need to implement a boundary condition, h2 must be 0.
    h2(1,2)=0; 
    h2(1,1)=0; 
    
    % The 3point stencil structure of the FDM requires a separate
    % initialization for the 2nd row compared to the remaining matrix.
    for j=2:2
        r=j*delr;
        h(j,j)=(-1/2/mh(j))*(-5/2/delr/delr)+Ue(j,2)+VXh(j)+1/2/mh(2)*l*(l+1)/r^2; 
        h(j,j+1)=(-1/2/mh(j))*(4/3/delrsq+(2/r)*2/3/delr);  
        h(j,j+2)=(-1/2/mh(j))*(-1/12/delrsq-(2/r)*1/12/delr); 
        h(j,j-1)=(-1/2/mh(j))*(4/3/delrsq-(2/r)*2/3/delr);
        h2(j,j)=0; 
        h2(j,j+1)=-(1/2)*(-1/mh(j)/mh(j))*((mh(j+1)-mh(j-1))/2/delr)*2/3/delr;  
        h2(j,j+2)=-(1/2)*(-1/mh(j)/mh(j))*((mh(j+1)-mh(j-1))/2/delr)*(-1/12/delr);
        h2(j,j-1)=-(1/2)*(-1/mh(j)/mh(j))*((mh(j+1)-mh(j-1))/2/delr)*(-2/3/delr);
    end
    % The remaining FDM matrix:
    for j=3:lengsim
        r=j*delr;
        h(j,j)=(-1/2/mh(j))*(-5/2/delrsq)+Ue(j,2)+VXh(j)+1/2/mh(j)*l*(l+1)/r^2; 
        h(j,j+1)=(-1/2/mh(j))*(4/3/delrsq+(2/r)*2/3/delr);  
        h(j,j+2)=(-1/2/mh(j))*(-1/12/delrsq-(2/r)*1/12/delr);
        h(j,j-1)=(-1/2/mh(j))*(4/3/delrsq-(2/r)*2/3/delr);
        h(j,j-2)=(-1/2/mh(j))*(-1/12/delrsq+(2/r)*1/12/delr);
        h2(j,j)=0; 
        h2(j,j+1)=-(1/2)*(-1/mh(j)/mh(j))*(-1/12*mh(j+2)+2/3*mh(j+1)-2/3*mh(j-1)+1/12*mh(j-2))/delr*(2/3/delr);  
        h2(j,j+2)=-(1/2)*(-1/mh(j)/mh(j))*(-1/12*mh(j+2)+2/3*mh(j+1)-2/3*mh(j-1)+1/12*mh(j-2))/delr*(-1/12/delr);
        h2(j,j-1)=-(1/2)*(-1/mh(j)/mh(j))*(-1/12*mh(j+2)+2/3*mh(j+1)-2/3*mh(j-1)+1/12*mh(j-2))/delr*(-2/3/delr);
        h2(j,j-2)=-(1/2)*(-1/mh(j)/mh(j))*(-1/12*mh(j+2)+2/3*mh(j+1)-2/3*mh(j-1)+1/12*mh(j-2))/delr*(1/12/delr);
    end
    
    % Sums the "normal" and variable mass correctd Hamiltonian matrices:
    q=h2(1:lengsim,1:lengsim)+h(1:lengsim,1:lengsim);
    
    % Sparse matrix eigenvale and eigenvector calculator.
    [V, D] = eigs(sparse(q), 15, 0);
    
    % MATLAB doesn't necessary present the eigenvalues and vectors in order. 
    % This part of the code addresses that:
    [d_notserp, idx] = sort(diag(D), 'ascend');
    Eens = diag(d_notserp);
    Vs = V(:, idx);
    
    % 1st excited state eigenvalue and eigenvector:
    EH2S=Eens(1,1);
    WF_H1SE2S=Vs(:,1)/Vs(1,1); %this is done to flip the WF over if it is negative
    
    %Normalization of the 2S state:
    norm1=0.0;
    for j=1:lengsim 
        norm1=norm1+(j*delr)^2*WF_H1SE2S(j)^2*delr;
    end
    WF_H1SE2S=WF_H1SE2S/sqrt(norm1);
    save WF_H1SE2S WF_H1SE2S;
end



function [blank]=Poisson_e(eps,lengsim,delr) 
%This is the electric field from the electron 1S and 2S states

delrsq=delr*delr;
load WF_E1SH1S; load WF_E2SH1S; 

%Traditional expression for the charge density:
for j=1:lengsim
    rho1(j)=1/4/pi*WF_E1SH1S(j)^2;
    rho2(j)=1/4/pi*WF_E2SH1S(j)^2;
end
% Some housekeeping on vector directions for the solver:
rho1=rho1';rho2=rho2'; 

% Matrix initializations:
final=sparse(lengsim,lengsim);
h=zeros([lengsim lengsim+2]); h2=zeros([lengsim lengsim+2]);

% Buiding the FDM matrix. h2 accounts for the variable dielectric
for j=1:1
    r=delr;
    h(j,j)=eps(j)*(-2/delrsq); 
    h(j,j+1)=eps(j)*(1/delrsq+(2/r)*1/2/delr);  
    h2(j,j)=0;
    h2(j,j+1)=((eps(j+1)-eps(j))/1/delr)*1/2/delr;
end
for j=2:lengsim
    r=j*delr;
    h(j,j)=eps(j)*(-2/delrsq); 
    h(j,j+1)=eps(j)*(1/delrsq+(2/r)*1/2/delr);  
    h(j,j-1)=eps(j)*(1/delrsq-(2/r)*1/2/delr);
    h2(j,j)=0; 
    h2(j,j+1)=((eps(j+1)-eps(j-1))/2/delr)*1/2/delr;  
    h2(j,j-1)=((eps(j+1)-eps(j-1))/2/delr)*-1/2/delr;
end

% Composing the total FDM matrix with the correct dimensionality
final=h2(1:lengsim,1:lengsim)+h(1:lengsim,1:lengsim);

% The following accounts for the fact that the solver only gets the shape
% of the potential correct. As such, it needs to be scaled to match a
% Coulomb potential at long range. The first part creates the correct ~1/r
% potential, which will be used to scale Ue.
for i=1:lengsim 
    rax=i*delr; 
    sig(i)=-1/eps(lengsim)/rax;
end

% The solutin to the Poisson equation:
tempUe=final\rho1;

% Scaling the potentials to the Coulomb field far from the QD center
Ue(:,1)=tempUe*(sig(lengsim-1)-sig(lengsim))/(tempUe(lengsim-1)-tempUe(lengsim));
Ue(:,1)=Ue(:,1)-(Ue(lengsim,1)-sig(lengsim));
tempUe=final\rho2;
Ue(:,2)=tempUe*(sig(lengsim-1)-sig(lengsim))/(tempUe(lengsim-1)-tempUe(lengsim));
Ue(:,2)=Ue(:,2)-(Ue(lengsim,2)-sig(lengsim));
save Ue Ue;
end


function [blank]=Poisson_eP(eps,lengsim,delr) 
%This is to solve the fields from the 1P electron.

delrsq=delr*delr; load WF_E1PH1S; 

% The carge density. This is a non-standard form, which is used due to the
% difficulties in calculating the quadrupole field. The quadrupole is
% ultimately scaled so this isn't important
for j=1:lengsim
    rho(j)=WF_E1PH1S(j)^2; 
end
% Some housekeeping on vector directions for the solver:
rho=rho'; 

% Matrix initializations:
final=sparse(lengsim,lengsim);finalq=sparse(lengsim,lengsim);
h=zeros([lengsim lengsim+2]); hq=zeros([lengsim lengsim+2]);
h2=zeros([lengsim lengsim+2]);

% Buiding the FDM matrix. h2 accounts for the variable dielectric
for j=1:1
    r=j*delr; rp=(j+1)*delr; rm=(j-1)*delr;
    h(j,j)=eps(j)*(-2/delrsq); 
    h(j,j+1)=eps(j)*(1/delrsq+(2/r)*1/2/delr);  
    h2(j,j)=0;
    h2(j,j+1)=((eps(j+1)-eps(j))/1/delr)*1/2/delr;
    hq(j,j)=eps(j)*((-2/delrsq)-6/r^2);
    hq(j,j+1)=eps(j)*(1/delrsq+(2/rp)*1/2/delr);
end
for j=2:lengsim
    r=j*delr; rp=(j+1)*delr; rm=(j-1)*delr; %This is necessary for hq to work
    h(j,j)=eps(j)*(-2/delrsq); 
    h(j,j+1)=eps(j)*(1/delrsq+(2/r)*1/2/delr);  
    h(j,j-1)=eps(j)*(1/delrsq-(2/r)*1/2/delr);
    h2(j,j)=0; 
    h2(j,j+1)=((eps(j+1)-eps(j-1))/2/delr)*1/2/delr;  
    h2(j,j-1)=((eps(j+1)-eps(j-1))/2/delr)*-1/2/delr;
    hq(j,j)=eps(j)*((-2/delrsq)-6/r^2); 
    hq(j,j+1)=eps(j)*(1/delrsq+(2/rp)*1/2/delr);  
    hq(j,j-1)=eps(j)*(1/delrsq-(2/rm)*1/2/delr);
end

% Composing the total FDM matrix with the correct dimensionality
final=h2(1:lengsim,1:lengsim)+h(1:lengsim,1:lengsim);
finalq=h2(1:lengsim,1:lengsim)+hq(1:lengsim,1:lengsim);

% The following accounts for the fact that the solver only gets the shape
% of the potential correct. As such, it needs to be scaled to match a
% Coulomb or Quadrupole potential at long range. The first part creates 
% the correct ~1/r^3 and ~1/r potentials, which will be used to scale UeP 
% and UePQ.
Qzz=0;
for i=1:lengsim
    Qzz=Qzz+1/2*rho(i)*(i*delr)^4*delr; 
end 
for j=1:lengsim
    sig1(j)=-1/eps(lengsim)/(j*delr);
    sig2(j)=2*-1/2*Qzz/eps(lengsim)/(j*delr)^3; %Mult by 2 b/c of 3cos^2-1=2 zz 
end

% The solutin to the Poisson equation for the electric potential from a
% p-symmetry source:
tempUeP=final\rho;
% Scaling the potentials to the Coulomb field far from the QD center
UeP=tempUeP*(sig1(lengsim-1)-sig1(lengsim))/(tempUeP(lengsim-1)-tempUeP(lengsim));
UeP=UeP-(UeP(lengsim)-sig1(lengsim));

% The solutin to the Poisson equation for the quadrupolar field from a 
% p-symmetry source:
tempUePQ=finalq\rho;
% Scaling the potentials to the Coulomb field far from the QD center
UePQ=tempUePQ*(sig2(lengsim-1)-sig2(lengsim))/(tempUePQ(lengsim-1)-tempUePQ(lengsim));
UePQ=UePQ-(UePQ(lengsim)-sig2(lengsim));

save UeP UeP;
save UePQ UePQ;
end

function [blank]=Poisson_h(eps,lengsim,delr) 
%This is the electric field from the hole's 1S state fields 

delrsq=delr*delr;
load WF_H1SE1S; load WF_H1SE1P.mat; load WF_H1SE2S.mat; 
%Traditional expression for the charge density:

for j=1:lengsim
    rax=j*delr; 
    rho1(j)=-1/4/pi*WF_H1SE1S(j)^2; 
    rho2(j)=-1/4/pi*WF_H1SE1P(j)^2; 
    rho3(j)=-1/4/pi*WF_H1SE2S(j)^2; 
end
% Some housekeeping on vector directions for the solver:
rho1=rho1';rho2=rho2';rho3=rho3';

% Matrix initializations:
final=sparse(lengsim,lengsim);
h=zeros([lengsim lengsim+2]);
h2=zeros([lengsim lengsim+2]);

% Buiding the FDM matrix. h2 accounts for the variable dielectric
for j=1:1
    r=delr;
    h(j,j)=eps(j)*(-2/delrsq); 
    h(j,j+1)=eps(j)*(1/delrsq+(2/r)*1/2/delr);  
    h2(j,j)=0;
    h2(j,j+1)=((eps(j+1)-eps(j))/1/delr)*1/2/delr;
end
for j=2:lengsim
    r=j*delr;
    h(j,j)=eps(j)*(-2/delrsq); 
    h(j,j+1)=eps(j)*(1/delrsq+(2/r)*1/2/delr);  
    h(j,j-1)=eps(j)*(1/delrsq-(2/r)*1/2/delr);
    h2(j,j)=0; 
    h2(j,j+1)=((eps(j+1)-eps(j-1))/2/delr)*1/2/delr;  
    h2(j,j-1)=((eps(j+1)-eps(j-1))/2/delr)*-1/2/delr;
end

% Composing the total FDM matrix with the correct dimensionality
final=h2(1:lengsim,1:lengsim)+h(1:lengsim,1:lengsim);

for i=1:lengsim 
    rax=i*delr; 
    sig(i)=1/eps(lengsim)/rax;
end


% The solutin to the Poisson equation:
tempUh=final\rho1;
% Scaling the potentials to the Coulomb field far from the QD center
Uh(:,1)=tempUh*(sig(lengsim-1)-sig(lengsim))/(tempUh(lengsim-1)-tempUh(lengsim));
Uh(:,1)=Uh(:,1)-(Uh(lengsim,1)-sig(lengsim));
% The solutin to the Poisson equation:
tempUh=final\rho2;
% Scaling the potentials to the Coulomb field far from the QD center
Uh(:,2)=tempUh*(sig(lengsim-1)-sig(lengsim))/(tempUh(lengsim-1)-tempUh(lengsim));
Uh(:,2)=Uh(:,2)-(Uh(lengsim,2)-sig(lengsim));
% The solutin to the Poisson equation:
tempUh=final\rho3;
% Scaling the potentials to the Coulomb field far from the QD center
Uh(:,3)=tempUh*(sig(lengsim-1)-sig(lengsim))/(tempUh(lengsim-1)-tempUh(lengsim));
Uh(:,3)=Uh(:,3)-(Uh(lengsim,3)-sig(lengsim));

save Uh Uh;

end
