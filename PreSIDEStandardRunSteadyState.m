%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                                                             %%%%%%%
%%%%%%%       Precambrian Silica Induced DolomitE model             %%%%%%%
%%%%%%%                         PreSIDE                             %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                                                             %%%%%%%
%%%%%%%       Yihang Fang, University of Wisconsin-Madison          %%%%%%%
%%%%%%%                      March - April 2022                     %%%%%%%
%%%%%%%                                                             %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all all;

%Simulation start and end time
tStart = 0;
tFinal = 30e6;
%Time step. Penman and Rooney (2019) suggest too coarse may cause CO2 to go negative if air-sea gas exchange is large, 
%suggest <1.
tStep = .05; 

%Boolean for controlling global glaciation
snowball=0;

%round up so that duration is an even multiple of tStep
if(mod(tFinal-tStart,tStep)>0)
    tFinal = tFinal + (tStep-mod((tFinal-tStart),tStep)); 
end

%Write final tracer concentrations to file SVSTRT.txt if true
SVSTRT=0;
%Open saved tracer concentrations from SVSTRT.txt if true
RSTRT=0;

%CONSTANTS THAT PROBABLY SHOULD NOT CHANGE
%Reservoir volumes
vsurf = 3.58e16; %m^3
vdeep = 1.3e18; %m^3
%Area of ocean, m^2
aoc = vsurf/100;
%Vertical mixing coefficient, m/yr, Yool & Tyrell
v = 3;
%Equilibrium pCO2 in ppm
pCO20 = 1000.;

%Temperature, salinity used in carbonate chemistry eq. constants
%(degrees C)
surfT0 = 10;
surfT = surfT0;
deepT0 = 2.;
deepT = deepT0;
%Temperature relaxation timescales, years
Trelaxsurf=20.;
Trelaxdeep=1000.;

surfSal = 35;
deepSal = 35;
CaC = 10.3e-3;% [Ca2+] concentration in mol/kg, modern
MgC = 50.3e-3;
pH = 8.3;
aH = 10^(-pH);
surfpH = vsurf*aH;

%Silicate & carbonate weathering strength exponents
nSi = 0.5;
nc = 0.5; 
fCarbw0 = 12.e12; %following LOSCAR
%air-sea gas exchange constant in units mol/uatm/m2/yr
kas = 0.06; % from LOSCAR
ppmtomolC = 2.2e15/12.;% from LOSCAR: 1ppm atm CO2 = 2.2GtC in atmosphere
m3tokg=1.025e3; %to convert from mol/kg to mol/m3
Tsens=3.; %temperature sensitivity to doubling of pCO2 (relative pCO20)

% initial (steady state) tracer concentrations (all in mol/m^3) 
% obtained by running model for 10e7 years w/o perturbation
% these will change as the model runs
surfSi = 1.250746; %mol/m^3, modern = 0.0077
deepSi = 1.251305; %mol/m^3, modern = .1103
surfDIC = 6.448753;
deepDIC = 6.451919;
surfTA = 6.954368;
deepTA = 6.954368;
atmCO2 = pCO20*ppmtomolC; %atmospheric CO2, in total mol C

if(RSTRT)
    disp('using restart values from SVSTRT.txt')
    fileID = fopen('SVSTRT.txt');
    RSTRTvalues = fscanf(fileID,'%f');
    surfSi = RSTRTvalues(2);
    deepSi = RSTRTvalues(3);
    surfDIC = RSTRTvalues(5);
    deepDIC = RSTRTvalues(6);
    surfTA = RSTRTvalues(7);
    deepTA = RSTRTvalues(8);
    atmCO2 = RSTRTvalues(11)*ppmtomolC;
    CaC = RSTRTvalues(12);
    MgC=RSTRTvalue(13);
end

%Initial fluxes, molesyear
fdust = 0.5e12; %dust input
fHydro = 0.6e12; %hydrothermal Si input
fSilw = 5.e12;
fCarbw = 12.e12;
fVolc = 5.e12; %volcanic CO2 degassing

%declare variables used in for loop, these values will be calculated each
%time step
%f_ for individual fluxes, in mol/yr
fSiBurial = 0; %SiO2 burial
fCaCO3Burial = 0;
fCaMgCO3Burial = 0;
fCorg = 0;
fvSDSi = 0; %Surface-Deep vertical Si mixing
fvDSSi = 0; %Deep-Surface vertical Si mixing
fvSDDIC = 0;
fvDSDIC= 0;
fvSDTA = 0;
fvDSTA = 0;
pCO2atm = 280.;
pCO2surf = 0;
CO2aq = 0;
CO3 = 0;
%d_ for integrated fluxes (dX/dt), one per tracer per reservoir, in mol/yr
dsurfSi = 0; 
ddeepSi = 0;
dsurfDIC = 0; 
ddeepDIC = 0;
dsurfTA = 0;
ddeepTA = 0;
datmCO2 = 0;
CO32 = 0;
HCO3 = 0;
%Carb chem constants
kh = 0;
k1 = 0;
k2 = 0;
ksp = 0;
kw = 0;
kb = 0;
bor = 0;
h = 1e-8;

%this sets how often results are recorded - don't want millions of rows
%e.g. if recordStep = 10, then every 10th time step is recorded
recordStep=1;
while((tFinal-tStart)/tStep/recordStep>10000.)
    recordStep=recordStep*10;
end

%Initialize output matrix:
%column 1 = time (years)
%column 2 = surfSi (mol/m^3)
%column 3 = deepSi (mol/m^3)
%column 4 = Si burial rate (mol/yr)
%column 5 = surf DIC (mol/m^3)
%column 6 = deep DIC (mol/m^3)
%column 7 = surf TA (mol/m^3)
%column 8 = deep TA (mol/m^3)
%column 9 = surface omega (calcite saturation state)
%column 10 = carbonate burial rate (mol/yr)
%column 11 = pCO2 (ppmv)
results = zeros (((tFinal-tStart)/tStep/recordStep),11);

tic;

presnowball=1;

%time loop
for t = tStart:tStep:tFinal 
     
    %record time each recordStep interval
    if(mod(t,tStep*recordStep))==0
        %otherwise results matrix indices may not be integers!
        row=int32(t/tStep/recordStep); 
        results(row+1,1) = t;
    end
    
    %update progress on each 10%
    if(mod(t-tStart,(tFinal-tStart)/10)<tStep)
        t
        pCO2atm
        toc
    end
    
    %convert atmCO2 (mol) to pCO2 (ppmv)
    pCO2atm = atmCO2/ppmtomolC; 
    
    %Snowball after 1myr
    if(t>15.e6 && presnowball)
        snowball=1;
        presnowball=0;
    end
    
    %check for deglaciation!
    if(snowball)
        if(pCO2atm>1.2e5)%threshold 0.12 bar
            snowball=0;
        end
    end
    
    %Calculate silicate weathering flux
    fSilw = fVolc * (pCO2atm/pCO20)^nSi;
    if(snowball)
        fSilw=0.;
    end
    
    %Calculate carbonate weathering flux
    fCarbw = fCarbw0 * (pCO2atm/pCO20)^nc;
    if(snowball)
        fCarbw=0.;
    end
    
    %calculate T's as a function of pCO2
    %with different relaxation timescales for surface and deep
    deltaT=log2(pCO2atm/pCO20)*Tsens;
    surfT = surfT + (surfT0 + deltaT - surfT)/Trelaxsurf;
    deepT = deepT + (surfT0 + deltaT - surfT)/Trelaxsurf;
    if(snowball)
        %Temperatures relax to 0 degrees
        surfT = surfT + (0. - surfT)/Trelaxsurf;
        deepT = deepT + (0. - deepT)/Trelaxdeep;
    end
    
    %No aeolian Si input during snowball
    if(snowball)
        fdust=0.;
    else
        fdust = 0.5e12;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Carbonate Chemistry Calculations %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %temperature in Kelvin
    tk=surfT+273.15;
    
    %kh = Henry's law constant of CO2 solubility in seawater (~[CO2]/pCO2)
    %from Weiss (1974), following LOSCAR
    kh=exp((9345.17/tk - 60.2409 + 23.3585*log(tk/100.))+surfSal*(0.023517-0.00023656*tk+0.0047036e-4*tk*tk));
    
    %k1 = first dissociation constant of carbonic acid (~=HCO3/H2CO3)
    %Mehrbach et al. (1973) refit by Lueker et al. (2000)
    k1=10^(-((3633.86/tk - 61.2172 + 9.6777*log(tk) - 0.011555*surfSal)+0.0001152*surfSal*surfSal));
    
    %k2 = second dissociation constant of carbonic acid (~=CO3/HCO3)
    %Mehrbach et al. (1973) refit by Lueker et al. (2000)
    k2=10^(-((471.78/tk + 25.9290 - 3.16967*log(tk) - 0.01781*surfSal)+0.0001122*surfSal*surfSal));
    
    %kb = dissociation constant of boric acid
    %(Dickson, 1990 in Dickson and Goyet, 1994, Chapter 5, p. 14)
    tmp1  = -8966.90-2890.53*sqrt(surfSal)-77.942*surfSal;
    tmp1 = tmp1 + 1.728*surfSal^1.5-0.0996*surfSal*surfSal;
    tmp2  =  148.0248+137.1942*sqrt(surfSal)+1.62142*surfSal;
    tmp3  = -24.4344-25.085*sqrt(surfSal)-0.2474*surfSal;
    tmp3 = tmp3 * log(tk);
    tmp4  = tmp1/tk + tmp2 + tmp3 + 0.053105*sqrt(surfSal)*tk; 
    kb = exp(tmp4);
    bor=4.16e-4*surfSal/35;
    
    %kw = ion product of water
    %Millero (1995)(in Dickson and Goyet (1994, Chapter 5, p.18))
    tmp1  = -13847.26/tk + 148.96502 - 23.6521*log(tk);
    tmp2  =  118.67/tk - 5.977 + 1.0495*log(tk);
    tmp2 = tmp2 * sqrt(surfSal);
    tmp1 = tmp1 + tmp2 - 0.01615*surfSal;	
    kw  = exp(tmp1); 
    
    %ksp = solubility product of calcite (~=[Ca][CO3])
    %Mucci 1983, units (mol/kg-soln)^2
    ksp=10^((-171.9065-0.077993*tk+2839.319/tk+71.595*log10(tk))+((-0.77712+0.0028426*tk+178.34/tk)*sqrt(surfSal))+(-0.07711*surfSal+0.0041249*surfSal^1.5));
    %Benezeth 2018, dolomite solubility product
    kspdol=10^((17.502-0.0689*tk-4220.119/tk)+((-0.77712+0.0028426*tk+178.34/tk)*sqrt(surfSal))+(-0.07711*surfSal+0.0041249*surfSal^1.5));
    %TA and DIC need to be in mol/kg for carb chem calculations
    surfTA=surfTA/m3tokg; 
    surfDIC=surfDIC/m3tokg;
    
    %Calculations for full carbonate chemistry from TA and DIC from csys.m
    p5  = -1.;
    p4  = -surfTA-kb-k1;
    p3  = surfDIC*k1-surfTA*(kb+k1)+kb*bor+kw-kb*k1-k1*k2;
    tmp = surfDIC*(kb*k1+2.*k1*k2)-surfTA*(kb*k1+k1*k2)+kb*bor*k1;
    p2  = tmp+(kw*kb+kw*k1-kb*k1*k2);
    tmp = 2.*surfDIC*kb*k1*k2-surfTA*kb*k1*k2+kb*bor*k1*k2;
    p1  = tmp+(+kw*kb*k1+kw*k1*k2);
    p0  = kw*kb*k1*k2;
    p   = [p5 p4 p3 p2 p1 p0];
    r   = roots(p);
    h   = max(real(r));
    s = surfDIC / (1.+k1/h+k1*k2/h/h);
    
    HCO3 = surfDIC/(1+h/k1+k2/h);
    
    CO3 = surfDIC/(1+h/k2+h*h/k1/k2);
    
    fCO2 = s/kh;
    R    = 83.14510;        % mol bar deg-1 
    Pstd = 1.01325;
    delC = (57.7 - 0.118.*tk);
    B = -1636.75 + 12.0408.*tk - 0.0327957.*tk.^2 + 3.16528.*0.00001.*tk.^3;
    p2f = exp((B + 2.*delC).*Pstd./(R*tk));
    
    pCO2surf = fCO2/p2f*1.e6;
    
    %translate TA and DIC back to mol/m3
    surfTA=surfTA*m3tokg;
    surfDIC=surfDIC*m3tokg;
    
    %Calculate Omega as function of Ksp
    OmegaClc = (CaC*CO3)/ksp; %CO3 needs units of mol/kg
    OmegaDol= (CaC*MgC*CO3)/kspdol;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Enough Carbonate Chemistry For Now  %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Inorganic carbonate precipitation following Zeebe's "Strangelove Ocean"
    %With critical saturation for precipitation at omega=10. 
    %Very fast above that. 
    e=2.71828182846;
    fCaCO3Burial=(fCarbw0)*(OmegaClc-10.)^3;
    fCaMgCO3Burial=(fCarbw0+fSilw)*(e^((surfSi-1))-1)*(OmegaDol-10.)^3;
    %fCaMgCO3Burial=((OmegaDol*(e^((surfSi-1))-1))-10.)^3;
    %%fCaMgCO3Burial=(fVolc+fCarbw0)*(OmegaDol-10.)^3;
    if(fCaCO3Burial<0)
        fCaCO3Burial=0;
    %   fCaMgCO3Burial=0;
    elseif(fCaMgCO3Burial<0)
        fCaMgCO3Burial=0;
    %elseif(surfSi<1.5)
        fCaMgCO3Burial=0;
    %elseif(pH<7)
    %    fCaMgCO3Burial=0;
    %    fCaCO3Burial=0;
    end
    
    %Corg export some constant (a fraction of equilibrium carbonate flux)
    fCorg = 0.2*(fVolc+fCarbw0);
    %But I figure it'd shut down under ice
    if(snowball)
        fCorg = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% BEGIN atmCO2 FLUXES %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %air-sea CO2 gas exchange following LOSCAR
    fASCO2=kas*aoc*(pCO2atm-pCO2surf);
    if(snowball)
        %hard or soft snowball?
        fASCO2=0;
    end
    
    %flux balance for atmCO2
    datmCO2 = fVolc - fASCO2 - 2*fSilw - fCarbw;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% END atmCO2 FLUXES %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% BEGIN Dissolevd Si FLUXES %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % opal burial is 0 if < saturation (1200mol/m3), then linear above that
    % balances Si inputs at 2000mol/m3
    fSiBurial = 6.7e12/0.55*(surfSi-.75);
    if(fSiBurial<0)
        fSiBurial=0;
    end
    
    fvSDSi = surfSi*v*aoc; %surface to deep flux
    fvDSSi = deepSi*v*aoc; %deep to surface flux
    
    %flux balance for surface Si (mol/yr)
    % =riverine input + aeolian input - diatom uptake + surface dissolution
    % + mixing with deep
    dsurfSi = fSilw + fdust + fvDSSi - fvSDSi - fSiBurial;
   
    
    
 
    
    %flux balance for deep (mol/yr)
    % = hydrothermal input + deep/sed dissolution 
    ddeepSi = fHydro + fvSDSi - fvDSSi;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% END Dissolevd Si FLUXES %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %oceanic CaCO3 dissolution during snowball
    snowballCaCO3diss=0;
    snowballCaMgCO3diss=0;
    if(snowball)
        %Obviously this will only occur while Omega<1)
        if(OmegaClc<1)
            %It should be scaled to ss terrestrial CaCO3 weathering?
            snowballCaCO3diss=(1-OmegaClc)*fCarbw0;
        elseif(OmegaDol<1)
            snowballCaMgCO3diss=(1-OmegaDol)*fCarbw0;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% BEGIN DIC FLUXES %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fvSDDIC = surfDIC*v*aoc; %surface to deep flux
    fvDSDIC = deepDIC*v*aoc; %deep to surface flux

    %Flux balances for surface and deep DIC
    dsurfDIC = 2*fSilw + 2*fCarbw +snowballCaCO3diss+ 2*snowballCaMgCO3diss+ fASCO2 - fCaCO3Burial - fCorg + fvDSDIC - fvSDDIC;
    ddeepDIC = fvSDDIC + fCorg - fvDSDIC;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% END DIC FLUXES %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% BEGIN pH FLUXES %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   % CO32 = (surfDIC*(10^-9.13)/(surfpH))*((1+((10^-9.13)/surfpH))^-1);
   % HCO3 = surfDIC-CO32;
   % dsurfpH = dsurfDIC+fCaCO3Burial+2*fCaMgCO3Burial-CO32;
   %if (pH<6)
   %    dsurfpH=-CO32;
   % end
   
   pH = 9.3-log10(surfTA);
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% END pH FLUXES %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% BEGIN TA FLUXES %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fvSDTA = surfTA*v*aoc; %surface to deep flux
    fvDSTA = deepTA*v*aoc; %deep to surface flux
    
    dsurfTA = 2*fSilw + 2*fCarbw + 2*snowballCaCO3diss +2*snowballCaMgCO3diss- 2*fCaCO3Burial + fvDSTA - fvSDTA;
    ddeepTA = fvSDTA - fvDSTA;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% END TA FLUXES %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%% Ca2+ FLUXES %%%%%%%%%%%
    dCaC = 0.5*fSilw + 0.5*fCarbw + 0.5*snowballCaCO3diss+0.5*snowballCaMgCO3diss- fCaCO3Burial;
    dMgC = 0.5*fSilw + 0.5*fCarbw + 0.5*snowballCaMgCO3diss-8*fHydro;
  
    if(snowball)
        if (snowball==0)
       dCaC = fSilw + fCarbw + snowballCaCO3diss+snowballCaMgCO3diss- fCaCO3Burial-fCaMgCO3Burial;
       dMgC = fSilw + fCarbw + snowballCaMgCO3diss - fCaMgCO3Burial-6*fHydro;
       dsurfTA = 2*fSilw + 2*fCarbw + 2*snowballCaCO3diss +2*snowballCaMgCO3diss-2*fCaMgCO3Burial- 2*fCaCO3Burial + fvDSTA - fvSDTA;
       dsurfDIC = 2*fSilw + 2*fCarbw +snowballCaCO3diss+ 2*snowballCaMgCO3diss+ fASCO2 - fCaCO3Burial - 2*fCaMgCO3Burial-fCorg + fvDSDIC - fvSDDIC;
       dsurfSi = fSilw + fdust + fvDSSi - fvSDSi - fSiBurial-(2*10^-6)*fCaMgCO3Burial;
 
        end
    end
    
    %%%%%%%% END Ca2+ FLUXES %%%%%%%%%
    
    
    
    %now integrate the fluxes (mol/m^3)
    surfSi = surfSi + dsurfSi*tStep/vsurf;
    deepSi = deepSi + ddeepSi*tStep/vdeep;
    surfDIC = surfDIC + dsurfDIC*tStep/vsurf;
    deepDIC = deepDIC + ddeepDIC*tStep/vdeep;
    surfTA = surfTA + dsurfTA*tStep/vsurf;
    deepTA = deepTA + ddeepTA*tStep/vdeep;
    atmCO2 = atmCO2 + datmCO2*tStep; %mol C
  

    %CaC is in mol/kg for some stupid reason
    CaC = CaC + dCaC*tStep/(vdeep+vsurf)/m3tokg;
    MgC = MgC + dMgC*tStep/(vdeep+vsurf)/m3tokg;
    if(mod(t,tStep*recordStep))==0
        %record results
        %column 1 = time (years)
        %column 2 = surfSi (mol/m^3)
        %column 3 = deepSi (mol/m^3)
        %column 4 = Si burial rate (mol/yr)
        %column 5 = surf DIC (mol/m^3)
        %column 6 = deep DIC (mol/m^3)
        %column 7 = surf TA (mol/m^3)
        %column 8 = deep TA (mol/m^3)
        %column 9 = surface omega (calcite saturation state)
        %column 10 = carbonate burial rate (mol/yr)
        %column 11 = pCO2 (ppmv)
        %column 12, [Ca2+], **mol/kg**
        results(row+1,2) = surfSi;
        results(row+1,3) = deepSi;
        results(row+1,4) = fSiBurial;
        results(row+1,5) = surfDIC;
        results(row+1,6) = deepDIC;
        results(row+1,7) = surfTA;
        results(row+1,8) = deepTA;
        results(row+1,9) = OmegaClc;
        results(row+1,10) = fCaCO3Burial;
        results(row+1,11) = pCO2atm;
        results(row+1,12) = CaC;
        results(row+1,13) = MgC;
        results(row+1,14) = OmegaDol;
        results(row+1,15) = fCaMgCO3Burial;
        results(row+1,16) = CO3;
        results(row+1,17) = kspdol;
        results(row+1,18) = ksp;
        results(row+1,19) = surfSal;
        results(row+1,20) = pH;
        results(row+1,21) = fSilw;
        results(row+1,22) = fCarbw;
    end
    
    %break loop, don't svstrt, print error if tracers or burial is negative
    if (surfSi < 0)
        fprintf('surface Si negative!\n')
        SVSTRT=0;
        break
    end
    
    if (deepSi < 0)
        fprintf('deep Si negative!\n')
        SVSTRT=0;
        break
    end
    if (fSiBurial < 0)
        fprintf('Si burial negative!\n')
        SVSTRT=0;
        break
    end
    if (pCO2atm < 0)
        fprintf('atm pCO2 negative!\n')
        SVSTRT=0;
        break
    end
    if (pCO2surf < 0)
        fprintf('surface pCO2 negative!\n')
        SVSTRT=0;
        break
    end
    
    
end

if(SVSTRT)
    disp('saving endstate values to snowball.txt')
    fileID = fopen('snowball.txt','w');
    fprintf(fileID,'%f %f  %f %f  %f %f  %f %f  %f %f  %f %f', t, surfSi, deepSi, fSiBurial, surfDIC, deepDIC, surfTA, deepTA, OmegaClc, fCaCO3Burial, pCO2atm, CaC);
    fclose(fileID);
end

%clean up results array
results=results(1:row+1,:);

%calculate excess SiO2 burial. NOTE: because not every tStep is recorded, 
%this calculation comes with rounding errors (<1%)
excessSiO2burial=zeros(1,length(results));
excessSiO2burial(1)=((results(1,4)-(fVolc+fdust+fHydro))*(results(2,1)-results(1,1)));
for i=2:1:(length(results))
    %results+=delta_t*(siO2burial-sio2burial_0)
    excessSiO2burial(i)=excessSiO2burial(i-1)+((results(i,4)-(fVolc+fdust+fHydro))*(results(i,1)-results(i-1,1)));
end


%calculate excess CaCO3 burial NOTE: because not every tStep is recorded, 
%this calculation comes with rounding errors (<1%)
excessCaCO3burial=zeros(1,length(results));
excessCaCO3burial(1)=((results(1,10)-(fVolc+fCarbw0))*(results(2,1)-results(1,1)));
for i=2:1:(length(results))
    %results+=delta_t*(CaCO3burial-CaCO3burial_0)
    excessCaCO3burial(i)=excessCaCO3burial(i-1)+((results(i,10)-(fVolc+fCarbw0))*(results(i,1)-results(i-1,1)));
end

excessCaMgCO3burial=zeros(1,length(results));
excessCaMgCO3burial(1)=((results(1,15)-(fVolc+fCarbw0))*(results(2,1)-results(1,1)));
for i=2:1:(length(results))
    %results+=delta_t*(CaCO3burial-CaCO3burial_0)
    excessCaMgCO3burial(i)=excessCaMgCO3burial(i-1)+((results(i,15)-(fVolc+fCarbw0))*(results(i,1)-results(i-1,1)));
end


%downsample for plotting efficiency
while (length(results)>100000)
    results = results(1:10:end,:);
end

%save results matrix
save('results.mat','results')

%plot results
%surface and deep silicic acid concentrations
subplot(8,1,1)
plot(results(:,1),results(:,2),results(:,1),results(:,3))
%set(gca,'YLim',[0,0.5]);
ylabel('[Si(OH)_4] (mol/m^3)','FontSize',15)
set(gca,'FontSize',10)
legend('surface','deep')

%omega
subplot(8,1,2)
plot(results(:,1),results(:,9),'Color','k')
%set(gca,'YLim',[6.e12,14.e12]);
ylabel('Omega')
xlabel('Time (years)')
set(gca,'FontSize',10)

%DIC
subplot(8,1,3)
plot(results(:,1),results(:,5),results(:,1),results(:,6))
ylabel('DIC (mol/m^3)')
xlabel('Time (years)')
set(gca,'FontSize',10)

%TA
subplot(8,1,4)
plot(results(:,1),results(:,7),results(:,1),results(:,8))
ylabel('TA (mol/m^3)')
xlabel('Time (years)')
set(gca,'FontSize',10)

%pCO2
subplot(8,1,5)
plot(results(:,1),results(:,11),'Color','k')
ylabel('pCO2 (ppmv)')
xlabel('Time (years)')
set(gca,'FontSize',10)

%pCO2
subplot(8,1,6)
plot(results(:,1),results(:,12),'Color','k')
ylabel('[Ca2+] (mol/kg)')
xlabel('Time (years)')
set(gca,'FontSize',10)

subplot(8,1,7)
plot(results(:,1),results(:,13),'Color','k')
ylabel('[Mg2+] (mol/kg)')
xlabel('Time (years)')
set(gca,'FontSize',10)

subplot(8,1,8)
plot(results(:,1),results(:,14),'Color','k')
ylabel('[Omega Dol] (mol/kg)')
xlabel('Time (years)')
set(gca,'FontSize',10)

figure;
%CaCO3 burial rate
subplot(3,1,1)
plot(results(:,1),results(:,10))
ylabel('CaCO3 burial (mol/yr)')
xlabel('Time (years)')
set(gca,'FontSize',10)

%SiO2 burial rate
subplot(3,1,2)
plot(results(:,1),results(:,4),'Color','k')
ylabel('SiO2 burial (mol/yr)')
xlabel('Time (years)')
set(gca,'FontSize',10)

%SiO2 burial rate
subplot(3,1,3)
plot(results(:,1),results(:,15),'Color','k')
ylabel('CaMgCO32 burial (mol/yr)')
xlabel('Time (years)')
set(gca,'FontSize',10)

%excess burial plots
figure
plot(results(:,1),excessCaCO3burial,results(:,1),excessSiO2burial,results(:,1),excessCaMgCO3burial)
ylabel('Excess burial (mol)')
xlabel('Time (years)')
legend('CaCO3','SiO2','CaMgCO32')

% Model dimenstions:
% Surface ocean: 100m
% Deep ocean: 3631m
% Atmosphere: dimensionless

% Tracers (all as concentrations in mol/m3): 
% dissolved silica concentration: surfSi, deepSi
% dissolved inorganic carbon (DIC): surfDIC, deepDIC, Catm (mol C)
% total alkalinity (TA): surfTA, deepTA
%
