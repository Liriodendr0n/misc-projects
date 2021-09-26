function [configurations] = evalConfig(configurations)
%evalConfig checks if each configuration meets constraints/requirements
%   Parameters given in each configuration are taken in, and combined with
%   known constants (total energy, TOFL) to flesh out each configuration
%   the configs are then checked against the requirements

Nconfigs = length(configurations);
%% main loop
for i = 1:Nconfigs
    
    %% config data
    S = configurations(i).S;
    b = configurations(i).b;
    Nbox = configurations(i).Nbox;
    CLmax = configurations(i).CLmax;
    EWF = configurations(i).EWF;
    P = configurations(i).P;
    
    AR = b^2/S;
    configurations(i).AR = AR;
    e0 = 1/(1.05 + 0.007*pi*AR);
    configurations(i).e0 = e0;
    
    A = 2.18; % prop area
    
    W = 0.5*Nbox/(1-EWF);
    WS = W/S;
    m = W/32.2;
    
    configurations(i).W = W;
    
    % CD0
    rho = 0.00238;  %slug ft^-3
    mu = 3.62*10^-7;   %lbf s ft^-2
    
    Vcrmin = 10*Nbox/2;
    ReL = rho*Vcrmin/mu;

    L_fus = Nbox*(3.5/12)+0.5;
    Swet_fus = 4*(L_fus*0.5);
    Cf_fus = 0.455./(log10(ReL*L_fus).^2.58);
    A_fus = 0.25;
    Cd0_fus = Swet_fus./S .*Cf_fus + A_fus./S * 0.8;
    
    Cf_wing = 0.455./(log10(ReL*S/b).^2.58);
    Cd0_wing = 2*S/S *(Cf_wing);
    
    CD0 = Cd0_fus + Cd0_wing + 0.01;
    configurations(i).CD0 = CD0;
    
    % takeoff (thrust from TW, drag from drag at 71% VTO)
    
    CLTO = CLmax*0.8;
    VTO = sqrt(2*W/(rho*S*CLTO)); 
    CDTO = CD0 + CLTO^2/(pi*e0*AR);
    configurations(i).VTO = VTO;
    LDTO = CLTO/CDTO;
    configurations(i).LDTO = LDTO;
    
    TTO = 0.6*(2*P^2*rho*A)^(1/3);
    DTO = 0.5*rho*(sqrt(2)/2 *VTO)^2*S*(CLmax^2/(pi*e0*AR)) + 0.5*rho*(sqrt(2)/2 * VTO)^2*S*CD0;
    %TTO = W*TW - 0.5*rho*(sqrt(2)/2 *VTO)^2*S*(CLmax^2/(pi*e0*AR)) - 0.5*rho*(sqrt(2)/2 * VTO)^2*S*CD0;
    
    FTO = TTO - DTO;
    
    TOFL = 0.5*(FTO/m)*(VTO/(FTO/m))^2;
    configurations(i).TOFL = TOFL;
    configurations(i).TOTW = FTO/W;
    
    % M3 cruise
    VmaxR = (4 * (WS)^2 * 1/(rho^2) * 1/CD0 * 1/(pi*e0*AR))^(1/4);
    Vavg = 1.5*max([10 * 0.5*Nbox, 30, VmaxR]);
    configurations(i).Vavg = Vavg/1.5;
    CLcr = (2*W)/(rho*Vavg^2*S);
    CDcr = CD0 + CLcr^2/(pi*e0*AR);
    LDtr = CLcr/CDcr;
    configurations(i).LDcr = LDtr;
    Preq = 0.5*rho*Vavg^3*S*CD0 + (2*W^2)/(rho*Vavg*S) * (1/(pi*e0*AR));
    
    Range = max([3*3000, Nbox*3000]);
    treq = min([Range/(Vavg), 10*60]);
    configurations(i).treq = treq;
    
    Ereq = Preq*treq * 0.0003766161; % convert from lbf*ft to Wh
    
    % M2 dash cruise (max power cruise)
    R = 9000;
    %Etot = 100/0.0003766161 * 0.7;
    
    a1 = 0.5*rho*S*CD0;
    b1 = 2*W^2/(rho*S*pi*e0*AR);
    
    Vdash = (0.7*P/a1)^(1/3);
    configurations(i).Vdash = Vdash;
    tM2 = R/Vdash;
    CLda = (2*W)/(rho*Vdash^2*S);
    CDda = CD0 + CLda^2/(pi*e0*AR);
    LDda = CLda/CDda;
    configurations(i).LDda = LDda;
    
    Ereqdash = 0.0003766161*tM2*0.5*rho*Vdash^3*S*CD0 + (2*W^2)/(rho*Vdash*S) * (1/(pi*e0*AR));
    
    
    configurations(i).Ereq = Ereq;
    configurations(i).Ereqdash = Ereqdash;
    
    configurations(i).Range = Range;
    
    Ereq = max([Ereq, Ereqdash]);
    
    
    configurations(i).M2score = 10*Nbox/tM2 + rand(1)*0.001;
    configurations(i).M3score = Nbox + rand(1)*0.001;


    %disp("working config: " + i + "/" + Nconfigs)
    
    %% validate configuration
    isTOFL = TOFL < 25 && 0 < TOFL; % takeoff in less than 25 feet
    isEtot = Ereq < 100; % requires less than 260000 ft*lb (100 Wh) and efficiency
    isW = W < 55; % weighs less than 55 pounds
    
    configurations(i).isValid = all([isTOFL, isEtot, isW]);
    if configurations(i).isValid == 0
        configurations(i).M2score = 0;
        configurations(i).M3score = 0;
    end

%% end main loop
end


end

