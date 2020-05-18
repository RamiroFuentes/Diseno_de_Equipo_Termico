import numpy as np
from colorama import init, Fore, Back, Style

# Calculo del flujo de calor
def heat_flux(m,cp,T1,T2):
    Q = m*cp*abs(T2-T1)*(1/3600)
    return Q

# Calculo de diferencia de temperatura logaritmica 
def TempML(T1,T2,t1,t2):
    Tml = ((T1-t2)-(T2-t1))/np.log((T1-t2)/(T2-t1))
    return Tml

### Calculo cambio de temperatura real 

# Comenzamos con el calculo de los factores R y S
def FactorR(T1,T2,t1,t2):
    R = (T1-T2)/(t2-t1)
    return R

def FactorS(T1,T2,t1,t2):
    S = (t2-t1)/(T1-t1)
    return S

# Procedemos con el calculo del factor de corrección 
def CorrectionFT(R,S,Np):
    "Por facilidad interpretaremos la ecuacion por componentes"
    fact = (((R**2)+1)**(.5))
    
    top = fact*np.log((1-S)/(1-R*S))
    bottom = (R-1) * np.log( (Np-S*(R+1-fact)) / (Np-S*(R+1+fact) ))
    
    Ft = top/bottom
                           
    return Ft

# Area de transferencia de calor
def TotalA(Q,Up,Tr):
    "En esta ecuación hay una corrección de magnitud de *10^-3"
    At = (Q*10**3)/(Up*Tr)
    return At

# Para el area de un solo tubo 
def TubeArea(do,L):
    "En esta ecuación hay una corrección de magnitud de *10^-3"
    TbA = 2*np.pi*(do/2)*10**-3*L
    return TbA

#Calculo del flujo másico del segundo servicio 
def MassFlux(Q,Cp,t2,t1):
    m2 = Q/(Cp*(t2-t1))
    return  m2

# velocidad del fluido en cada tubo 
def u_Tube(mTube,rho,di):
    "En esta ecuación hay una corrección de magnitud de *10^-3"
    A = (np.pi*(di*10**-3)**2)/4
    u_tube = mTube / (rho*A )
    return u_tube

# Numeros adimensionales para los tubos 
def Reynolds(rho,u,di,Mu):
    "En esta ecuación hay una corrección de magnitud de *10^-3"
    Re = rho*u*(di*10**-3)/Mu
    return Re

def Prandtl(Cp,Mu,k):
    "En esta ecuación hay una corrección de magnitud de *10^3"
    Pr = Cp*10**3*Mu/k
    return Pr

# Para determinar el coeficiente de trasferencia de calor para los tubos 
"Se hara uso de un diagrama el cual se espera posteriormente encontrar su polinomio 12.23"

def hT_coefficient(JH,Re,Pr,kf,di):
    hT = ( JH * Re * ((Pr)**(1/3)) * kf ) / (di*10**-3)
    return hT

# Caida de presion para los tubos 
"Se hara uso de un diagrama el cual se espera posteriormente encontrar su polinomio 12.24"

def pressure_drop(Np,JF,L,di,rho,u):
    DP = Np * (8*JF*(L/(di*(10**-3)))+2.5)*((rho*(u**2))/2) 
    return DP

# Area de la coraza
def As_shell(Pt,do,Ds,Lb):
    As = (Pt-(do*10**-3))*Ds*Lb/Pt
    return As


def Solve_ICTyC(Services_Data,Physical_properties,Exchanger_Data,Reynolds_Factors):
    # Limpieza de información
    
    m = Services_Data["m"] # flujo masico del agua 
    t1 = Services_Data["t1"] # Temperaturas
    t2 = Services_Data["t2"]
    T1= Services_Data["T1"]
    T2= Services_Data["T2"]
    Up = Services_Data["Up"] # Coeficiente global de transferencia de calor propuesto

    Cp = Physical_properties["Cp"] # Calor específico
    Mu = Physical_properties["Mu"] # Viscocidad
    k = Physical_properties["k"] # Coeficiente de trasferencia de calor 
    rho = Physical_properties["rho"] # Densidad
    had = Physical_properties["had"]
    hid = Physical_properties["hid"]

    do = Exchanger_Data["do"] # Diametro exterior
    t = Exchanger_Data["t"] # Espesor del tubo 
    di = Exchanger_Data["di"] # Diametro interno
    L = Exchanger_Data["L"]# Longitud de un tubo
    Np = Exchanger_Data["Np"] # numero de pasos totales por la coraza 
    K = Exchanger_Data["K"]
    
    # Cálculos simplificados:
    Q = heat_flux(m,Cp,T1,T2)
    print(f'Cantidad de calor transferido: {Q} KW')

    Tml = TempML(T1,T2,t1,t2)
    print(f'Temperatura media logaritmica: {Tml} ºC')

    R = FactorR(T1,T2,t1,t2)
    print(f'Factor R = {R}')

    S = FactorS(T1,T2,t1,t2)
    print(f'Factor S = {S}')

    Ft = CorrectionFT(R,S,Np)
    print(f'Factor Ft = {Ft}')

    # Para la temperatura real
    Tr = Tml*Ft
    print(f'Temperatura real = {Tr} ºC')

    # Area total de transferencia
    At = TotalA(Q,Up,Tr)
    print(f'Area total = {At} m^2') 

    #Area de transferencia por tubo 
    TbA = TubeArea(do,L)
    print(f'Area por un solo tubo = {TbA} m^2')

    # Calculo del numero total de tubos 
    NTb = At/TbA
    print(f'Numero total de tubos = {NTb}')

    #Flujo masico segundo servicio
    m2 = MassFlux(Q,Cp,t2,t1)
    print(f'Flujo masico 2 = {m2} Kg/s')

    #Flujo masico por tubo 
    mTube = m2/(NTb/Np)
    print(f'Flujo masico por tubo = {mTube} Kg/s')

    # Velocidad en los tubos
    ut = u_Tube(mTube,rho,di)
    print(f'Velocidad en cada tubo= {ut} m/s')

    # Reynolds en los tubos 
    Tube_Re = Reynolds(rho,ut,di,Mu)
    print(f'Reynolds en los tubos = {Tube_Re}')

    # Prandlt en los tubos
    Pr = Prandtl(Cp,Mu,k)
    print(f'Valor de Prandtl = {Pr}')

    "Tubos - Se hara uso de un diagrama el cual se espera posteriormente encontrar su polinomio 12.23"
    JH = Reynolds_Factors["JH"]

    # Calculo del coeficiente de transferencia del calor en los tubos
    hT = hT_coefficient(JH,Tube_Re,Pr,k,di)
    print( f"Coeficiente de tranferencia de calor para tubos = {hT_coefficient(JH,Tube_Re,Pr,k,di)} W/m^2ºC")

    "Tubos - Se hara uso de un diagrama el cual se espera posteriormente encontrar su polinomio 12.24"
    JF = Reynolds_Factors["JF"]

    # Calclo de caídas de presión
    DP = pressure_drop(Np,JF,L,di,rho,ut)
    print(f'Caida de presion en los tubos = {DP} Pa')

    # Para determinar el coeficiente de trasferencia de calor para la coraza y la caida de presion
    # se necesitan calcular unos datos antes

    # Pitch
    Pt = 1.25*do*10**-3
    print(f'Pitch = {Pt} m')

    # Diametro del banco de tubos
    "Se usaran datos obtenidos de tablas"
    k1 = .249
    m1 = 2.207

    Db = (do*10**-3)*(NTb/k1)**(1/m1)
    print(f'Diametro de banco de tubos={Db} m')

    # Diametro interno de la coraza
    "De la grafica bundle diameter vs shell inside se obtiene un valor de BD"
    BD = 65*10**-3

    Ds = Db+BD
    print(f'Diametro Coraza = {Ds} m')

    # Separacion de los bafles
    Lb = 0.4*Ds
    print(f'Separacion de baffles = {Lb} m')

    # Area de la coraza
    As = As_shell(Pt,do,Ds,Lb)
    print(f'Area de la coraza = {As} m^2')

    # Velocidad masica
    "Conversion de horas a segundos en el flujo masico"
    Gs = m*(1/3600)/As
    print(f'Velocidad másica = {Gs} Kg/m^2s')

    # Velocidad lineal coraza
    us = Gs/rho
    print(f"Velocidad lineal coraza = {us} m/s")

    # Diametro equivalente para el arregl triangular
    de = ( 1.10 / (do*10**-3) ) * ((Pt**2)-(0.917*(do*10**-3)**2))
    print(f"Diametro equivalente = {de}")

    # Calculo del reynolds para coraza
    ReS = (Gs*de)/Mu
    print(f"Reynolds coraza = {ReS}")

    # Calculo del Prandtl Coraza
    PrS = Cp*10**3*Mu/k
    print(f"Prandtl coraza = {PrS}")

    # Para determinar el coeficiente de trasferencia de calor para la coraza
    " Coraza - Usando constantes dadas en tablas JH y JF"
    JHs = Reynolds_Factors["JHs"]
    JFs = Reynolds_Factors["JFs"]

    hs = JHs*ReS*(PrS**(1/3))*k/de
    print(f"Coeficiente de transferencia de calor coraza = {hs} W/m^2ºC")

    # Caida de presion para la coraza
    DPs = 8 * JFs * (Ds/de) * (L/Lb) * ((rho*(us**2)) / 2)
    print(f"Caida de presion en la coraza = {DPs} Pa")

    # Calculo del coeficiente global de transderencia de calor 
    "Se obtienen los valores del factor de ensuciamiento y conductividad termica"
    Denominador = (1/hs) + (1/had) + ((do/1000)*np.log(do/di)/(2 * K)) + (((do)/(di))*(1/hid)) + ((do/di)*(1/hT))

    U_Global = 1/Denominador

    print(f"Coeficiente global de transferencia de calor: {U_Global} W/m^2ºC")

    print("")

    print(f'Caida de presion en los tubos = {DP*10**-3} KPa')
    print(f'Reynolds en los tubos = {Tube_Re}')
    print(f"Caida de presion en la coraza = {DPs*10**-3} KPa")
    print(f"Reynolds coraza = {ReS}")
    print(f"Coeficiente global de transferencia de calor: {U_Global} W/m^2ºC")
    return