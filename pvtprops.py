import numpy as np

def critical_presssure_Sutton(sg_gas):
    Ppc = 756.8 - 131*sg_gas - 3.6*(sg_gas**2)
    return Ppc

def critical_temperature_Sutton(sg_gas):
    Tpc = 169.2 + 349*sg_gas - 74*(sg_gas**2)
    return Tpc

def carr_kobayashi_temperature(sg_gas, CO2, H2S, N2):
    temperature_corrected = critical_temperature_Sutton(sg_gas) - (80*CO2) + (130*H2S) - (250*N2)
    return temperature_corrected

def carr_kobayashi_pressure(sg_gas, CO2, H2S, N2):
    pressure_corrected = critical_presssure_Sutton(sg_gas) + (440*CO2) + (600*H2S) - (170*N2)
    return pressure_corrected


def redlich_kwong_eos(P, T, sg_gas, CO2, H2S, N2):#INSERT T AS FAHRENHEIT
    R = 10.73
    a = (0.42727*(R**2)*(carr_kobayashi_temperature(sg_gas, CO2, H2S, N2)**2.5))/(carr_kobayashi_pressure(sg_gas, CO2, H2S, N2))
    b = 0.08664*(R*carr_kobayashi_temperature(sg_gas, CO2, H2S, N2))/(carr_kobayashi_pressure(sg_gas, CO2, H2S, N2))
    A = (a*P)/((R**2)*((T+460)**2.5))
    B = (b*P)/(R*(T+460))
    eqn = np.poly1d([1,-1,(A-B-(B**2)), -A*B])
    for i in eqn.roots:
        if np.imag(i) == 0 and i > 0:
            solution = i
    return solution.real

def Lee_Gonzales_Eakin(P,T,sg_gas,CO2,H2S,molecular_weight,method, N2):#INSERT T AS FAHRENHEIT
    density = molecular_weight*P/(redlich_kwong_eos(P,T,sg_gas,CO2,H2S,method,N2)*10.73*(T+460))
    X = 3.5 + (986/(T+460)) + 0.01*molecular_weight
    Y = 2.4 - 0.2*X
    K = (9.4+0.02*molecular_weight)*((T+460)**1.5)*(10**-4)/(209 + 19*molecular_weight + (T+460))
    viscosity = K*np.exp(X*((density/62.4)**Y))
    return viscosity


def calculate_Pb(Rsb,sg_gas,T,API):#StandingCorrelation
    Pb = 18.2*(((Rsb/sg_gas)**0.83)*10**(0.00091*T-0.0125*API)-1.4)
    return Pb

def calculate_GOR(Rsb,sg_gas,P,API,T):#Standing Correlation
    if P < calculate_Pb(Rsb,sg_gas,T,API):
        GOR = sg_gas*((((P/18.2)+1.4)*10**(0.0125*API-0.00091*T))**1.2048)
    else:
        GOR = Rsb
    return GOR

def calculate_FVF(Rsb,sg_gas,P,API,T):#Standing Correlation
    GOR = calculate_GOR(Rsb,sg_gas,P,API,T)
    FVF = 0.972+0.000147*((GOR*((sg_gas/calculate_sg_oil(API))**0.5)+1.25*T)**1.175)
    return FVF

def calculate_sg_oil(API):
    sg_oil = 141.5/(API+131.5)
    return sg_oil

