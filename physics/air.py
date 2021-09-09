def rho_air(Pabs,T):
    #P in Pa, T in C
    Rspecific = 287.058 #J/(kg*K) specific gas constant for dry air
    TK = T + 273.15
    return Pabs/(Rspecific*TK)

def Q(mdot, rho_air): #CMS
    return mdot / rho_air

def mdot(Q, rho_air):
    return Q * rho_air

