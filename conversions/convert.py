def psi_to_Pa(psi):
    Pa = psi*6894.757
    return Pa

def Pa_to_psi(Pa):
    psi = Pa/6894.757
    return psi

def m_to_ft(m):
    return m * 3.2808

def ft_to_m(ft):
    return ft / 3.2808

def in_to_m(L):
    Lm = (L/12)/3.2808
    return Lm
    
def H_feet_to_psi(H):
    Pw = H*62.4/144 #psi, head must be in feet 
    return Pw

def H_m_to_Pa(H):
    Pw = H*9806.38 #Pa, head must be in meters
    return Pw

def pressure_to_absolute(pressure, atm_pressure=101325):
    return pressure + atm_pressure

def pressure_to_gauge(pressure, atm_pressure=101325):
    return pressure - atm_pressure

def CMS_to_CMM(Q_cms):
    return Q_cms * 60

def CMM_to_CMS(Q_cmm):
    return Q_cmm / 60

def CFS_to_CFM(Q_cfs):
    return Q_cfs * 60

def CFM_to_CFS(Q_cfm):
    return Q_cfm / 60

def CFS_to_CMS(Q_cfs):
    return Q_cfs / (3.2808**3)

def CMS_to_CFS(Q_cms):
    return Q_cms * (3.2808**3)

def CFM_to_CMM(Q_cfm):
    return Q_cfm / (3.2808**3)

def CMM_to_CFM(Q_cmm):
    return Q_cmm * (3.2808**3)

def ft_to_m(ft):
    return ft / 3.2808

def C_to_F(C):
    return (C * (9/5)) + 32

def F_to_C(F):
    return (F - 32) * (5/9)
    
    