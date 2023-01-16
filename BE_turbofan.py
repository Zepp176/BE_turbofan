import matplotlib.pyplot as plt
import numpy as np

def compute(Tt4=1600, pi_OPR=40, lam=11, pi_f=1.4):
    
    # Data
    
    gamma = 1.4
    gamma_s = 1.33
    r = 287
    r_s = 291.6
    cp = gamma*r/(gamma-1)
    cp_s = gamma_s*r_s/(gamma_s-1)
    P = 21000
    pi_CHP = 22
    pi_CBP = pi_OPR/(pi_CHP*pi_f)
    PCI = 42800000
    M0 = 0.78
    p0 = 22700
    T0 = 217
    
    pc_02 = 0.98
    eta_c = 0.90
    eta_f = 0.92
    eta_comb = 0.99
    pc_34 = 0.95
    eta_m = 0.98
    eta_THP = 0.89
    eta_TBP = 0.90
    pc_tuyere = 0.98
    
    F_rayon=0.3

    # 0
    
    Tt0 = T0*(1 + (gamma-1)/2 * M0**2)
    pt0 = p0*(1 + (gamma-1)/2 * M0**2)**(gamma/(gamma-1))
    V0 = M0*np.sqrt(gamma*r*T0)
    
    #print("etat 0 \t", pt0, Tt0)
    
    # 2
    
    Tt2 = Tt0
    M2 = 0.6
    T2 = Tt2/(1 + (gamma-1)/2 * M2**2)
    pt2 = pc_02*pt0
    p2 = pt2/((1 + (gamma-1)/2 * M2**2)**(gamma/(gamma-1)))
    V2 = M2*np.sqrt(gamma*r*T2)
    
    #print("etat 2 \t", pt2, Tt2)
    
    # 21
    
    pt21 = pi_f*pt2
    Tt21 = Tt2*pi_f**((gamma-1)/(gamma*eta_f))
    
    #print("etat 21\t", pt21, Tt21)
    
    # 25
    
    pt25 = pi_CBP*pt21
    Tt25 = Tt21*pi_CBP**((gamma-1)/(gamma*eta_c))
    
    #print("etat 25\t", pt25, Tt25)
    
    # 3
    
    pt3 = pi_CHP*pt25
    Tt3 = Tt25*pi_CHP**((gamma-1)/(gamma*eta_c))
    
    #print("etat 3 \t", pt3, Tt3)
    
    # 4 - combustion
    
    pt4 = pt3*pc_34
    alpha = (cp_s*Tt4-cp*Tt3)/(PCI*eta_comb - Tt4*cp_s)
    
    #print("etat 4 \t", pt4, Tt4)
    
    # 45
    
    Tt45 = Tt4 - cp/((1+alpha)*cp_s*eta_m) * (Tt3 - Tt25)
    pt45 = pt4*(Tt45/Tt4)**(gamma_s/((gamma_s-1)*eta_THP))
    
    #print("etat 45\t", pt45, Tt45)
    
    # 5
    
    Tt5 = Tt45 - cp/((1+alpha)*cp_s*eta_m) * ((1+lam)*(Tt21 - Tt2) + (Tt25 - Tt21))
    pt5 = pt45*(Tt5/Tt45)**(gamma_s/((gamma_s-1)*eta_TBP))
    
    #print("etat 5 \t", pt5, Tt5)
    
    # 9
    
    pt9 = pt5*pc_tuyere
    Tt9 = Tt5
    M9 = np.sqrt( 2/(gamma_s-1) * ((pt9/p0)**((gamma_s-1)/gamma_s) - 1) )
    T9 = Tt9 * (1 + (gamma_s-1)/2*M9**2)**-1
    V9 = M9*np.sqrt(gamma_s*r_s*T9)
    F9_sp = ((1+alpha)*V9 - V0)/(lam+1)
    
    #print("etat 9 \t", pt9, Tt9)
    
    # 19 - tuyère flux secondaire
    
    pt19 = pt21*pc_tuyere
    Tt19 = Tt21
    M19 = np.sqrt( 2/(gamma-1) * ((pt19/p0)**((gamma-1)/gamma) - 1) )
    T19 = Tt19 * (1 + (gamma-1)/2*M19**2)**-1
    V19 = M19*np.sqrt(gamma*r*T19)
    F19_sp = lam*(V19-V0)/(1+lam)
    
    #print("etat 19\t", pt19, Tt19)
    
    # débits
    
    m_dot = P/(F9_sp + F19_sp)
    m_p = m_dot/(1+lam)
    m_s = m_dot*lam/(1+lam)
    m_k = alpha*m_p
    TSFC = m_k/P
    
    #print(m_dot)
    #print("TSFC : {:.5} g/Kn*s".format(TSFC*1e6))
    
    # Entrée d'air
    
    rho2 = p2/(r*T2)
    S2 = m_dot/(rho2*V2)
    rmax = np.sqrt(S2/(np.pi*(1-F_rayon**2)))
    Dmax = 2*rmax
    
    # Performances
    
    P_chim = m_k*PCI
    P_cy = 1/2*(m_p*((1+alpha)*V9**2-V0**2) + m_s*(V19**2-V0**2))
    P_pr = P*V0
    
    eta_th = P_cy/P_chim
    eta_pr = P_pr/P_cy
    
    eta_glob = eta_th * eta_pr
    
    return eta_th, eta_pr

#%%

pi_OPR = np.linspace(20, 70, 100)
Tt4 = np.array([1300, 1400, 1500, 1600, 1700])
eta_glob = np.empty((len(Tt4), len(pi_OPR)))
eta_th = np.empty((len(Tt4), len(pi_OPR)))
eta_pr = np.empty((len(Tt4), len(pi_OPR)))

for j, T in enumerate(Tt4):
    for i, pi in enumerate(pi_OPR):
        eta_th[j, i], eta_pr[j, i] = compute(T, pi)    
        eta_glob[j, i] = eta_th[j, i]*eta_pr[j, i]
   
#%%

plt.figure(figsize=(8,6))
for j in range(len(Tt4)):
    plt.plot(pi_OPR, 100*eta_glob[j, :])
plt.grid()
plt.title("Rendement global en fonction du rapport de compression")
plt.xlabel("$\pi_{OPR}$ [-]")
plt.ylabel("$\eta_{global}$ [%]")
plt.legend(["T_t4 = {}".format(T) for T in Tt4])
plt.tight_layout()
plt.savefig("rend_glob.png", dpi=300)

plt.figure(figsize=(8,6))
for j in range(len(Tt4)):
    plt.plot(pi_OPR, 100*eta_th[j, :])
plt.grid()
plt.title("Rendement thermique en fonction du rapport de compression")
plt.xlabel("$\pi_{OPR}$ [-]")
plt.ylabel("$\eta_{thermique}$ [%]")
plt.legend(["T_t4 = {}".format(T) for T in Tt4])
plt.tight_layout()
plt.savefig("rend_ther.png", dpi=300)

plt.figure(figsize=(8,6))
for j in range(len(Tt4)):
    plt.plot(pi_OPR, 100*eta_pr[j, :])
plt.grid()
plt.title("Rendement propulsif en fonction du rapport de compression")
plt.xlabel("$\pi_{OPR}$ [-]")
plt.ylabel("$\eta_{propulsif}$ [%]")
plt.legend(["T_t4 = {}".format(T) for T in Tt4])
plt.tight_layout()
plt.savefig("rend_prop.png", dpi=300)

#%%

lam = np.linspace(5, 20, 100)
pi_f = np.array([1.2, 1.4, 1.6, 1.8])
eta_glob = np.empty((len(pi_f), len(lam)))
eta_th = np.empty((len(pi_f), len(lam)))
eta_pr = np.empty((len(pi_f), len(lam)))

for j, pi in enumerate(pi_f):
    for i, l in enumerate(lam):
        eta_th[j, i], eta_pr[j, i] = compute(1600, 40, l, pi)    
        eta_glob[j, i] = eta_th[j, i]*eta_pr[j, i]

#%%

plt.figure(figsize=(8,6))
for j in range(len(pi_f)):
    plt.plot(lam, 100*eta_glob[j, :])
plt.grid()
plt.title("Rendement global en fonction du taux de dilution")
plt.xlabel("BPR [-]")
plt.ylabel("$\eta_{global}$ [%]")
plt.legend(["pi_f = {}".format(pi) for pi in pi_f])
plt.tight_layout()
plt.savefig("rend2_glob.png", dpi=300)

plt.figure(figsize=(8,6))
for j in range(len(pi_f)):
    plt.plot(lam, 100*eta_th[j, :])
plt.grid()
plt.title("Rendement thermique en fonction du taux de dilution")
plt.xlabel("BPR [-]")
plt.ylabel("$\eta_{thermique}$ [%]")
plt.legend(["pi_f = {}".format(pi) for pi in pi_f])
plt.tight_layout()
plt.savefig("rend2_ther.png", dpi=300)

plt.figure(figsize=(8,6))
for j in range(len(pi_f)):
    plt.plot(lam, 100*eta_pr[j, :])
plt.grid()
plt.title("Rendement propulsif en fonction du taux de dilution")
plt.xlabel("BPR [-]")
plt.ylabel("$\eta_{propulsif}$ [%]")
plt.legend(["pi_f = {}".format(pi) for pi in pi_f])
plt.tight_layout()
plt.savefig("rend2_prop.png", dpi=300)