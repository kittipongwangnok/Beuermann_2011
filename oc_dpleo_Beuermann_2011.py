#Kittipong Wangnok
#School of Physics, Institute of Science, Suranaree University of Technology
#Import module
import sys
import os
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
from matplotlib import pyplot as plt

from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif',size=14)

#Please change the input file
DP_Leo_bjd = open('oc_dpleo_Beuermann_2011.inp','r').readlines()
#DP_Leo_bjd = open('T_mid_Beuermann_Schwope.in','r').readlines()
N_dpleo_bjd = len(DP_Leo_bjd)

T_obs = []
T_obs_err = []
#Please change the input file
for line in open("oc_dpleo_Beuermann_2011.inp"):
#for line in open("T_mid_Beuermann_Schwope.in"):
    li=line.strip()
    if not li.startswith("#"):
        T_obs.append(float(li.split(" ")[0]))
        T_obs_err.append(float(li.split(" ")[1]))
        
#Ephemeris (From equation 1)
T0_bjd = 2454914.8322920
T0_bjd_err = 0.0000020
P0_day = 0.06236285648
P0_day_err = 0.00000000090

#calculate cycle (E)
BJD_time_a = [i for i in range(N_dpleo_bjd)]
Delta_aT = [i for i in range(N_dpleo_bjd)]
Delta_aT_err = [i for i in range(N_dpleo_bjd)]
E_af = [i for i in range(N_dpleo_bjd)] #float number
E_aj = [i for i in range(N_dpleo_bjd)] #integer number
P_aE = [i for i in range(N_dpleo_bjd)]
P_err_aE = [i for i in range(N_dpleo_bjd)]
T_aC = [i for i in range(N_dpleo_bjd)]
T_aO = [i for i in range(N_dpleo_bjd)]
limit = 1500
for i in range (0,N_dpleo_bjd):
    BJD_time = np.array(T_obs) - 2450000
    BJD_time_a[i] = BJD_time
    Delta_T = np.array(T_obs) - np.array(T0_bjd)
#print (Delta_T)
    Delta_aT[i] = Delta_T #arrays
    Delta_T_err = np.sqrt((np.array(T_obs_err)/np.array(T_obs))**2 + (np.array(T0_bjd_err)/np.array(T0_bjd))**2)
#print (Delta_T_err)
    Delta_aT_err[i] = Delta_T_err #arrays
    E_f = Delta_T / P0_day
    E_af[i] = E_f #arrays
    E_j = Delta_T // P0_day
    E_aj[i] = E_j #arrays
#print (E_j)
    if E_j[i] != 0:
        P_E_day = Delta_T[i] / E_j[i]
    else:
        P_E_day = 0.06236285648
    P_aE[i] = P_E_day #arrays
#    print(P_E_day)
    if E_j[i] != 0:
        P_E_err_day = np.abs((np.array(T_obs_err[i]) - np.array(T0_bjd_err)) / E_j[i])
#print(P_E_err_day)
    else:
        P_E_err_day = 0
    P_err_aE[i] = P_E_err_day #arrays
#    print (P_E_err_day)
#    P_E_err_day = np.abs((np.array(Delta_T_err) - np.array(T0_bjd_err))/np.array(E_j))
#    P_err_aE[i] = P_E_err_day #arrays
#print (P_E_err_day)
#Calculate T from calculation: C = T_0+P_0*E
    T_C = T0_bjd + P0_day*E_j
    T_aC[i] = T_C #arrays
#print (T_C)
#Calculate T from observation: O = T_0+P_E*E
    T_O = T0_bjd + P_aE*E_j
    T_aO[i] = T_O #arrays
#print (T_O)    
    
#Calculate O-C
#O-C = T_O - T_C
OC_a = [i for i in range(N_dpleo_bjd)]
OC_err_a = [i for i in range(N_dpleo_bjd)]
OC_sa = [i for i in range(N_dpleo_bjd)]
OC_err_sa = [i for i in range(N_dpleo_bjd)]


num = 1000
#for i in range (0,N_dpleo_bjd):
OC = np.array(T_O) - np.array(T_C)
OC_a[i] = OC #arrays
n = len(OC)
OC_s = (np.array(T_O) - np.array(T_C))*24*60*60
OC_sa[i] = OC_s #arrays
#print(n)
for i in range (0,n):
    if OC_s[i] >= num:
            E_j = (Delta_T // P0_day) + 1
            E_aj[i] = E_j #arrays
            if E_j[i] != 0:
                P_E_day = Delta_T[i] / E_j[i]
            else:
                P_E_day = 0.06236285648
            P_aE[i] = P_E_day #arrays
            if E_j[i] != 0:
                P_E_err_day = np.abs((np.array(T_obs_err[i]) - np.array(T0_bjd_err)) / E_j[i])
            else:
                P_E_err_day = 0
            P_err_aE[i] = P_E_err_day #arrays
            #Calculate T from calculation: C = T_0+P_0*E
            T_C = T0_bjd + P0_day*E_j
            print (T_C)
            T_aC[i] = T_C #arrays
            #print (T_C)
            #Calculate T from observation: O = T_0+P_E*E
            T_O = T0_bjd + P_aE*E_j
            T_aO[i] = T_O #arrays
            #print (T_O)    
            OC = np.array(T_O) - np.array(T_C)
            OC_a = OC
            OC_err = np.abs(np.sqrt((np.array(P_err_aE)**2 + (np.array(P0_day_err)**2))) * np.array(E_j))
            OC_err_a = OC_err
            OC_s = (np.array(T_O) - np.array(T_C))*24*60*60
            OC_sa = OC_s
            OC_s_err = OC_err*24*60*60
            OC_err_sa = OC_s_err
    else:
        OC = np.array(T_O) - np.array(T_C)
        OC_a = OC
        OC_err = np.abs(np.sqrt((np.array(P_err_aE)**2 + (np.array(P0_day_err)**2))) * np.array(E_j))
        OC_err_a = OC_err
        OC_s = (np.array(T_O) - np.array(T_C))*24*60*60
        OC_sa = OC_s
        OC_s_err = OC_err*24*60*60
        OC_err_sa = OC_s_err


#print(OC_s)
#print(OC_s_err)

#type of txt file to save
np.savetxt('oc_dpleo_Beuermann_2011.out', np.c_[BJD_time, E_j, OC_s, OC_s_err], comments='#BJD, Cycle, OC_s, OC_s_err')

#Plotgraph
fig=plt.figure(figsize=(10, 5))
#plt.errorbar(E_j, OC_s, yerr=OC_s_err, fmt='o', color='limegreen')
plt.errorbar(BJD_time, OC_s, yerr=OC_s_err, fmt='o', color='limegreen')
#plt.ylim(-100.0,300.00)
plt.xlim(4850,5450)
plt.ylim(-5,5)
#Please change the number to generate x-axis
plt.xlabel('BJD - 2450000')
plt.ylabel('O-C (sec)')
plt.grid(linestyle='dotted')
#plt.title('O-C diagram: DP Leo')
output_filename = os.path.splitext(__file__)[0] + '.png'
plt.savefig(output_filename, dpi=1000)
plt.show()


