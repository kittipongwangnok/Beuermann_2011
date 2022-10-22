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

#Read datat
Cycle = []
T_obs = []
T_obs_err = []
#Please change the input file
for line in open("oc_dpleo_Beuermann_2011.inp"):
#for line in open("T_mid_Beuermann_Schwope.in"):
    li=line.strip()
    if not li.startswith("#"):
        Cycle.append(float(li.split(" ")[0]))
        T_obs.append(float(li.split(" ")[1]))
        T_obs_err.append(float(li.split(" ")[2]))
        
#Linear phemeris equation(From equation 1)
T0_bjd = 2454914.8322920
T0_bjd_err = 0.0000020
P0_day = 0.06236285648
P0_day_err = 0.00000000090


#calculate cycle (E)
#Array
BJD_time_a = [i for i in range(N_dpleo_bjd)]
Delta_aT = [i for i in range(N_dpleo_bjd)]
Delta_aT_err = [i for i in range(N_dpleo_bjd)]
E_af = [i for i in range(N_dpleo_bjd)] #float number
E_aj = [i for i in range(N_dpleo_bjd)] #integer number
P_aE = [i for i in range(N_dpleo_bjd)]
P_err_aE = [i for i in range(N_dpleo_bjd)]
T_aC = [i for i in range(N_dpleo_bjd)]
T_aO = [i for i in range(N_dpleo_bjd)]


OC_linear = []
print ('-----------------------------------------------------------------------------')
print ('Cycle \t\t T_O \t   T_C \t\t BJD - 2450000 \t OC \t OC_err')
print ('-----------------------------------------------------------------------------')
for i in range (0,N_dpleo_bjd):
    BJD_time = np.array(T_obs) - 2450000
    BJD_time_a[i] = BJD_time
    Delta_T = np.array(T_obs) - np.array(T0_bjd)
    E_k = Cycle
    E_aj[i] = E_k #arrays
##print (Delta_T)
    Delta_aT[i] = Delta_T #arrays
    Delta_T_err = np.sqrt((np.array(T_obs_err)/np.array(T_obs))**2 + (np.array(T0_bjd_err)/np.array(T0_bjd))**2)
    P_E_day = Delta_T[i] / E_k[i]
    P_aE[i] = P_E_day
    P_E_err_day = np.abs((np.array(T_obs_err[i]) - np.array(T0_bjd_err)) / E_k[i])
    P_err_aE[i] = P_E_err_day
    T_O = T0_bjd + P_aE[i]*E_k[i]
    T_aO[i] = T_O #arrays
#    print ('%0.6f' %(T_O))
    BJD_time = np.array(T_obs) - 2450000
    BJD_time_a[i] = BJD_time
    Delta_T = np.array(T_obs) - np.array(T0_bjd)
    #print (Delta_T)
    Delta_aT[i] = Delta_T #arrays
    Delta_T_err = np.sqrt((np.array(T_obs_err)/np.array(T_obs))**2 + (np.array(T0_bjd_err)/np.array(T0_bjd))**2)
    Delta_aT_err[i] = Delta_T_err #arrays
#    print (Delta_T_err[i])
    E_f = Delta_T / P0_day                      #Calculate cycle with float number
##    print (E_f)                                 #print cycle with float number
    E_af[i] = E_f #arrays
    E_j = np.round(Delta_T / P0_day)           #Calculate cycle with integer number
#print (E_j)                                #print cycle with integer number
    if  E_j[i] != 0:
        Delta_T[i] / E_j[i]
        P_E_day = Delta_T[i] / E_j[i]
        P_aE[i] = P_E_day
        P_E_err_day = np.abs((np.array(T_obs_err[i]) - np.array(T0_bjd_err)) / E_j[i])
#        print (P_E_err_day)
        P_err_aE[i] = P_E_err_day
        T_C = T0_bjd + P0_day*E_j[i]
        T_aC[i] = T_C #arrays
#    print (T_O, T_C)
        OC = np.array(T_O) - np.array(T_C)
        OC_s = (np.array(T_O) - np.array(T_C))*24*60*60
        OC_err = np.abs(np.sqrt((np.array(P_err_aE[i])**2 + (np.array(P0_day_err)**2))) * np.array(E_j[i]))
#        print (OC_err)
        OC_s_err = OC_err*24*60*60
    else:
        P_E_day = Delta_T[i] / E_k[i]
        P_aE[i] = P_E_day
        P_E_err_day = np.abs((np.array(T_obs_err[i]) - np.array(T0_bjd_err)) / E_k[i])
        P_err_aE[i] = P_E_err_day
        T_C = T0_bjd + P0_day*E_j[i]
        T_aC[i] = T_C #arrays
#    print (T_O, T_C)
        OC = np.array(T_O) - np.array(T_C)
        OC_s = (np.array(T_O) - np.array(T_C))*24*60*60
        OC_err = np.abs(np.sqrt((np.array(P_err_aE[i])**2)) *np.array(E_k[i]))
#        print (OC_err)
        OC_s_err = OC_err*24*60*60
#    print (P_E_day, P_E_err_day)
#    print (OC_s, OC_s_err)
#    print (Cycle[i], OC_s)
    print ('%0.0f\t%0.6f\t%0.6f\t%0.6f\t%0.2f\t%0.2f' %(Cycle[i], T_O, T_C, BJD_time[i], OC_s, OC_s_err))
    OC_linear.append('%0.0f\t%0.6f\t%0.6f\t%0.6f\t%0.2f\t%0.2f' %(Cycle[i], T_O, T_C, BJD_time[i], OC_s, OC_s_err))
    
rerults = OC_linear
f = open('oc_linear_dpleo_Beuermann_2011_rev1.txt', 'w')
#for upper_result in upper_result:
for i in range(len(rerults)):
    f.write(str(rerults[i])+ '\n')
f.close()


#Plot O-C vs BJD
InputFileInput  = 'oc_linear_dpleo_Beuermann_2011_rev1.txt'
Data   = np.genfromtxt(InputFileInput)
Cycle = Data[:,0]
T_O = Data[:,1]
T_C = Data[:,2]
BJD_time = Data[:,3]
OC_s = Data[:,4]
OC_s_err = Data[:,5]

##Plotgraph
fig=plt.figure(figsize=(10, 5))
plt.errorbar(BJD_time, OC_s, yerr=OC_s_err, fmt='o', color='limegreen')
plt.xlim(4850,5450)
plt.ylim(-5,5)
plt.xlabel('BJD - 2450000')
plt.ylabel('O-C (sec)')
plt.grid(linestyle='dotted')
######plt.title('O-C diagram: DP Leo')
output_filename = os.path.splitext(__file__)[0] + '.png'
plt.savefig(output_filename, dpi=1000)
plt.show()
#
#
