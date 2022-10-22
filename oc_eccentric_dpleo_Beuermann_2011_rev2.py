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
DP_Leo_bjd = open("oc_dpleo_Beuermann_Schwope_best.inp",'r').readlines()
N_dpleo_bjd = len(DP_Leo_bjd)

#Read datat
Cycle = []
T_obs = []
T_obs_err = []
#Please change the input file
for line in open("oc_dpleo_Beuermann_Schwope_best.inp"):
    li=line.strip()
    if not li.startswith("#"):
        Cycle.append(float(li.split(" ")[0]))
        T_obs.append(float(li.split(" ")[1]))
        T_obs_err.append(float(li.split(" ")[2]))

#Linear phemeris equation(From equation 1)
#T0_bjd = 2454914.8322920
#T0_bjd_err = 0.0000020
#P0_day = 0.06236285648
#P0_day_err = 0.00000000090


#Binary DP Leo ab:
T0_bjd = 2448773.21461
T0_bjd_err = 0.00009
P0_day = 0.0623628426
P0_day_err = 0.0000000006

#Giant Planet DP Leo (ab)c
P_c = 28                        #Orbital period (years)
a_c = 8.19                      #Semi-major axis (years)
e_c = 0.39                      #Eccentricity
Omega_c = -78                     #Longitude of periastron (Degree)
T_c = 2453025                   #Time of periastron passage
Kbin_c = 33.7                  #Amplitude of the LTT effect (s)
sini_c_abin_c = 1.01*10**12     #Semi-major axis (cm)
sini_c_M_c = 6.05               #Mass sin ic Mc
C = 3*10**8                     #unit m/s
c = C/(1*10**(-2))              #unit cm/s

print('Orbital period (years):', P_c)
print('Semi-major axis (years):', a_c)
print('Eccentricity:', e_c)
print('Longitude of periastron (Degree):', Omega_c)
print('Time of periastron passage:', T_c)
print('Amplitude of the LTT effect (s):', Kbin_c)
print('The speed of light (cm/s):', c)

#Kbin_c = sini_c_abin_c/c
#print (Kbin_c)

T_c_bjd = 2453025.004756446
P_c_yr = T_c*365.25
#Calculate the mean anomoly: M
#M = (2*(np.pi)/P_c_yr)*(T_c_bjd - T0_bjd)
M = (2*(np.pi)/P_c_yr)*(T_c - T0_bjd)
#M = (2*(np.pi)/P_c)*(T_c_bjd)
M_deg = (M*360)/(2*np.pi)
#M = T_c_bjd - T0_bjd
#M = (2*(360)/P_c)*(T_c - T0_bjd)
print('The Mean anomoly (rad):', M)
print('The Mean anomoly (Deg):', M_deg)

##################################################################
#The true anomaly: f
#InputFileInput  = 'EAno.txt'
#Data   = np.genfromtxt(InputFileInput)
#EAno = Data[:,0]
#print (EAno[0])
EAno = 0.00280
f1 = np.arccos((np.cos(EAno) - e_c)/(1 - e_c*np.cos(EAno)))
f2 = 2*np.arctan(np.sqrt((1+e_c)/(1-e_c))*np.tan(EAno/2))
print ('The true anomaly:', f1)
print ('The true anomaly:', f2)

#Light Travel Time (LTT) effect
LTT = Kbin_c*(1 - (e_c*e_c))/(1 + np.cos(f1))*np.sin(f1 - Omega_c)
print ('LTT value1 (s):', LTT)

LTT2 = Kbin_c*(1 - (e_c*e_c))/(1 + np.cos(f2))*np.sin(f2 - Omega_c)
print ('LTT value2 (s):', LTT2)

LTT3 = (sini_c_abin_c/c)*(1 - (e_c*e_c))/(1 + np.cos(f1))*np.sin(f1 - Omega_c)
print ('LTT value3 (s):', LTT3)

#Change LTT (s) to LTT (day)
LTT_day = LTT/(24*60*60)
print ('LTT value (day):', LTT_day)

#calculate cycle (E)
#Array
BJD_time_a = [i for i in range(N_dpleo_bjd)]
Delta_aT = [i for i in range(N_dpleo_bjd)]
Delta_aT_err = [i for i in range(N_dpleo_bjd)]
E_af = [i for i in range(N_dpleo_bjd)] #float number
E_ak = [i for i in range(N_dpleo_bjd)] #integer number
E_aj = [i for i in range(N_dpleo_bjd)] #integer number
P_aE = [i for i in range(N_dpleo_bjd)]
P_err_aE = [i for i in range(N_dpleo_bjd)]
T_aC_linear = [i for i in range(N_dpleo_bjd)]
T_aO_linear = [i for i in range(N_dpleo_bjd)]
T_aC_eccentric = [i for i in range(N_dpleo_bjd)]
T_aO_eccentric = [i for i in range(N_dpleo_bjd)]


OC_linear = []
#print ('-----------------------------------------------------------------------------')
#print ('Cycle \t\t T_O \t   T_C \t\t BJD - 2450000 \t OC_lin OC_err_Lin OC_occ')
#print ('-----------------------------------------------------------------------------')
for i in range (0,N_dpleo_bjd):
    BJD_time = np.array(T_obs) - 2450000
    BJD_time_a[i] = BJD_time
    Delta_T = np.array(T_obs) - np.array(T0_bjd)
    Delta_aT[i] = Delta_T #arrays
    Delta_T_err = np.sqrt((np.array(T_obs_err)/np.array(T_obs))**2 + (np.array(T0_bjd_err)/np.array(T0_bjd))**2)
    E_k = Cycle
    E_ak[i] = E_k #arrays
    #    print (Delta_T_err[i])
    E_f = Delta_T / P0_day                      #Calculate cycle with float number
    ##    print (E_f)                                 #print cycle with float number
    E_af[i] = E_f #arrays
    E_j = np.round(Delta_T / P0_day)           #Calculate cycle with integer number
##print (Delta_T)
    if  E_j[i] != 0:
        P_E_day = Delta_T[i] / E_j[i]
        P_aE[i] = P_E_day
        P_E_err_day = np.abs((np.array(T_obs_err[i]) - np.array(T0_bjd_err)) / E_j[i])
        P_err_aE[i] = P_E_err_day
        T_O_linear = T0_bjd + P_aE[i]*E_j[i]               #Linear
        T_O_eccentric = T0_bjd + P_aE[i]*E_j[i] + LTT_day      #Eccentric
        T_aO_linear[i] = T_O_linear #arrays
        T_aO_eccentric[i] = T_O_eccentric #arrays
    else:
        E_k[i] = 1
        P_E_day = Delta_T[i] / E_k[i]
#        print (P_E_day)
        P_aE[i] = P_E_day
        P_E_err_day = np.abs((np.array(T_obs_err[i]) - np.array(T0_bjd_err)) / E_k[i])
        P_err_aE[i] = P_E_err_day
        T_O_linear = T0_bjd + P_aE[i]*E_k[i]               #Linear
        T_O_eccentric = T0_bjd + P_aE[i]*E_k[i] + LTT_day        #Eccentric
        T_aO_linear[i] = T_O_linear #arrays
        T_aO_eccentric[i] = T_O_eccentric #arrays
#    print ('%0.6f' %(T_O))
#print (E_j)                                #print cycle with integer number
    if  E_j[i] != 0:
        P_E_day = Delta_T[i] / E_j[i]
        P_aE[i] = P_E_day
        P_E_err_day = np.abs((np.array(T_obs_err[i]) - np.array(T0_bjd_err)) / E_j[i])
#        print (P_E_err_day)
        P_err_aE[i] = P_E_err_day
        T_C_linear = T0_bjd + P0_day*E_j[i]              #Linear
        T_C_eccentric = T0_bjd + P0_day*E_j[i] + LTT_day     #Eccentric
        T_aC_linear[i] = T_C_linear #arrays
        T_aC_eccentric[i] = T_C_eccentric #arrays
#    print (T_O, T_C)
        OC = np.array(T_O_linear) - np.array(T_C_linear)
        OC_s = (np.array(T_O_linear) - np.array(T_C_linear))*24*60*60
        OC_err = np.abs(np.sqrt((np.array(P_err_aE[i])**2 + (np.array(P0_day_err)**2))) * np.array(E_j[i]))
#        print (OC_err)
        OC_s_err = OC_err*24*60*60
        OC_ecc = np.array(T_O_eccentric) - np.array(T_C_eccentric)
        OC_s_ecc = (np.array(T_O_eccentric) - np.array(T_C_eccentric))*24*60*60
    else:
        P_E_day = Delta_T[i] / E_k[i]
        P_aE[i] = P_E_day
        P_E_err_day = np.abs((np.array(T_obs_err[i]) - np.array(T0_bjd_err)) / E_k[i])
        P_err_aE[i] = P_E_err_day
        T_C_linear = T0_bjd + P0_day*E_j[i]              #Linear
        T_C_eccentric = T0_bjd + P0_day*E_j[i] + LTT_day     #Eccentric
        T_aC_linear[i] = T_C_linear #arrays
        T_aC_eccentric[i] = T_C_eccentric #arrays
#    print (T_O, T_C)
        OC = np.array(T_O_linear) - np.array(T_C_linear)
        OC_s = (np.array(T_O_linear) - np.array(T_C_linear))*24*60*60
        OC_err = np.abs(np.sqrt((np.array(P_err_aE[i])**2)) *np.array(E_k[i]))
#        print (OC_err)
        OC_s_err = OC_err*24*60*60
#        OC_ecc = np.array(T_O_eccentric) - np.array(T_C_eccentric)
#        OC_s_ecc = (np.array(T_O_eccentric) - np.array(T_C_eccentric))*24*60*60
#        Res = (T1 - T2)
#        print ('%0.10f' %(Res))
#    print (P_E_day, P_E_err_day)
#    print (OC_s, OC_s_err)
 #   print ('%0.6f\t%0.6f\t%0.6f' %(T_O_linear, T_C_linear, OC_s))
#    print ('%0.6f\t%0.6f\t%0.6f' %(T_O_eccentric, T_C_eccentric, OC_s_ecc))
    print ('%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f' %(T_O_linear, T_C_linear, OC_s, T_O_eccentric, T_C_eccentric,  OC_s_ecc))
#    Res = (OC_s - OC_s_ecc)
#    print ('%0.10f' %(Res))
#    print ('%0.6f\t%0.6f' %(T_O_eccentric, T_C_eccentric))
#    T_OC = T_O_eccentric - T_C_eccentric
#    print (T_OC)
#    print (OC_ecc)
#    print (OC_s_ecc)
#    print (Cycle[i], OC_s)
#    print ('%0.0f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.5f\t%0.2f\t%0.5f' %(Cycle[i], T_O_linear, T_O_eccentric, T_C_linear, T_C_eccentric, BJD_time[i], OC_s, OC_s_err, OC_s_ecc))
    OC_linear.append('%0.0f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.2f\t%0.2f\t%0.2f' %(Cycle[i], T_O_linear, T_O_eccentric, T_C_linear, T_C_eccentric, BJD_time[i], OC_s, OC_s_err, OC_s_ecc))
 
    
    #Cylce:E
#    E = np.abs(Cycle[i])
#    E = Cycle[i]
#    print (E)
    #The eccentric anmoly: e
#    e = (M - 0.00280)/np.sin(0.00280)
#    e = np.sin(E)
#    print ('The eccentric anmoly:', e)
    


rerults = OC_linear
f = open('oc_eccentric_dpleo_Beuermann_2011_rev1.txt', 'w')
#for upper_result in upper_result:
for i in range(len(rerults)):
    f.write(str(rerults[i])+ '\n')
f.close()


#Plot O-C vs BJD
InputFileInput  = 'oc_eccentric_dpleo_Beuermann_2011_rev1.txt'
Data   = np.genfromtxt(InputFileInput)
#Cycle = Data[:,0]
#T_O = Data[:,1]
#T_C = Data[:,2]
#BJD_time = Data[:,3]
#OC_s = Data[:,4]
#OC_s_err = Data[:,5]

Cycle = Data[:,0]
T_O_linear = Data[:,1]
T_O_eccentric = Data[:,2]
T_C_linear = Data[:,3]
T_C_eccentric = Data[:,4]
BJD_time = Data[:,5]
OC_s = Data[:,6]
OC_s_err = Data[:,7]
OC_s_ecc = Data[:,8]


#Plotgraph
fig=plt.figure(figsize=(10, 5))
plt.errorbar(BJD_time, OC_s, yerr=OC_s_err, fmt='o', color='limegreen')
#plt.plot(BJD_time, OC_s_ecc, lw=2, color='black')
plt.xlim(-6000,6000)
plt.ylim(-40,60)
plt.xlabel('BJD - 2450000')
plt.ylabel('O-C (sec)')
plt.grid(linestyle='dotted')
#####plt.title('O-C diagram: DP Leo')
output_filename = os.path.splitext(__file__)[0] + '.png'
plt.savefig(output_filename, dpi=1000)
plt.show()
#
#
