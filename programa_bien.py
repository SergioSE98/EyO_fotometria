# -*- coding: utf-8 -*-
"""
Created on Sat Apr 10 13:19:34 2021

@author: Sergio
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 15:10:46 2021

@author: Sergi
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

#Leo los datos

J_data = np.loadtxt('J_trans.dat', skiprows=2)  
J_long = J_data[:,0]                                          #Long es de longitud de onda
J_trans = J_data[:,1]

H_data = np.loadtxt('H_trans.dat', skiprows=2)
H_long = H_data[:,0]
H_trans = H_data[:,1]

K_data = np.loadtxt('Kshort_trans.dat', skiprows=2)
K_long = K_data[:,0]
K_trans = K_data[:,1]

ukk5_data = np.loadtxt('ukk5iii.dat', skiprows=2)
ukk5_long_uds = ukk5_data[:,0]
ukk5_normflux = ukk5_data[:,1]  #Esto lo pasaré a fotones al final para graficar
ukk5_long=ukk5_long_uds/10000    #las uds de long vienen 10000 veces mayores que las del resto (paso a micrometros)

uka0_data = np.loadtxt('uka0v.dat', skiprows=2)
uka0_long_uds = uka0_data[:,0]
uka0_normflux = uka0_data[:,1]
uka0_long=uka0_long_uds/10000

vega_data = np.loadtxt('vegfluxtot.dat', skiprows=2)
vega_long = vega_data[:,0]
vega_absflux = vega_data[:,1]

trans_filts_array = np.array([J_trans, H_trans, K_trans])

#Calculo también el salto entre longitudes de onda que se da para cada filtro 

dist_J_long = J_long[1]-J_long[0]
dist_H_long = H_long[1]-H_long[0]
dist_K_long = K_long[1]-K_long[0]

dist_longs_array = np.array([dist_J_long, dist_H_long, dist_K_long])

#Ploteo los datos para ver que los he importado correctamente
"""
plt.figure()
plt.plot(J_long,J_trans)
plt.plot(H_long,H_trans)
plt.plot(K_long,K_trans)
plt.figure()
plt.plot(ukk5_long,ukk5_normflux)
plt.plot(uka0_long,uka0_normflux)
plt.figure()
plt.plot(vega_long,vega_absflux)
"""
#Datos
trans_sys = 0.3
diam = 10
h = 6.62607015e-34
c = 299792458e6   #Velocidad luz en uds de long correctas (no m/s, sino micrometros/s)
area=np.pi*(diam/2)**2
airmass = 1.5
ext_array = np.array([0.1,0.15,0.2])

m_array = np.array([10,15,20])

m_01 = m_array + airmass*ext_array[0] 
m_015 = m_array + airmass*ext_array[1] 
m_02 = m_array + airmass*ext_array[2] 

#m_ext_array = np.array([m_01, m_015, m_02])   #Un array con el array de magnitudes calculado para cada extinción


#Debo conseguir los datos de las estrellas a las mismas longitudes de onda que los datos de los filtros.
#Para ello, hago una interpolación.

interpol_ukk5 = interpolate.interp1d(ukk5_long, ukk5_normflux)

interpol_uka0 = interpolate.interp1d(uka0_long, uka0_normflux)

interpol_vega = interpolate.interp1d(vega_long, vega_absflux)

#Una vez hechas las interpolaciones, calculo los nuevos flujos para cada una de las longitudes de los filtros.

ukk5_normflux_interp_J = interpol_ukk5(J_long)
ukk5_normflux_interp_H = interpol_ukk5(H_long)
ukk5_normflux_interp_K = interpol_ukk5(K_long)

uka0_normflux_interp_J = interpol_uka0(J_long)
uka0_normflux_interp_H = interpol_uka0(H_long)
uka0_normflux_interp_K = interpol_uka0(K_long)

vega_absflux_interp_J = interpol_vega(J_long)
vega_absflux_interp_H = interpol_vega(H_long)
vega_absflux_interp_K = interpol_vega(K_long)

#Me piden también energía en fotones, calculo flujo en fotones 

vega_absflux_interp_J_ph = (vega_absflux_interp_J*J_long)/(h*c)
vega_absflux_interp_H_ph = (vega_absflux_interp_H*H_long)/(h*c)
vega_absflux_interp_K_ph = (vega_absflux_interp_K*K_long)/(h*c)

vega_absflux_wats_array = np.array([vega_absflux_interp_J,vega_absflux_interp_H,vega_absflux_interp_K])
vega_absflux_ph_array = np.array([vega_absflux_interp_J_ph,vega_absflux_interp_H_ph,vega_absflux_interp_K_ph])

#Flujo en fotones de UKK y UKA

ukk5_normflux_interp_J_ph = (ukk5_normflux_interp_J*J_long)/(h*c)
ukk5_normflux_interp_H_ph = (ukk5_normflux_interp_H*H_long)/(h*c)
ukk5_normflux_interp_K_ph = (ukk5_normflux_interp_K*K_long)/(h*c)

uka0_normflux_interp_J_ph = (uka0_normflux_interp_J*J_long)/(h*c)
uka0_normflux_interp_H_ph = (uka0_normflux_interp_H*H_long)/(h*c)
uka0_normflux_interp_K_ph = (uka0_normflux_interp_K*K_long)/(h*c)

print(len(J_long))
print(len(ukk5_normflux_interp_J))

#Calculo primero para VEGA (estrella de referencia)

def flux (mag, flux_d, trans_filter, dist_long, trans_sys):
    flux_mag_cero = sum(flux_d*trans_filter)*dist_long*trans_sys
    flux_tot_wats = 10**(-mag/2.5)*flux_mag_cero
    return flux_tot_wats

#Creo ahora unas listas donde meter mis resultados

flux_final_W = []
flux_final_ph = []
energy_W = []
energy_ph = []

#Calculo primero sin extinción
for i in range (len(trans_filts_array)):
  for m in m_array:
    flux_W = flux(m, vega_absflux_wats_array[i], trans_filts_array[i], dist_longs_array[i], trans_sys)
    flux_ph = flux(m, vega_absflux_ph_array[i], trans_filts_array[i], dist_longs_array[i], trans_sys)
    
    flux_final_W.append(flux_W)
    flux_final_ph.append(flux_ph)
    energy_W.append(flux_W*area)
    energy_ph.append(flux_ph*area)

#Devuelve array, con la primera triada para filtro J, [10,15,20], segunda triada para H [10,15,20], tercera triada para K [10,15,20]

#print(flux_final_W)
#print(flux_final_ph)
#print(energy_W)
#print(energy_ph)



#Calculo para las distintas extinciones

#Extinción 0.1

flux_final_W_01 = []
flux_final_ph_01 = []
energy_W_01 = []
energy_ph_01 = []

for i in range (len(trans_filts_array)):
  for m in m_01:
    flux_W = flux(m, vega_absflux_wats_array[i], trans_filts_array[i], dist_longs_array[i], trans_sys)
    flux_ph = flux(m, vega_absflux_ph_array[i], trans_filts_array[i], dist_longs_array[i], trans_sys)
    
    flux_final_W_01.append(flux_W)
    flux_final_ph_01.append(flux_ph)
    energy_W_01.append(flux_W*area)
    energy_ph_01.append(flux_ph*area)

#Devuelve array, con la primera triada para filtro J, [10,15,20], segunda triada para H [10,15,20], tercera triada para K [10,15,20]

#print(flux_final_W_01)
#print(flux_final_ph_01)
print(energy_W_01)
print(energy_ph_01)


#Extinción 0.15

flux_final_W_015 = []
flux_final_ph_015 = []
energy_W_015 = []
energy_ph_015 = []

for i in range (len(trans_filts_array)):
  for m in m_015:
    flux_W = flux(m, vega_absflux_wats_array[i], trans_filts_array[i], dist_longs_array[i], trans_sys)
    flux_ph = flux(m, vega_absflux_ph_array[i], trans_filts_array[i], dist_longs_array[i], trans_sys)
    
    flux_final_W_015.append(flux_W)
    flux_final_ph_015.append(flux_ph)
    energy_W_015.append(flux_W*area)
    energy_ph_015.append(flux_ph*area)

#Devuelve array, con la primera triada para filtro J, [10,15,20], segunda triada para H [10,15,20], tercera triada para K [10,15,20]

#print(flux_final_W_015)
#print(flux_final_ph_015)
print(energy_W_015)
print(energy_ph_015)


#Extinción 0.2

flux_final_W_02 = []
flux_final_ph_02 = []
energy_W_02 = []
energy_ph_02 = []

for i in range (len(trans_filts_array)):
  for m in m_02:
    flux_W = flux(m, vega_absflux_wats_array[i], trans_filts_array[i], dist_longs_array[i], trans_sys)
    flux_ph = flux(m, vega_absflux_ph_array[i], trans_filts_array[i], dist_longs_array[i], trans_sys)
    
    flux_final_W_02.append(flux_W)
    flux_final_ph_02.append(flux_ph)
    energy_W_02.append(flux_W*area)
    energy_ph_02.append(flux_ph*area)

#Devuelve array, con la primera triada para filtro J, [10,15,20], segunda triada para H [10,15,20], tercera triada para K [10,15,20]

#print(flux_final_W_02)
#print(flux_final_ph_02)
print(energy_W_02)
print(energy_ph_02)


#Con esto tengo los apartados a y b, queda el c, representar graficamente

#Voy a pintar el espectro de los objetos Pickles dados por profe y comparar con mis datos. Usaré magnitud 10 sin extinción.


flux_W_J = flux_final_W[0]
flux_W_H = flux_final_W[3]
flux_W_K = flux_final_W[6]

flux_ph_J = flux_final_ph[0]
flux_ph_H = flux_final_ph[3]
flux_ph_K = flux_final_ph[6]

"""

print(flux_W_J)
print(flux_W_H)
print(flux_W_K)

print(flux_ph_J)
print(flux_ph_H)
print(flux_ph_K)
"""

#Calculo los flujos absolutos de K5III para cada filtro a partir de sus flujos normalizados a 5556 Angstrom.

#Primero creo una imagen ara representar la transmisividad de los filtros

imag1 = plt.figure("Transmisión filtros en tanto por uno") 

plt.plot(J_long,J_trans, "b", label = "J")
plt.plot(H_long,H_trans, "g", label = "H")
plt.plot(K_long,K_trans, "r", label = "K")
plt.legend()
plt.xlabel("Longitud de onda $[\mu m]$")
plt.ylabel("Flujo (tanto por uno) ")
plt.title("Transmisión de filtros en tanto por uno")
plt.grid(linestyle="--", linewidth = 0.5)
imag1.savefig("Transmision_tanto_por_uno.png")


#Empiezo con ukk5 en watios

ukk5_integ_norm_J = 0.3*(sum(ukk5_normflux_interp_J*J_trans)*dist_J_long)
ukk5_C_J = flux_W_J/ukk5_integ_norm_J
ukk5_absflux_J_W = ukk5_C_J*ukk5_normflux_interp_J


ukk5_integ_norm_H = 0.3*(sum(ukk5_normflux_interp_H*H_trans)*dist_H_long)
ukk5_C_H = flux_W_H/ukk5_integ_norm_H
ukk5_absflux_H_W = ukk5_C_H*ukk5_normflux_interp_H


ukk5_integ_norm_K = 0.3*(sum(ukk5_normflux_interp_K*K_trans)*dist_K_long)
ukk5_C_K = flux_W_K/ukk5_integ_norm_K
ukk5_absflux_K_W = ukk5_C_K*ukk5_normflux_interp_K


#Calculo también estos flujos en photons

#Debo convertir el flujo normalizado a fotones, paso a fotones ukk5_normflux
ukk5_normflux_ph = (ukk5_normflux*ukk5_long)/(h*c)

#Esta normalizado al flujo de 5550 A, debo volver a normalizar en fotones como estaba normalizado ya en watts para luego pasar a flujo absoluto
ukk5_normflux_ph_renorm = ukk5_normflux_ph/((1.029639*5550)/(h*c))

ukk5_normflux_ph_renorm_J = ukk5_normflux_interp_J_ph/((1.029639*5550)/(h*c))
ukk5_normflux_ph_renorm_H = ukk5_normflux_interp_H_ph/((1.029639*5550)/(h*c))
ukk5_normflux_ph_renorm_K = ukk5_normflux_interp_K_ph/((1.029639*5550)/(h*c))

#Ya conozco todos los flujos normalizados en photons, calculo los flujos absolutos en photons


ukk5_integ_norm_J_ph = 0.3*(sum(ukk5_normflux_ph_renorm_J*J_trans)*dist_J_long)
ukk5_C_J_ph = flux_ph_J/ukk5_integ_norm_J_ph
ukk5_absflux_J_ph = ukk5_C_J_ph*ukk5_normflux_ph_renorm_J

ukk5_integ_norm_H_ph = 0.3*(sum(ukk5_normflux_ph_renorm_H*H_trans)*dist_H_long)
ukk5_C_H_ph = flux_ph_H/ukk5_integ_norm_H_ph
ukk5_absflux_H_ph = ukk5_C_H_ph*ukk5_normflux_ph_renorm_H

ukk5_integ_norm_K_ph = 0.3*(sum(ukk5_normflux_ph_renorm_K*K_trans)*dist_K_long)
ukk5_C_K_ph = flux_ph_K/ukk5_integ_norm_K_ph
ukk5_absflux_K_ph = ukk5_C_K_ph*ukk5_normflux_ph_renorm_K


#Ahora imprimo actuación por filtros de ukk5iii en watts y photons

imag2 = plt.figure("Transmisión UKK5III Watts")
plt.plot(J_long, ukk5_absflux_J_W, "b", label = "J")
plt.plot(H_long, ukk5_absflux_H_W, "g", label = "H")
plt.plot(K_long, ukk5_absflux_K_W, "r", label = "K")
plt.legend()
plt.xlabel("Longitud de onda $[\mu m]$")
plt.ylabel("Flujo $[W \, m^{-2} \, \mu m^{-1}]$ ")
plt.title("Transmisión de UKK5III [Watts]")
plt.grid(linestyle="--", linewidth = 0.5)
imag2.savefig("Transmision_ukk5iii_watts.png")

imag3 = plt.figure("Transmisión UKK5III Photons")
plt.plot(J_long, ukk5_absflux_J_ph, "b", label = "J")
plt.plot(H_long, ukk5_absflux_H_ph, "g", label = "H")
plt.plot(K_long, ukk5_absflux_K_ph, "r", label = "K")
plt.legend()
plt.xlabel("Longitud de onda $[\mu m]$")
plt.ylabel("Flujo $[ph/s \; m^{-2} \, \mu m^{-1}]$ ")
plt.title("Transmisión de UKK5III [photons]")
plt.grid(linestyle="--", linewidth = 0.5)
imag3.savefig("Transmision_ukk5iii_photons.png")



#Repito todo lo anterior para UKA0V


#Empiezo con uka0v en watios

uka0_integ_norm_J = 0.3*(sum(uka0_normflux_interp_J*J_trans)*dist_J_long)
uka0_C_J = flux_W_J/uka0_integ_norm_J
uka0_absflux_J_W = uka0_C_J*uka0_normflux_interp_J

uka0_integ_norm_H = 0.3*(sum(uka0_normflux_interp_H*H_trans)*dist_H_long)
uka0_C_H = flux_W_H/uka0_integ_norm_H
uka0_absflux_H_W = uka0_C_H*uka0_normflux_interp_H

uka0_integ_norm_K = 0.3*(sum(uka0_normflux_interp_K*K_trans)*dist_K_long)
uka0_C_K = flux_W_K/uka0_integ_norm_K
uka0_absflux_K_W = uka0_C_K*uka0_normflux_interp_K


#Calculo también estos flujos en photons

#Debo convertir el flujo normalizado a fotones, paso a fotones uka0_normflux
uka0_normflux_ph = (uka0_normflux*uka0_long)/(h*c)

#Esta normalizado al flujo de 5550 A, debo volver a normalizar en fotones como estaba normalizado ya en watts para luego pasar a flujo absoluto
uka0_normflux_ph_renorm = uka0_normflux_ph/((1.011062*5550)/(h*c))

uka0_normflux_ph_renorm_J = uka0_normflux_interp_J_ph/((1.011062*5550)/(h*c))
uka0_normflux_ph_renorm_H = uka0_normflux_interp_H_ph/((1.011062*5550)/(h*c))
uka0_normflux_ph_renorm_K = uka0_normflux_interp_K_ph/((1.011062*5550)/(h*c))

#Ya conozco todos los flujos normalizados en photons, calculo los flujos absolutos en photons


uka0_integ_norm_J_ph = 0.3*(sum(uka0_normflux_ph_renorm_J*J_trans)*dist_J_long)
uka0_C_J_ph = flux_ph_J/uka0_integ_norm_J_ph
uka0_absflux_J_ph = uka0_C_J_ph*uka0_normflux_ph_renorm_J

uka0_integ_norm_H_ph = 0.3*(sum(uka0_normflux_ph_renorm_H*H_trans)*dist_H_long)
uka0_C_H_ph = flux_ph_H/uka0_integ_norm_H_ph
uka0_absflux_H_ph = uka0_C_H_ph*uka0_normflux_ph_renorm_H

uka0_integ_norm_K_ph = 0.3*(sum(uka0_normflux_ph_renorm_K*K_trans)*dist_K_long)
uka0_C_K_ph = flux_ph_K/uka0_integ_norm_K_ph
uka0_absflux_K_ph = uka0_C_K_ph*uka0_normflux_ph_renorm_K


#Ahora imprimo actuación por filtros de uka0v en watts y photons

imag4 = plt.figure("Transmisión uka0III Watts")
plt.plot(J_long, uka0_absflux_J_W, "b", label = "J")
plt.plot(H_long, uka0_absflux_H_W, "g", label = "H")
plt.plot(K_long, uka0_absflux_K_W, "r", label = "K")
plt.legend()
plt.xlabel("Longitud de onda $[\mu m]$")
plt.ylabel("Flujo $[W \, m^{-2} \, \mu m^{-1}]$ ")
plt.title("Transmisión de UKA0V [Watts]")
plt.grid(linestyle="--", linewidth = 0.5)
imag4.savefig("Transmision_uka0v_watts.png")

imag5 = plt.figure("Transmisión uka0v Photons")
plt.plot(J_long, uka0_absflux_J_ph, "b", label = "J")
plt.plot(H_long, uka0_absflux_H_ph, "g", label = "H")
plt.plot(K_long, uka0_absflux_K_ph, "r", label = "K")
plt.legend()
plt.xlabel("Longitud de onda $[\mu m]$")
plt.ylabel("Flujo $[ph/s \; m^{-2} \, \mu m^{-1}]$ ")
plt.title("Transmisión de UKA0V [fotones]")
plt.grid(linestyle="--", linewidth = 0.5)
imag5.savefig("Transmision_uka0v_photons.png")



