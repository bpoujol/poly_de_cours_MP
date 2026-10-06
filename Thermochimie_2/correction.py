# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

### Constantes physiques
R = 8.314
#Loi d'Arrhénius
A = 1.19e10  # en L/mol/s
Ea = 1.33e5 # en J/mol
Delta_r_H = -226000. #J/mol
#Coefficients de Laplace
gamma_NO = 1.4
gamma_CO = 1.4
gamma_NO2 = 1.3
gamma_CO2 = 1.3
gamma_air = 1.4


### Caracteristiques de l'etat initial
#Volume
V = 4.   # en L
#Quantite de matiere en reactifs
n0 = 0.05 # en mol
#Quantite de matiere en air
n_air = 4 #en mol
#Temperature initiale
T0 = 700 #K

#Parametres de la simulation
#Choix du pas de temps
dt = .1 #s
#Duree de simulation
t_max = 100 #s
#Liste des temps
les_t = np.arange(0,t_max,dt)
#Nombre de pas de temps
N = len(les_t)

#initialisation des températures et avancements
les_T = np.zeros((N))
les_xi = np.zeros((N))

#Calcul de la capacité thermique totale en fonction de l'avancement
def Cv(xi):
  return (n0-xi)*R/(gamma_NO2-1) + (n0-xi)*R/(gamma_CO-1) + xi*R/(gamma_CO2-1) + xi*R/(gamma_NO-1) + n_air*R/(gamma_air-1)

#Vitesse de réaction en fonction de la température et de xi :
def v(T,xi):
  #Loi d'Arrhénius
  k = A*np.exp(-Ea/R/T)
  #Loi de vitesse
  v = k*((n0-xi)/V)**2
  return v

#Etat initial
#Température
T = T0
#Avancement
xi = 0

#Boucle sur tous les temps de la simulation
for i in range(N):
  #On enregistre la température
  les_T[i] = T
  les_xi[i] = xi

  #METHODE D'EULER
  #Variation infinitesimale de xi et T
  d_xi = v(T,xi)*dt*V
  d_T = -Delta_r_H*d_xi/Cv(xi)
  #On incrémente xi et T
  xi = xi+d_xi
  T = T+d_T

  #On vérifie qu'on n'a pas dépassé l'avancement maximal
  if xi>n0:
    xi = n0



# Affichage des résultats
plt.figure(figsize=(8,12))

plt.subplot(2,1,1)
plt.title('Température',fontsize=16)
plt.plot(les_t,les_T-273)
plt.xlabel('Temps [s]',fontsize=16)
plt.ylabel('Température [°C]',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlim(np.min(les_t),np.max(les_t)) #Limites des axes en x
plt.ylim(0,1000) # Limites des axes en y

plt.subplot(2,1,2)
plt.title('Avancement',fontsize=16)
plt.plot(les_t,les_xi)
plt.xlabel('Temps [s]',fontsize=16)
plt.ylabel('Avancement [mol]',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlim(np.min(les_t),np.max(les_t)) #Limites des axes en x
plt.ylim(0,1.2*n0) # Limites des axes en y

plt.tight_layout()

#Affichage
plt.show()
