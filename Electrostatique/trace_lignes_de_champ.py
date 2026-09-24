#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Auteur : Basile Poujol
Créé le : 24 septembre 2025
Dernière modification : 23 septembre 2026
Code adapte a partir du code de Etienne Thibierge disponible sur :
https://www.etienne-thibierge.fr/cours_2025/13_potentiel_simu-HCl.py
"""

import numpy as np
import matplotlib.pyplot as plt

### Constantes physiques
eps0 = 8.85e-12  # en F.m-1

### Caracteristiques de la molecule
a = 127.4   # en pm
q = 0.2*1.6e-19 # en coulomb

#Coordonnées de l'atome de chlore
x_Cl = 0
y_Cl = -a/2
#Cordonnées de l'atome d'hydrogène
x_H = 0
y_H = a/2

### Caracteristiques de la grille
h = 450     # hauteur de la grille en pm
Nx = 100    # taille de la grille
Ny = 100
dx = h/Nx   #espace entre deux points de grille
dy = h/Ny


### Construction de la grille symetrique
xx = np.linspace(-h/2 , h/2, num=Nx+1) # liste de coordonnees x
yy = np.linspace(-h/2 , h/2, num=Ny+1) # liste de coordonnees y
x, y = np.meshgrid(xx, yy) #x et y sont deux matrices qui contiennent les coordonnees x et y sur la grille


### Calcul du potentiel :

#Remarque pour la suite : On peut faire des opérations directement avec les matrices numpy en python

# Distance à l'atome d'hydrogène en mètres
r_H = np.sqrt((x-x_H)**2+(y-y_H)**2)*1e-12
#Potentiel créé par l'atome d'hydrogène
V_H = q/(4*np.pi*eps0*r_H)

# Distance à l'atome de chlore en mètres
r_Cl = np.sqrt((x-x_Cl)**2+(y-y_Cl)**2)*1e-12
#Potentiel créé par l'atome de chlore
V_Cl = -q/(4*np.pi*eps0*r_Cl)

#Potentiel total : on applique le principe de superposition
V = V_Cl + V_H


### Calcul du champ
Ey, Ex = np.gradient(-V)
### np.gradient renvoie une **liste** dont les elements sont des np.array,
### On ne peut pas mettre de signe - devant une liste
### le signe - est donc a mettre a l'interieur de la fonction


### Traces :

# Ce qui est important ...
plt.figure(figsize=(12,12))

#On choisit de dessiner 500 equipotentielles (peut etre ajuste)
plt.contour(x,y,V, 500)
#Et on trace les lignes de champ avec la fonction streamplot
plt.streamplot(x,y,Ex,Ey, color='k', linewidth=.5,broken_streamlines=False,density=0.5)
#La fonction broken_streamlines=False permet de s'assurer qu'on a des lignes continues.

#On place les atomes :
#    de chlore
plt.plot([0],[-a/2], "o", markersize = 20, markeredgecolor='k',color='green')
plt.text(0,-a/2,'Cl', ha='center', va='center')
#    d'hydrogène
plt.plot([0],[a/2], "o", markersize = 20, markeredgecolor='k',color='white')
plt.text(0,a/2,'H', ha='center', va='center')

#Pour décorer
plt.axis('scaled') #Axes à l'échelle
plt.xlim(np.min(xx),np.max(xx)) #Limites des axes en x
plt.ylim(np.min(yy),np.max(yy)) # Limites des axes en y
plt.xlabel('x [pm]',fontsize=14)
plt.ylabel('y[pm]',fontsize=14)

#Affichage
plt.show()
