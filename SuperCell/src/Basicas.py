import matplotlib.pyplot as plt
import numpy as np
from matplotlib import collections  as mc
import math
import os
plt.rcParams['figure.figsize'] = (10, 10)

'''-------------------------------------------------------------------------------------------------------------'''
#Funciones Básicas
def sumaV(x,y):
    '''Suma 2 vectores dados'''
    a,b=x
    c,d=y
    return (a+c,b+d)

def multV(n,a):
    '''Multiplica un vector a por una constante n'''
    a1,a2 = a
    return (a1*n,a2*n)

def m2V(n,m,s):
    '''Multiplica la matriz [n,m] por el vector s regresando la suma de n*s1+m*s2'''
    a,b=s
    return sumaV(multV(a,n),multV(b,m))

def rota(vect, theta):
    '''Rota un vector en un ángulo de Theta grados'''
    ang = float((theta/180.0)*math.pi)
    (x,y) = vect
    vr = ((math.cos(ang)*x)-(math.sin(ang)*y), (math.sin(ang)*x)+(math.cos(ang)*y))
    return vr

def dist(a, b):
    '''Calcula la distancia entre 2 puntos en el plano'''
    a1, a2 = a
    b1, b2 = b
    return math.sqrt(((b1-a1)**2)+((b2-a2)**2))

def long(v):
    '''Calcula la longitud de un vector'''
    return dist((0,0),v)

def cAng(u,v):
    '''Calcula el ángulo interno entre 2 vectores'''
    (u1,u2),(v1,v2) = u, v
    return math.degrees(math.acos(((u1*v1)+(u2*v2))/(dist((0,0),u)*dist((0,0),v))))

def cRot(v):
    '''Calcula el ángulo de un vector con respecto al eje x'''
    return cAng((1,0),v)

def getLim(u,v,m,n):
    '''Calcula los valores mínimo y máximo en "x" y "y" del rombo formado por m*u y n*v'''
    (a1,a2) = multV(m,u)
    (b1,b2) = multV(n,v)
    (c1,c2) = m2V(u,v,(m,n))
    xma = max(a1,b1,c1,0)
    xmi = min(a1,b1,c1,0)
    yma = max(a2,b2,c2,0)
    ymi = min(a2,b2,c2,0)
    return [xmi-1, xma+1], [ymi-1, yma+1]
