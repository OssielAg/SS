import matplotlib.pyplot as plt
import numpy as np
from matplotlib import collections  as mc
import math
import random
import re
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
    x = round(((u1*v1)+(u2*v2))/(dist((0,0),u)*dist((0,0),v)),5)
    return math.degrees(math.acos(x))

def cRot(v):
    '''Calcula el ángulo de un vector con respecto al eje x'''
    return cAng((1,0),v)

def vTm(v1,v2):
    '''
    Transforma dos vectores de dimencion 2 a una matríz de 2x2
    '''
    (a1,a2), (b1,b2) = v1, v2
    return [[a1,b1],[a2,b2]]

def mTv(m):
    '''
    Transforma una Matríz de 2x2 a un par de vectores de 2 dimenciones
    '''
    [[a1,b1],[a2,b2]] = m
    return (a1,a2), (b1,b2)
    
def det(t):
    '''
    Calcula el determinante de una matriz t de 2x2
    '''
    return (t[0][0]*t[1][1])-(t[0][1]*t[1][0])

def inv2x2(m):
    '''
    Calcula la matriz inversa de una matriz de 2x2
    '''
    [[a,b],[c,d]] = m
    de = det(m)
    return [[d/de, -b/de], [-c/de, a/de]]

def m2M(m1,m2):
    '''
    Regresa la Matriz resultante de multiplica las 2 matrices de 2x2 m1 y m2
    '''
    [[a,b],[c,d]] = m1
    [[x,y],[z,w]] = m2
    return [[(a*x+b*z),(a*y+b*w)],[(c*x+d*z),(c*y+d*w)]]

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

def pC(a,b):
    '''
    Calcula el producto Crux de los vectores a y b en R3
    '''
    a1,a2,a3 = a
    b1,b2,b3 = b
    return ((a2*b3-a3*b2),(a3*b1-a1*b3),(a1*b2-a2*b1))

def pP(a,b):
    '''
    Calcula el producto punto de los vectores a y b en R3
    '''
    a1,a2,a3 = a
    b1,b2,b3 = b
    return a1*b1+a2*b2+a3*b3

def to2D(v):
    '''
    Recibe un vector en R3 y regresa su proyección en R2
    '''
    x,y,z=v
    return (x,y)

def cRecip(a,b,c):
    '''
    Calcula los vectores de la Red en el espacio reciproco de la red formada por los vectores a,b,c en R3
    '''
    (a1,a2,a3) = a
    (b1,b2,b3) = b
    (c1,c2,c3) = c
    e = math.pi/(pP(a,pC(b,c)))
    (x,y,z) = pC(b,c)
    u = (x*e, y*e, z*e)
    (x,y,z) = pC(a,c)
    v = (x*e, y*e, z*e)
    (x,y,z) = pC(a,b)
    w = (x*e, y*e, z*e)
    return [u,v,w]

def esRotacion(a,b,c,d,eps=0.001):
    '''
    Determina si el par de vectores a,b es una rotación del par de vectores c,d
    '''
    if abs(cAng(a,c)-cAng(b,d))<eps:
        return True
    if abs(cAng(a,d)-cAng(b,c))<eps:
        return True
    return False
    
def transforma2v(u,v,M):
    '''
    Obtiene los vectores ut y vt a partir de transformar u y v con una Matríz M 
    '''
    mv = vTm(u,v)
    m = m2M(mv,M)
    ut, vt = mTv(m)
    return ut, vt