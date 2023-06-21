import matplotlib.pyplot as plt
import numpy as np
from matplotlib import collections  as mc
from matplotlib.patches import Polygon
import copy
import math
import random
import re
import os
plt.rcParams['figure.figsize'] = (10, 10)

'''-------------------------------------------------------------------------------------------------------------'''
#Funciones Básicas

#--------------------Funciones para Vectores--------------------#
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
    '''Suma el el vector n multiplicado con el primer elemento del par s con el vector m multiplicado por el segundo elemento del par s'''
    a,b = s
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
    '''Calcula el ángulo formado entre los vectores 'u' y 'v'.'''
    (u1,u2),(v1,v2) = u, v
    x = round(((u1*v1)+(u2*v2))/(dist((0,0),u)*dist((0,0),v)),5)
    th = math.degrees(math.acos(x))
    if (u1*v2-u2*v1)<0: th=-th
    return th

def cRot(v):
    '''Calcula el ángulo de un vector con respecto al eje x'''
    return cAng((1,0),v)

def pC(a,b):
    '''
    Calcula el producto Cruz de los vectores a y b en R3
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
    
#--------------------Funciones para matrices--------------------#
def sumaM(A,B):
    '''
    Suma las matrices A y B
    '''
    C = np.array(A)+np.array(B)
    return C.tolist()

def multM(M,x):
    '''
    Multiplica la matriz M por una constante x
    '''
    xM=x*np.array(M)
    return xM.tolist()

    
def vTm(v1,v2):
    '''
    Transforma dos vectores de dimensión 2 a una matriz de 2x2
    '''
    (a1,a2), (b1,b2) = v1, v2
    return [[a1,b1],[a2,b2]]

def mTv(m):
    '''
    Transforma una Matriz de 2x2 a un par de vectores de 2 dimensiones
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

def transfVs(u,v,t):
    '''
    Transforma los vectores u y v al multiplicar la matriz [[u1,v1],[u2,v2]] por la matriz [[t1,t2],[t3,t4]]
    '''
    m,n,p,q = t
    return m2V(u,v,(m,p)), m2V(u,v,(n,q))


#--------------------Funciones misceláneas--------------------#
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
    Obtiene los vectores ut y vt a partir de transformar u y v con una Matriz M 
    '''
    mv = vTm(u,v)
    m = m2M(mv,M)
    ut, vt = mTv(m)
    return ut, vt

def acomoda(x, valor, es, tamMax):
    '''
    Acomoda un objeto 'x' en la lista ordenada 'es' de acuerdo a su 'valor' respetando que
    el tamaño de 'es' no sea mayor a tamMax.
    x     -> Elemento que queremos ingresar a la lista
    valor -> Valor dado al elemento 'x'
    es    -> Lista en que trabajamos
    tamMax-> Tamaño máximo aceptable para 'es'
    '''
    if len(es)<tamMax:
        es.append([x,valor])
        return es
    for i in range(len(es)):
        if valor<es[i][1]:
            es.insert(i,[x,valor])
            if len(es)>tamMax:
                trash=es.pop()
            return es
    return es

#--------------------Funciones para el espacio reciproco--------------------
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

def buscaEsquina(m, j, e1, e2):
    '''
    **Función auxiliar para calVerticesFBZ
    Indica con qué 'fila' de la matriz en que se opera se va a intercambiar.
    m  -> Matriz en que operamos
    j  -> Columna que estamos iterando
    e1 -> Primera dimensión de la base actual
    e2 -> Segunda dimensión de la base actual
    '''
    ind = -1
    val = (10.0)**300
    for i in range(8):
        if not(i==e1 or i==e2):
            if (m[i][j]>0):
                c = m[i][10]/m[i][j]
                if c < val:
                    val = c
                    ind = i
    return ind

def opera(m,e,s):
    '''
    **Función auxiliar para calcVerticesFBZ
    Opera en la matriz para dejar en la columna 'e' tengamos ceros excepto en m[s][e] donde tenemos 1
    m -> Matriz que operamos
    e -> Columna de la celda pivote
    s -> Fila de la celda pivote
    '''
    m[s] = m[s]/(m[s][e])
    for i in range(8):
        if i!=s:
            piv=m[i][e]
            m[i] = m[i]-(m[s]*piv)

def pmat(m):
    '''
    Imprime en pantalla la matriz m.
    '''
    for i in range(len(m)):
        lin = ""
        for j in range(len(m[0])):
            lin = lin+"\t"+str(round(m[i][j],2))
        print(lin)
    print("")
    
def dameVecinos(red):
    '''
    Regresa los 8 puntos de red del espacio reciproco de 'red' más cercanos al origen 
    '''
    b1, b2 = to2D(red.reciprocalVectors[0]), to2D(red.reciprocalVectors[1])
    d1 = long(m2V(b1,b2,(1,1)))
    d2 = long(m2V(b1,b2,(-1,1)))
    if d1<d2:
        n = round(d2/d1)
    else:
        n = round(d1/d2)
    vec = []
    for i in range(-(2*n),(2*n)+1):
        for j in range(-(2*n),(2*n)+1):
            if (i,j)!=(0,0):
                v = m2V(b1,b2,(i,j))
                vec.append([v,long(v),(i,j)])
    vec=sorted(vec, key=lambda op : op[1])
    vecinos=[]
    vC=[]
    ind=0
    while len(vecinos)<8 and ind<len(vec):
        if califica(vec[ind][2],vC):
            #vecinos.append([vec[ind][0],vec[ind][1],vec[ind][2]])
            vecinos.append(vec[ind][0])
            vC.append(vec[ind][2])
        ind+=1
    return vecinos

def califica(pos, lista):
    '''
    Califica de forma booleana si una posición 'pos' dada es vista desde el origen sin ser tapada por las posiciones
    guardadas en la lista dada.
    '''
    for i in range(len(lista)):
        (px,py) = pos
        (lx,ly) = lista[i]
        if px!=0:
            k=lx/px
            if k>0:
                if ly==k*py:
                    #print("({},{}) es tapado por ({},{})".format(px,py,lx,ly))
                    return False
        elif py!=0:
            k=ly/py
            if k>0:
                if lx==k*px:
                    #print("({},{}) es {} veces ({},{})".format(px,py,k,lx,ly))
                    return False
    return True
    
def calcVerticesFBZ(l):
    '''
    Calcula la posición de los vértices frontera de la FBZ a partir de los vectores recíprocos a* y b*
    proyectados en el plano XY.
    b1 -> Vector reciproco a* de la red
    b2 -> Vector reciproco b* de la red
    '''
    #b1,b2 = to2D(l.reciprocalVectors[0]), to2D(l.reciprocalVectors[1])
    pts = []
    vecinos=dameVecinos(l)
    '''[m2V(b1,b2,( 1, 0)),m2V(b1,b2,( 1, 1)),
             m2V(b1,b2,( 0, 1)),m2V(b1,b2,(-1, 1)),
             m2V(b1,b2,(-1, 0)),m2V(b1,b2,(-1,-1)),
             m2V(b1,b2,( 0,-1)),m2V(b1,b2,( 1,-1))]'''
    # Creamos el espacio para la matriz en que operaremos
    v = np.array(vecinos)
    eq = np.zeros((8,11))
    et = ["V1","V2","V3","V4","V5","V6","V7","V8"]
    xy = [10,10]
    cruce = [-1,-1]
    # Damos los valores correspondientes
    for i in range(8):
        eq[i][i] = 1
        eq[i][8] = 2*v[i][0]
        eq[i][9] = 2*v[i][1]
        eq[i][10] = (v[i][0]**2) + (v[i][1]**2)
    #pmat(eq) #---Imprime la matriz para hacer seguimiento
    # Entra valor de X
    ind = buscaEsquina(eq,8,xy[0],xy[1])
    xy[0]=ind
    cruce[0]=ind
    opera(eq,8,ind)
    # Entra valor de Y
    ind = buscaEsquina(eq,9,xy[0],xy[1])
    xy[1]=ind
    cruce[1]=ind
    opera(eq,9,ind)
    p = (round(eq[xy[0]][10],10),round(eq[xy[1]][10],10))
    pts.append(p)
    #pmat(eq) #---Imprime la matriz para hacer seguimiento
    cont=0
    while True:
        cont+=1
        e = cruce[0]
        s = buscaEsquina(eq,e,xy[0],xy[1])
        opera(eq,e,s)
        #pmat(eq) #---Imprime la matriz para hacer seguimiento
        cruce[0] = cruce[1]
        cruce[1] = s
        p2 = (round(eq[xy[0]][10],10),round(eq[xy[1]][10],10))
        if p2!=p:
            if p2==pts[0]:
                break
            p=p2
            pts.append(p)
        if cont>10:
            break
    
    return pts, eq

def reciprocalBackgroundMesh(l,vl,t):
    '''
    Calcula los puntos de Red de la red reciproca y la maya formada por las FBZ de estos.
    Regresa 2 listas con las posiciones X y Y de los puntos y una 'colección' de lineas usadas por
    'Polygon' para cargar un polígono en pantalla usando la función 'add_patch'.
    l  -> Red para la que se lleva a cabo la operación
    vl -> Vértices de la FBZ asociada a esta red
    t  -> Tamaño con el que se pintarán las lineas.
    '''
    xs = []
    ys = []
    enls=[]
    sra, srb = to2D(l.reciprocalVectors[0]), to2D(l.reciprocalVectors[1])
    for i in range(-10,10):
        for j in range(-10,10):
            (px,py) = m2V(sra,srb,(i,j))
            #Calcula la posición de cada centro
            xs.append(px)
            ys.append(py)
            o = sumaV((px,py),vl[len(vl)-1])
            for p in vl:#calcula las aristas de la FBZ de cada centro
                f = sumaV((px,py), p)
                enls.append([o, f])
                o = f
    xs = np.array(xs)
    ys = np.array(ys)
    linkList = mc.LineCollection(np.array(enls), colors='silver', linewidths=(t/10))
    return xs, ys, linkList