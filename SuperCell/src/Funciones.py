from .Lattice import *

# Metodos derivados

def isitin(r,cent, sr, slvl):
    '''
    Verifica que átomos pertenecientes a una celda de la Red "r" con centro en "cen" se encuentran en la celda principal de la Red "sr".
    Todos los átomos que lo estén se agregan a los átomos de sr
    '''
    er = 1/(10**7) #Error de calculo aceptable en la frontera de la Celda unitaria de sr
    (u1,u2) = r.a
    (v1,v2) = r.b
    (p1,p2) = sr.a
    (q1,q2) = sr.b
    eq0 = (p2*q1)-(p1*q2)
    eq1 = (q1*u2)-(q2*u1)
    eq2 = (q1*v2)-(q2*v1)
    eq3 = (p2*u1)-(p1*u2)
    eq4 = (p2*v1)-(p1*v2)
    for c in r.atms:# Iteramos en cada lista de Atomos de la Red r
        nc=[]
        for a in c:# Iteramos en cada atomo de la lista c
            (x,y) = sumaV(cent,a.pos)
            # Calculamos la posición del átomo con respecto a las coordenadas expresadas en los vectores de sr
            nx = (eq1*x+eq2*y)/eq0 
            ny = (eq3*x+eq4*y)/eq0
            # Evaluamos si se encuentra dentro de la Celda minima de sr
            if (nx<(1+er) and nx>(0-er)) and (ny<(1+er) and ny>(0-er)):
                aPosZ = ((a.posZ*r.detachment)+slvl)/sr.detachment
                nAtm = Atomo((nx,ny),posZ = aPosZ ,color=a.color,sig=a.sig)
                #print("Agregado atomo",nAtm)
                nAtm.clasifica(sr.atms)
    for e in r.enls:
        (x1,y1) = sumaV(cent,e[0])
        (x2,y2) = sumaV(cent,e[1])
        (ox,oy) = ((eq1*x1+eq2*y1)/eq0,(eq3*x1+eq4*y1)/eq0)
        f = ((eq1*x2+eq2*y2)/eq0,(eq3*x2+eq4*y2)/eq0)
        if (ox<1 and ox>0) and (oy<1 and oy>0):
            sr.enls.append([(ox,oy),f])
    return 1

def megeCut(mo, sm, lvl=0):
    '''
    Método auxiliar para el superMesh.
    '''
    (u1,u2), (v1,v2) = mo.a, mo.b
    (p1,p2), (q1,q2) = sm.a, sm.b
    np1 = ((p2*v1)-(p1*v2))/((u2*v1)-(u1*v2))
    np2 = ((p1*u2)-(p2*u1))/((u2*v1)-(u1*v2))
    nq1 = ((q2*v1)-(q1*v2))/((u2*v1)-(u1*v2))
    nq2 = ((q1*u2)-(q2*u1))/((u2*v1)-(u1*v2))
    npq1 = np1+nq1
    npq2 = np2+nq2
    lu = [round(min(np1,nq1,npq1,0)-1),round(max(np1,nq1,npq1,0)+1)]
    lv = [round(min(np2,nq2,npq2,0)-1),round(max(np2,nq2,npq2,0)+1)]
    for i in range(lu[1]-lu[0]):
        a = i+lu[0]
        for j in range(lv[1]-lv[0]):
            b = j+lv[0]
            isitin(mo,(a,b),sm,lvl)
    return 1

def superMesh(sa,sb,loLs):
    '''Crea una SuperMalla en bace a una lista de Redes "loLs" con los vectores sa y sb'''
    sR = Red(sa,sb)
    sR.enls = []
    detachment = 0
    for l in loLs:
        detachment = detachment + l.detachment
    sR.detachment = detachment
    sR.prof=0
    newName="SuperLattice"
    i=0
    for m in loLs:
        megeCut(m, sR, lvl=i)
        sR.prof=sR.prof+m.prof
        newName=newName+" ["+m.name+"]"
        i = i + m.detachment
    sR.name = newName
    return sR

def hexa6(p,atms=['C','C'],name=''):
    '''Crea una Red exagonal s6 con constante de red P'''
    u,v=(p,0.0),(-p/2,math.sqrt(3)*(p/2))
    p1,p2,p3,p4 = (1/3,2/3),(2/3,1/3),(1/3,-1/3),(4/3,2/3)
    ats = [Atomo(p1, sig = atms[0]),Atomo(p2, sig = atms[1])]
    return Red(u,v,atms=ats,name=name,enls=[(p1,p2),(p2,p3),(p2,p4)])

def hexa3(p,atms=['C','C'],name=''):
    '''Crea una Red exagonal s3 con constante de red P'''
    u,v=(p,0.0),(-p/2,math.sqrt(3)*(p/2))
    p1,p2,p3,p4 = (0.0,0.0),(1/3,2/3),(0,1),(1,1)
    ats = [Atomo(p1, sig = atms[0]),Atomo(p2, sig = atms[1])]
    return Red(u,v,atms=ats,name=name,enls=[(p1,p2),(p2,p3),(p2,p4)])

def rectMesh(p1,p2,atms='C',name=''):
    '''Crea una Red cuadrada con constantes de red p1 y p2'''
    u,v = (p1,0.0),(0.0,p2)
    p1,p2,p3 = (1/2,1/2),(3/2,1/2),(1/2,3/2)
    ats = [Atomo(p1,sig = atms)]
    return Red(u,v,atms=ats,name=name,enls=[(p1,p2),(p1,p3)])
    
def grafeno():
    '''Crea una red de Grafeno con ambos átomos en el centro'''
    return hexa6(2.44, name='Grafeno')
    
def grafeno3():
    '''Crea una red de Grafeno con un átomo en un extremo'''
    return hexa3(2.44, name='Grafeno(s3)')

def blackPhospho():
    m1=rectMesh(3.3061099052,4.552418232)
    m1.name='Black-Phosphorene'
    p1,p2=(0.000000000,0.913483083),(0.500000000,0.579813302)
    p3,p4=(0.000000000,0.079836130),(0.500000000,0.413437814)
    ats = [Atomo(p1,sig='P',posZ=0.266835123),Atomo(p2,sig='P',posZ=0.266945183),Atomo(p3,sig='P',posZ=0.181006327),Atomo(p4,sig='P',posZ=0.181094214)]
    m1.atms[0] = ats
    #m1.showNM(1,1)
    return m1

def limpia(loa):
    '''Quita atomos repetidos de una lista de atomos'''
    err=1/(10**8)
    res = []
    for i in range(len(loa)):
        pasa=True
        x1,x2 = loa[i].pos
        if abs(x1-1) < err:
            for j in range(len(loa)):
                y1,y2=loa[j].pos
                if (abs(y1) < err) and (abs(x2-y2) < err):
                    pasa = pasa and False
            if pasa:
                res.append(loa[i])
        elif abs(x2-1) < err:
            pasa=True
            for j in range(i+1,len(loa)):
                y1,y2=loa[j].pos
                if abs(x1-y1) < err:
                    pasa = pasa and False
        if pasa:
            res.append(loa[i])
    return res

def transfVs(u,v,t):
    '''Transforma los vectores u y v al multiplicar la matriz [[u1,v1],[u2,v2]] por la matriz [[t1,t2],[t3,t4]]'''
    m,n,p,q = t
    return m2V(u,v,(m,p)), m2V(u,v,(n,q))

def buscaSVect(vectU,vectV, th, rango=15, limDelta=0.1, show=True):
    lim = limDelta
    f1, f2 = 0, 0
    res = [[],[]]
    rmin = [0,0,0,0,0.0]
    rmin2 = [0,0,0,0,0.0]
    ang = math.radians(th)
    cos = math.cos(ang)
    sen = math.sin(ang)
    ru, rv = rota(vectU,th), rota(vectV,th)
    (u1,u2) = vectU
    (v1,v2) = vectV
    ax1 = (u2*v1)-(u1*v2)
    ax2 = (u1*v1)+(u2*v2)
    ax3 = (v1**2)+(v2**2)
    ax4 = (u1**2)+(u2**2)
    delta=0.0
    for k in range(1,(2*rango)+1):
        for i in range(k+1):
            j = k-i
            if(i<(rango+1) and j<(rango+1)):
                # Buscando en a+
                a,b = i,-j
                c = (a*(ax1*cos-ax2*sen)/ax1)-(b*(ax3*sen)/ax1)
                d = (b*(ax1*cos+ax2*sen)/ax1)+(a*(ax4*sen)/ax1)
                r1 = sumaV(multV(a,vectU),multV(b,vectV))
                r2 = sumaV(multV(round(c),ru),multV(round(d),rv))
                delta = dist((0,0),r1)/dist((0,0),r2)
                err = dist(r1,r2)*(abs(delta))
                if (err<limDelta):
                    if(abs(1-delta)<0.03):
                        res[0].append([[a,b],[round(c),round(d)],delta])
                        print(">({},{})-({},{}): Delta={}%".format(a,b,round(c),round(d),delta*100),":",dist(r1,r2))
                # Buscando en a-
                if j!=0:
                    a,b = i,j
                    c = (a*(ax1*cos-ax2*sen)/ax1)-(b*(ax3*sen)/ax1)
                    d = (b*(ax1*cos+ax2*sen)/ax1)+(a*(ax4*sen)/ax1)
                    r1 = sumaV(multV(a,vectU),multV(b,vectV))
                    r2 = sumaV(multV(round(c),ru),multV(round(d),rv))
                    delta = dist((0,0),r1)/dist((0,0),r2)
                    err = dist(r1,r2)*(abs(delta))
                    if (err<limDelta):
                        if(abs(1-delta)<0.03):
                            res[1].append([[a,b],[round(c),round(d)],delta])
                            print(">>({},{})-({},{}): Delta={}%".format(a,b,round(c),round(d),delta*100,":",dist(r1,r2)))
    return res