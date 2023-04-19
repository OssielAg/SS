from .Lattice import *
from .analyzer import *

# Métodos derivados

def isitin(r,cent, sr, slvl):
    '''
    Verifica que átomos pertenecientes a una celda de la Red "r" con centro en "cen" se encuentran en la celda principal de la Red "sr".
    Todos los átomos que lo estén se agregan a los átomos de sr
    '''
    er = 1/(10**3) #Error de calculo aceptable en la frontera de la Celda unitaria de sr
    (u1,u2) = r.a
    (v1,v2) = r.b
    (p1,p2) = sr.a
    (q1,q2) = sr.b
    eq0 = (p2*q1)-(p1*q2)
    eq1 = (q1*u2)-(q2*u1)
    eq2 = (q1*v2)-(q2*v1)
    eq3 = (p2*u1)-(p1*u2)
    eq4 = (p2*v1)-(p1*v2)
    for c in r.atms:# Iteramos en cada lista de Átomos de la Red r
        nc=[]
        for a in c:# Iteramos en cada átomo de la lista c
            (x,y) = sumaV(cent,a.pos)
            # Calculamos la posición del átomo con respecto a las coordenadas expresadas en los vectores de sr
            nx = (eq1*x+eq2*y)/eq0 
            ny = (eq3*x+eq4*y)/eq0
            # Evaluamos si se encuentra dentro de la Celda mínima de sr
            if (nx<(1+er) and nx>(0-er)) and (ny<(1+er) and ny>(0-er)):
                aPosZ = ((a.posZ*r.detachment)+slvl)/sr.detachment
                nAtm = Atomo((nx,ny),posZ = aPosZ ,color=a.color,sig=a.sig)
                #print("Agregado átomo",nAtm)
                nAtm.clasifica(sr.atms)
    for e in r.enls:
        (x1,y1) = sumaV(cent,e[0])
        (x2,y2) = sumaV(cent,e[1])
        (ox,oy) = ((eq1*x1+eq2*y1)/eq0,(eq3*x1+eq4*y1)/eq0)
        (fx,fy) = ((eq1*x2+eq2*y2)/eq0,(eq3*x2+eq4*y2)/eq0)
        if (ox<1+er and ox>0-er) and (oy<1+er and oy>0-er):
            sr.enls.append([(ox,oy),(fx,fy)])
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

def superMesh(sa,sb,layerList):
    '''
    Crea una Súper Red en base a una lista de Redes "layerList" con los vectores primitivos sa y sb
    sa        -> Vector primitivo 'a' de la Súper Red 
    sb        -> Vector primitivo 'b' de la Súper Red
    layerList -> Lista de capas que forman el sistema al que identifica la Súper Red
    '''
    sR = Red(sa,sb)
    sR.enls = []
    detachment = 0
    for l in layerList:
        detachment = detachment + l.detachment
    sR.detachment = detachment
    sR.prof=0
    newName="SuperLattice"
    i=0
    for m in layerList:
        megeCut(m, sR, lvl=i)
        sR.prof=sR.prof+m.prof
        newName=newName+" ["+m.name+"]"
        i = i + m.detachment
    sR.name = newName
    sR.layerList = layerList
    return sR

def limpia(loa, acc=8):
    '''
    Quita átomos repetidos de una lista de átomos
    loa -> Lista de átomos que limpiaremos
    acc -> Exactitud con la que trabajaremos
    '''
    err=1/(10**acc)
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
    '''
    Transforma los vectores u y v al multiplicar la matriz [[u1,v1],[u2,v2]] por la matriz [[t1,t2],[t3,t4]]
    '''
    m,n,p,q = t
    return m2V(u,v,(m,p)), m2V(u,v,(n,q))

def buscaSVects(vectU,vectV, th, rango=15, limDelta=0.1, show=True):
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

def calculaPares(r1, r2, th = 0.0, maxIt=15, eps=0.1, show=False):
    '''
    Calcula los pares enteros (a,b) y (c,d) tales que si u,v son los vectores generadores de 'r1' y rp,rq los
    vectores generadores de 'r2' rotada en 'th' grados, entonces P1 = (au + bv) y P2 = (c(rp) + d(rq)) difieren en
    menos de 'eps', además a y b son a lo más 'maxIt'
    
    Regresa una lista doble con todos los resultados que cumplen lo anterior separados en los que tienen b positiva
    y b negativa (hace la búsqueda en los cuadrantes I y IV del plano cartesiano), además regresa también el promedio
    de los errores mínimos en ambas zonas.
    '''
    (u,v), (p,q) = r1.getVectors(), r2.getVectors()
    (u_1,u_2), (v_1,v_2) = u, v
    (p_1,p_2), (q_1,q_2) = rp, rq = rota(p,th), rota(q,th)
    if th==0.0:
        th=r2.theta
    res = []
    f = 1/2 # Factor de importancia de el tamaño del vector resultante
    rango = maxIt
    delta = 0.0
    minE1 = 100
    minE2 = 100
    eq0 = (p_2*q_1)-(p_1*q_2)
    eq1 = (q_1*u_2)-(q_2*u_1)
    eq2 = (q_1*v_2)-(q_2*v_1)
    eq3 = (p_2*u_1)-(p_1*u_2)
    eq4 = (p_2*v_1)-(p_1*v_2)
    for k in range(1,(2*rango)+1):
        for i in range(k+1):
            j = k-i
            if(i<(rango+1) and j<(rango+1)):
                # Buscando en b+
                a,b = i,-j
                c = ((eq1*a)+(eq2*b))/(eq0)
                d = ((eq3*a)+(eq4*b))/(eq0)
                # Vector esperado
                r1 = m2V(u,v,(a,b))
                #Vector aproximado
                r2 = m2V(rp,rq,(round(c),round(d)))
                # Calcula la proporción entre el tamaño del vector esperado y el aproximado
                div=long(r2)
                if div==0: div=10**(-5)
                delta = long(r1)/div
                err = dist(r1,r2)#/(long(r2)*f)
                #print("-({},{}),({},{}):Err={}".format(a,b,round(c),round(d),err))
                if err < minE1:
                    minE1 = err
                if (err < eps):
                    #print(1-delta)
                    if(abs(1-delta) < 0.07):
                        res.append([[a,b],[round(c),round(d)],delta,err])
                        if show:
                            print(">{:.3f}°:({},{})-({},{}): Delta={}%".format(th,a,b,round(c),round(d),delta*100),":",dist(r1,r2))
                # Buscando en b-
                if j!=0:
                    a,b = i,j
                    c = ((eq1*a)+(eq2*b))/(eq0)
                    d = ((eq3*a)+(eq4*b))/(eq0)
                    # Vector esperado
                    r1 = sumaV(multV(a,u),multV(b,v))
                    #Vector aproximado
                    r2 = sumaV(multV(round(c),rp),multV(round(d),rq))
                    # Calcula la proporción entre el tamaño del vector esperado y el aproximado
                    if((long(r1)!=0.0) & (long(r2)!=0.0)):
                        delta = long(r1)/long(r2)
                    err = dist(r1,r2)#/(long(r2)*f)
                    #print("--({},{}),({},{}):Err={}".format(a,b,round(c),round(d),err))
                    if err < minE2:
                        minE2 = err
                    if (err < eps):
                        if(abs(1-delta)<0.07):
                            res.append([[a,b],[round(c),round(d)],delta,err])
                            if show:
                                print(">>{:.3f}°:({},{})-({},{}): Delta={}%".format(th,a,b,round(c),round(d),delta*100),":",dist(r1,r2))
    if (minE1+minE2)/2 < 0.5:
        if show: print("----------\n{:.3f}°:{}\n\tdelta1={}\n\tdelta2={}\n----------".format(th,(minE1+minE2)/2,minE1,minE2))
    return res, ((minE1+minE2)/2)

#-------------------------Funciones Auxiliares para la función "Importa(File)"------------------------
def readFile(name):
    '''
    Lee el Archivo señalado y lo transforma en un arreglo de Strings de cada una de sus líneas
    '''
    file = open(name, 'r')
    count = 0
    lines = []
    while True:
        count+=1
        line = file.readline()
        if not line:
            break
        lines.append(format(line))
    file.close()
    return lines

def leeNumeros(linea):
    '''
    Lee un string y regresa un lista con los numeros en el
    '''
    s = []
    s = [float(l) for l in re.findall(r'-?\d+\.?\d*', linea)]
    return s
#----------------------------------------------------------------------------    
def importa(name):
    '''
    Importa un archivo vasp señalado por 'name' y lo transforma en una Red
    '''
    errormsg = '''El archivo no tiene el formato soportado por el programa.
Observe las características soportadas escribiendo <importa?>'''
    xyz = []
    try:
        # cargamos el nombre de la Red
        lines = readFile(name+".vasp")
        print("Se leerá el archivo {}".format(name+".vasp"))
        nameRed = lines[0].rstrip()
        #Cargamos los vectores primitivos de la Red
        vA = leeNumeros(lines[2])
        vB = leeNumeros(lines[3])
        vC = leeNumeros(lines[4])
        if (round(vA[2],5)!=0.0 or round(vB[2],5)!=0.0 or round(vC[0],5)!=0.0 or round(vC[1],5)!=0.0):
            msg = '''
            Vectores iniciales no soportados tal cómo están, se modificarán para tener la forma:
                a = (a1,a2,0), b = (b1,b2,0), c = (0,0,c3)
            Se deja discreción su funcionalidad.
            '''
            print(msg,errormsg)
            #raise SyntaxError('Documento no soportado')
        #Cargamos una lista con los tipos de Átomos en la Red y una con el numero de átomos de ese tipo
        aTipos = lines[5].split()
        aCant = leeNumeros(lines[6])
        if len(aTipos)!=len(aCant):
            print("No coincide el número de Tipos")
            raise SyntaxError('Documento no soportado')
        atomos = []
        ind = 8
        for i in range(len(aTipos)):
            col = '#'+''.join([random.choice('123456789ABCD') for i in range(6)])
            for j in range(round(aCant[i])):
                pA = leeNumeros(lines[ind+j])
                at = Atomo((pA[0], pA[1]), posZ=pA[2], color=col, sig=aTipos[i])
                atomos.append(at)
            ind = ind+round(aCant[i])
        leido = Red((vA[0],vA[1]),(vB[0],vB[1]),atms=atomos,name=nameRed,detachment=vC[2])
        print("--Red creada a partir del archivo '{}'--".format(name+".vasp"))
        return leido
    except FileNotFoundError:
        print("No se identifica el Archivo con el nombre '{}'".format(name))
    except IndexError:
        print(errormsg)
    except SyntaxError:
        print(errormsg)


#----------------------Redes prediseñadas------------------------------
def ejemplos():
    texto ='''
    Se cuenta con redes predefinidas, estas son:
hexa6(p,atms,name) -> Genera una Red hexagonal con constante de red 'p' y con los átomos de la lista atms.
    Si atms no se dá, entonces tendrá 2 átomos dentro de su base, generando una red hexagonal con 6 simetrías radiales.

hexa3(p,atms,name) -> Genera una Red hexagonal con constante de red 'p' y con los átomos de la lista atms.
    Si esta no se da entonces tendrá 2 átomos, uno dentro de su base y otro en un vértice, generando una red hexagonal con 3 simetrías radiales.

rectMesh(p1,p2,atms,name) -> Genera una Red rectangular con las constantes de red p1, p2 y los átomos señalados en la lista atms.
    Si esta no se dá, se generará con un solo átomo en el centro de su base.

grafeno() -> Genera una red de Grafeno, con constante de red 2.44 A y con sus átomos acomodados en el formato de hexa6

grafeno3() -> Genera una red de Grafeno, con constante de red 2.44 A y con sus átomos acomodados en el formato de hexa3

blackPhospho() -> Genera una Red de Fosforeno Negro con las constantes de red 3.3061099052 y 4.552418232.
'''
    print(texto)
    
def hexa6(p,atms=['C','C'],name=''):
    '''
    Genera una Red hexagonal con constante de red 'p' y con los átomos de la lista atms.
    Si esta no es dada entonces tendrá 2 átomos dentro de su base, generando una red hexagonal con 6 simetrías radiales.
    '''
    u,v=(p,0.0),(-p/2,math.sqrt(3)*(p/2))
    p1,p2,p3,p4 = (1/3,2/3),(2/3,1/3),(1/3,-1/3),(4/3,2/3)
    ats = [Atomo(p1, sig = atms[0]),Atomo(p2, sig = atms[1])]
    return Red(u,v,atms=ats,name=name,enls=[(p1,p2),(p2,p3),(p2,p4)])

def hexa3(p,atms=['C','C'],name=''):
    '''
    Genera una Red hexagonal con constante de red 'p' y con los átomos de la lista atms.
    Si esta no se da entonces tendrá 2 átomos, uno dentro de su base y otro en un vértice, generando una red hexagonal con 3 simetrías radiales.'''
    u,v=(p,0.0),(-p/2,math.sqrt(3)*(p/2))
    p1,p2,p3,p4 = (0.0,0.0),(1/3,2/3),(0,1),(1,1)
    ats = [Atomo(p1, sig = atms[0]),Atomo(p2, sig = atms[1])]
    return Red(u,v,atms=ats,name=name,enls=[(p1,p2),(p2,p3),(p2,p4)])

def rectMesh(p1,p2,atms='C',name=''):
    '''
    Genera una Red rectangular con las constantes de red p1, p2 y los átomos señalados en la lista atms.
    Si esta no se da, se generará con un solo átomo en el centro de su base.'''
    u,v = (p1,0.0),(0.0,p2)
    p1,p2,p3 = (1/2,1/2),(3/2,1/2),(1/2,3/2)
    ats = [Atomo(p1,sig = atms)]
    return Red(u,v,atms=ats,name=name,enls=[(p1,p2),(p1,p3)])
    
def grafeno():
    '''
    Genera una red de Grafeno, con constante de red 2.44 A y con sus átomos dentro de su base,
    generando una red hexagonal con 6 simetrías radiales.
    '''
    return hexa6(2.44, name='Grafeno')
    
def grafeno3():
    '''
    Genera una red de Grafeno, con constante de red 2.44 A y con uno de sus átomos dentro de su base y otro en un vértice,
    generando una red hexagonal con 3 simetrías radiales.
    '''
    return hexa3(2.44, name='Grafeno(s3)')

def blackPhospho():
    '''
    Genera una Red de Fosforeno Negro con las constantes de red 3.3061099052 y 4.552418232.
    '''
    a,b=(3.3061099052,0.0), (0.0,4.552418232)
    p1,p2=(0.000000000,0.913483083),(0.500000000,0.579813302)
    p3,p4=(0.000000000,0.079836130),(0.500000000,0.413437814)
    ats = [Atomo(p1,sig='P',posZ=0.266835123),Atomo(p2,sig='P',posZ=0.266945183),Atomo(p3,sig='P',posZ=0.181006327),Atomo(p4,sig='P',posZ=0.181094214)]
    enl = [(p1,p2),(p3,p4),(p2,p4),(p2,sumaV(p1,(1.0,0.0))),(p4,sumaV(p3,(1.0,0.0))),(p1,sumaV(p3,(0.0,1.0)))]
    return Red(a,b,atms=ats,name='Black-Phosphorene',enls=enl)
