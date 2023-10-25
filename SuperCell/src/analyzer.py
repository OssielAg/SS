from .Lattice import *
from prettytable import PrettyTable

#Funciones para analizar un par de redes y encontrar las mejores opciones de transformaciones para crear una Suer-Red que las contenga.

def calculaEM(r1, r2, th = 0.0, maxIt=15):
    '''
    Calcula los pares enteros (a,b) y (c,d) con a y b < maxIt tales que, si: u,v son los vectores primitivos de 'r1' y rp,rq los
    vectores primitivos de 'r2' rotada en 'th' grados; entonces: P1 = (au + bv), P2 = (c(rp) + d(rq)) son cercanos.
    
    Regresa tambien el promedio de los las diferencias mínimas en los cuadrantes I y IV del plano cartesiano que cumplen lo anterior.
    '''
    (u,v), (p,q) = r1.getVectors(), r2.getOV()
    (u_1,u_2), (v_1,v_2) = u, v
    (p_1,p_2), (q_1,q_2) = rp, rq = rota(p,th), rota(q,th)
    puntos= []
    rango = maxIt
    minim = 10**9
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
                err = dist(r1,r2)
                if err < 0.6:
                    puntos = acomoda(r1,err,puntos,2)
                # Buscando en b-
                if j!=0:
                    a,b = i,j
                    c = ((eq1*a)+(eq2*b))/(eq0)
                    d = ((eq3*a)+(eq4*b))/(eq0)
                    # Vector esperado
                    r1 = sumaV(multV(a,u),multV(b,v))
                    #Vector aproximado
                    r2 = sumaV(multV(round(c),rp),multV(round(d),rq))
                    err = dist(r1,r2)/2
                    if err < 0.6:
                        puntos = acomoda(r1,err,puntos,2)
    if len(puntos)>=2:
        minim = (puntos[0][1]+puntos[1][1])/2
    return minim

def explora(r1, r2, mIt=15, eMax=0.5, thI=0.0, thF=180.0, acc=1):
    '''
    Hace una exploracion para angulos entre 'thI' y 'thF' en un intervalo de 10^(-'acc') grados para obtener
    una estimación del error minimo para pares de enteros menores a 'mIt' del sistema de coordenadas con base
    en los vectores generadores de 'r1' y sus contrapartes enteras en el sistema de coordenadas con base en
    los vectores generadores de 'r2'.
    '''
    graf = []
    i = round(thI*(10**acc))
    f = round(thF*(10**acc))
    print("Analizando para theta en intervalo [{}°,{}°]".format(i/(10**acc),f/(10**acc)))
    print(".............")
    for t in range(i,f+1):
        theta = t/(10**acc)
        prc=round((t/(f-i))*100)
        b = "Analizando..."+str(prc)+"%"
        print(b,end="\r")
        m = calculaEM(r1, r2, th = theta, maxIt=mIt)
        if m < eMax:
            graf.append([theta,m])
    print(end="\r")
    print('**********Exploración finalizada**********')
    return np.array(graf)
   
def analiza(r1,r2,roAng=(0.0,180.0),erMax=0.005,mor=15,accuracy=2):
    '''
    Ejecuta la función 'explora' para las redes 'r1' y 'r2' para thetas en el intervalo 'roAng' con un mIt igual a 'mor'
    y un acc igual a accuracy.
    Analiza los resultados imprimiendo una gráfica de la relación Theta vs error_minimo_aproximado y regresa una lista
    con los ángulos cuyo error_mínimo_aproximado es menor a 'erMax'.
    Tambien regresa la lista resultante de la función 'explora'
    '''
    (i,f) = roAng
    graf = explora(r1, r2, mIt=mor, thI=i, thF=f, acc=accuracy)
    xs, ys = graf[:,0], graf[:,1]
    plt.plot(xs, ys)
    plt.show()
    resultados=[]
    analisis=np.r_[True, ys[1:] < ys[:-1]] & np.r_[ys[:-1:] < ys[1:], True]
    print("**********\nLos ángulos con los errores mínimos son:")
    for i in range(len(analisis)):
        if analisis[i]==True:
            if(ys[i]<erMax):
                resultados.append([xs[i],ys[i]])
                print("\tTh={:.3f} : Er.Aprx={:.5f}A".format(xs[i],ys[i]))
    return resultados,graf

def muestra(lor,r1,r2):
    '''
    Muestra una tabla con las caracteristicas de los mejores resultados para las trasnformaciones requeridas para la super Red que describa el sistema r2 sobre r1.
    lor -> Lista de resultados obtenidos de la funciòn 'analiza'
    r1  -> Red1
    r2  -> Red2
    '''
    leyend = "Mejores candidatos para Super-Red:'R1={}',R2='{}'\n".format(r1.name,r2.name)
    print(leyend)
    if round(r1.inAngle,3)==120.000 and round(r1.inAngle,3)==120.000:
        table = PrettyTable(["# de Atms","T1","T2","Tención/Rotación Red1","Tención/Rotación Red2", "Err"])
        for c in lor:
            noa = r1.nOAtms()*det(dameTH(c[0]))+r2.nOAtms()*det(dameTH(c[1]))
            sa,sb = mTv(m2M(vTm(r1.a,r1.b),dameTH(c[0])))
            md1,md2=ajusta(r1,r2,dameTH(c[0]),dameTH(c[1]),sa,sb)
            text1 = describeAjuste(r1,md1)
            text2 = describeAjuste(r2,md2)
            table.add_row([noa,
                           mtoStr(dameTH(c[0])),
                           mtoStr(dameTH(c[1])),
                           text1,
                           text2,
                           str(round(c[3],10))])
        print(table)
        return table
    return []

def ajusta(r1,r2,t1,t2,p1,p2):
    '''
    Calcula el ajuste requerido en los vectores primitivos de las redes dadas para que al transformarlos con su respectiva matriz 't' coincida con las posiciones p1 y p2 esperadas.
    r1 -> Red1
    r2 -> Red2
    t1 -> Transformación para Red1
    t2 -> Transformaciòn para Red2
    p1 -> Vector A de la super red que corresponde a Red2 sobre Red1
    p2 -> Vector B de la super red que corresponde a Red2 sobre Red1
    '''
    mr1i = inv2x2(vTm(r1.a,r1.b))
    mr2i = inv2x2(vTm(r2.a,r2.b))
    mres = vTm(p1,p2)
    mt1i = inv2x2(t1)
    mt2i = inv2x2(t2)
    md1 = m2M(mr1i,m2M(mres,mt1i))
    md2 = m2M(mr2i,m2M(mres,mt2i))
    return md1, md2  


#-----------------Funciones Auxiliares para la funciòn 'muestra'.-----------------
def strS(i,n=4):
    '''
    Combierte i en un string de al menos tamaño n rellenando con espacios en blanco
    '''
    if type(i)==float:
        s = str(round(i,5))
        n = max(8,n)
    else:
        s = str(i)
    while len(s)<n:
        s=" "+s
    return s

def dameTH(mn):
    '''
    Dada una lista con 2 valores m y n, genera la matriz de transformación nesesaria para un elemento de un sistema de 2 redes hexagonales
    '''
    return [[mn[0],-mn[1]],[mn[1],mn[0]-mn[1]]]

def mtoStr(t):
    '''
    Combierte una matriz de Transformacion t en texto.
    '''
    sm = ""
    for f in range(len(t)):
        l="|"
        for c in range(len(t[f])):
            if c!=0: l=l+" "
            l=l+strS(t[f][c])
        sm = sm + l + "|\n"
    return sm

def top(lor, n=5):
    '''
    Recibe la lista de resultados arrojada por la función 'calculaPares' y regresa una lista con los 'n' resultados con el menor valor de error.
    '''
    try:
        if len(lor[0][0])<n:
            lis = lor[0][0]+lor[0][1]
        lis = lor[0][0]
        if len(lis)<5: n = len(lis)
        lo = sorted(lis, key=lambda op : op[3])[0:n]
    except:
        print("La lista dada no proviene de la operación 'calculaPares'")
    return lo


#--------------------Funciones Auxiliares para 'ajuste'--------------------
def leeAjuste(red,ajuste):
    ajA,ajB = mTv(m2M(vTm(red.a,red.b),ajuste))
    ajA1,ajB1 = ((long(ajA)/long(red.a))*100)-100, ((long(ajB)/long(red.b))*100)-100
    rotA,rotB = cAng((1,0),ajA)-cAng((1,0),red.a),cAng((1,0),ajB)-cAng((1,0),red.b)
    return ajA1,ajB1,rotA,rotB

def describeAjuste(red,ajuste):
    aA,aB,rA,rB=leeAjuste(red,ajuste)
    text=""
    if round(aA,5)==round(aB,5):
        text=text+"{}%".format(round(aA,5))
    else:
        text=text+"({}%, {}%)".format(round(aA,5),round(aB,5))
    if round(rA,5)==round(rB,5):
        text=text+"\n{}°".format(round(rA,5))
    else:
        text=text+"\n({}°, {}°)".format(round(rA,5),round(rB,5))
    return text

