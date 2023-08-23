from .Atomo import *

class Red:
    def __init__(self, vA, vB, atms=[], enls=[], name='',detachment=10, prof=1):
        '''
        Inicializa una Red con los vectores, su rotacion, lista de Atomos, lista de sus Enlaces, Nombre y nivel que ocupa en el apilamiento.
        Sus datos son...
            detachment: Tamaño del vector C de la Red
            name: Nombre de la Red
            OriginalA: Vector A sin rotación
            OriginalB: Vector B sin rotación
            a: Vector A (con la rotación establecida)
            b: Vector B (con la rotación establecida)
            theta: Ángulo de rotacion de la Red
            prof: Profundidad con respecto a la pila de Redes dada
            inAngle: Ángulo entre ambos vectores generadores de la celda Unitaria
            atms: Lista de Átomos en la celda Unitaria
            enls: Lista de enlaces en la celda Unitaria
        '''
        lols = []
        self.enls=[]
        (u1,u2),(v1,v2) = vA,vB
        ang = abs(cAng(vA,vB))
        
        rAngle = cRot(vA)
        if u2<0: rAngle=-rAngle
        self.detachment = detachment
        
        self.a = vA
        self.b = vB
        self.OriginalA = rota(vA,-rAngle)
        self.OriginalB = rota(vB,-rAngle)
        self.prof = prof
        self.theta = cAng((1,0),vA)
        self.enls = enls
        self.inAngle = ang
        self.layerList = []
        a1,a2 = vA
        b1,b2 = vB
        self.reciprocalVectors = cRecip((a1,a2,0.0),(b1,b2,0.0),(0.0,0.0,detachment))
        
        for a in atms:
            a.clasifica(lols)
        self.atms = lols
        
        if name=='':
            self.name = "Lattice({:.3f}°)".format(ang)
        else:
            self.name = name
    
    def __str__(self):
        return self.showData()
    
    def get_pv(self):
        '''
        Regresa los vectores primitivos de la red
        '''
        return self.a, self.b
    
    def showXY(self, x, y, x0=0 ,y0=0, t=10,name=''):
        '''
        Muestra en pantalla la imagen resultante de repetir la celda unitaria de la Red 'x' veces en su vector 'a' y
        'y' veces en su vector 'b'
        - x0 y y0 corresponden a la posición del origen de la imagen y t el grosor de dibujo de los atomos.
        - Si se da un nombre se generará una imagen PNG con ese nombre. 
        '''
        red=[]
        lx1,ly1 = getLim(self.a,self.b,x,y)
        lx2,ly2 = getLim(self.a,self.b,x,y0)
        lx3,ly3 = getLim(self.a,self.b,x0,y)
        lx4,ly4 = getLim(self.a,self.b,x0,y0)
        lmsx = [min(lx1[0],lx2[0],lx3[0],lx4[0]),max(lx1[1],lx2[1],lx3[1],lx4[1])]
        lmsy = [min(ly1[0],ly2[0],ly3[0],ly4[0]),max(ly1[1],ly2[1],ly3[1],ly4[1])]
        difMax = max((lmsx[1]-lmsx[0]),(lmsy[1]-lmsy[0]))
        '''
        lx1,ly1 = getLim(self.a,self.b,x,y)
        mi,ma = min(lx1[0],ly1[0]), max(lx1[1],ly1[1])
        lx2,ly2 = getLim(self.a,self.b,x,y0)
        mi,ma = min(lx2[0],ly2[0],mi), max(lx2[1],ly2[1],ma)
        lx3,ly3 = getLim(self.a,self.b,x0,y)
        mi,ma = min(lx3[0],ly3[0],mi), max(lx3[1],ly3[1],ma)
        lx4,ly4 = getLim(self.a,self.b,x0,y0)
        mi,ma = (min(lx4[0],ly4[0],mi)-1), (max(lx4[1],ly4[1],ma)+1)'''
        ats=[]
        col=[]
        enls=[]
        #self.enls=[]
        for i in range(abs(x-x0)):
            a=i+x0
            for j in range(abs(y-y0)):
                b=j+y0
                #Carga contornos de celda
                red.append([m2V(self.a,self.b,(a,b)),m2V(self.a,self.b,(a+1,b))])
                red.append([m2V(self.a,self.b,(a,b)),m2V(self.a,self.b,(a,b+1))])
                red.append([m2V(self.a,self.b,(a+1,b)),m2V(self.a,self.b,(a+1,b+1))])
                red.append([m2V(self.a,self.b,(a,b+1)),m2V(self.a,self.b,(a+1,b+1))])
                #Carga los atomos
                for c in self.atms:
                    for at in c:
                        (pu,pv) = at.pos
                        na = m2V(self.a,self.b,(a+pu,b+pv))
                        ats.append(na)
                        col.append(at.color)
                #Carga los enlaces
                for (ei,ef) in self.enls:
                    (ei1,ei2) = ei
                    (ef1,ef2) = ef
                    o = m2V(self.a,self.b,(a+ei1,b+ei2))
                    f = m2V(self.a,self.b,(a+ef1,b+ef2))
                    enls.append([o,f])
        fig, maxs = plt.subplots()
        lis = np.array(ats.copy())
        xs, ys = lis[:,0], lis[:,1]
        #Cargamos los contornos de las celdas en una lista de lineas
        lr = mc.LineCollection(np.array(red), colors='silver', linewidths=(t/30))
        #Cargamos los enlaces en una lista de lineas
        lc = mc.LineCollection(np.array(enls), colors='black', linewidths=(t/20))
        #Dibuja los Enlaces
        maxs.add_collection(lc)
        #Dibuja los contornos de las celdas
        maxs.add_collection(lr)
        #Dibuja los Atomos
        maxs.scatter(xs,ys, color=col,s=t*3)
        maxs.axes.xaxis.set_visible(False)
        maxs.axes.yaxis.set_visible(False)
        medX = ((lmsx[0]+lmsx[1])/2)
        medY = ((lmsy[0]+lmsy[1])/2)
        maxs.set(xlim=(medX-+(difMax/2),medX+(difMax/2)), ylim = (medY-+(difMax/2),medY+(difMax/2)))
        #Dibuja la escala de la imagen
        s=max(round(difMax/4),1)
        f=s/difMax
        si=0.1
        sf=si+f
        plt.text((medX-(0.75*(difMax/2))),
                 (medY-(0.74*(difMax/2))),
                 "{} nm".format(s/10),
                 fontsize=14,
                 weight='bold',
                 c='royalblue',
                 backgroundcolor='white')
        plt.axhline(y=(medY-(0.8*(difMax/2))), xmin=si, xmax=sf, c='white',lw=8.0)
        plt.axhline(y=(medY-(0.8*(difMax/2))), xmin=si, xmax=sf, c='royalblue',lw=4.0)
        
        if name!='':
            plt.savefig(('imagenes/'+name),dpi=900, bbox_inches='tight')
        plt.show()
        return 1
        
    def showme(self):
        '''
        Muestra la celda unitaria de la red en pantalla
        '''
        self.showXY(1,1)
    
    def addAtms(self, lols):
        '''
        Agrega de forma ordenada una lista de Átomos a la Red
        '''
        for l in lols:
            for a in l:
                a.clasifica(self.atms)
        return 1
        
    def showData(self):
        '''
        Muestra la Red en un Texto con el formato de VASP en coordenadas fraccionales 
        '''
        ang = cAng((1,0),self.a)
        (u1,u2) = self.OriginalA
        (v1,v2) = self.OriginalB
        data=self.name+'''
1.0
        {:.10f}         {:.10f}         0.0000000000
        {:.10f}         {:.10f}         0.0000000000
        0.0000000000         0.0000000000         {:.10f}
'''.format(u1,u2,v1,v2,(self.detachment))
        atms1=""
        atms2=""
        atms3="Direct"
        for i in range(len(self.atms)):
            atms1=atms1+"\t{}".format(self.atms[i][0].sig)
            atms2=atms2+"\t{}".format(len(self.atms[i]))
            for atm in self.atms[i]:
                (x,y)=atm.pos
                atms3=atms3+"\n     {:.10f}         {:.10f}         {:.10f}".format(abs(x),abs(y),(atm.posZ + 0.0001))
        atms = atms1+"\n"+atms2+"\n"+atms3
        data = data+atms
        return data
    
    def aligned(self):
        '''
        Alinea una Red con el eje X
        '''
        self.a = self.OriginalA
        self.b = self.OriginalB
        self.theta = 0.0
        
    def rotate(self, th):
        '''
        Rota la Red tal cual está th grados.
        '''
        self.a = rota(self.OriginalA,self.theta + th)
        self.b = rota(self.OriginalB,self.theta + th)
        self.theta = self.theta + th
    
    def mAlig(self):
        '''
        Regresa una copia de la Red alineada al eje X
        ''' 
        na, nb = self.OriginalA, self.OriginalB
        natms, nenls = self.atms.copy(), self.enls.copy()
        nName = self.name+"(aligned)"
        mr = Red(na, nb, enls=nenls, name=nName, prof=self.prof)
        mr.atms = self.atms
        mr.detachment = self.detachment
        return mr
    
    def mRot(self, ang):
        '''
        Regresa una copia de la Red rotada en "ang" grados
        ''' 
        na, nb = rota(self.a,ang), rota(self.b,ang)
        natms, nenls = self.atms.copy(), self.enls.copy()
        nName = self.name+"(rot {}°)".format(ang)
        mr = Red(na, nb, enls=nenls, name=nName, prof=self.prof)
        mr.atms = self.atms
        mr.detachment = self.detachment
        return mr
    
    def exporta(self,name=''):
        '''
        Exporta la red a un archivo VASP. Si no se da un nombre utiliza el de la red.
        '''
        if name=='':
            name=self.name
        name = "ArchivosVASP/" + name + ".vasp"
        f = open(name,"w")
        f.write(self.showData())
        f.close()
        return 1
    
    def setNewVectors(self, newA, newB):
        '''
        Cambia los Vectores generadores de la Red
        '''
        th = cAng((1,0),newA)
        self.theta = th
        self.a, self.b = newA, newB
        self.OriginalA = rota(newA, -th)
        self.OriginalB = rota(newB, -th)
        self.inAngle = abs(cAng(newA, newB))
        a1,a2 = newA
        b1,b2 = newB
        self.reciprocalVectors = cRecip((a1,a2,0.0),(b1,b2,0.0),(0.0,0.0,self.detachment))
        
    def getVectors(self):
        '''
        Regresa los Vectores generadores de la Red
        '''
        return (self.a, self.b)
    
    def getOV(self):
        '''
        Regresa los Vectores generadores sin rotación de la Red
        '''
        return (self.OriginalA, self.OriginalB)
    
    def nOAtms(self):
        '''
        Regresa el número de átomos que tiene la red
        '''
        noa = 0
        for loa in self.atms:
            noa = noa + len(loa)
        return noa

    def printReciprocalSpace(self, t=10, border=1.0,prnt=False):
        '''
        Imprime en pantalla la FBZ de la red en el espacio reciproco, si esta red pertenece a un sistema multicapa
        imprime tambien la FBZ de cada capa.
        t      -> Valor con el que se determinarán los tamaños con el que se dibujarán puntos y líneas.
        border -> Límites de la gráfica dibujada.
        '''
        lol = self.layerList
        print("Calculando...")
        ax = plt.subplot()
        #"Pintamos" la FBZ de cada capa que forma la red si es que es un sistema multicapa
        if len(lol) > 0:
            colors=['#'+''.join([random.choice('3456789') for i in range(6)]) for j in range(len(lol))]
            for i in range(len(lol)):
                vl, eq = calcVerticesFBZ(lol[i])
                fbzLayer = Polygon(vl, alpha=0.4, color = colors[i], label = lol[i].name)
                ax.add_patch(fbzLayer)
                print("...Pintando capa {} ({})".format(i+1,lol[i].name))
        #"Pntamos" la FBZ de la red y el fondo dado por la función 'reciprocalBackgroundMesh'
        vl, eq = calcVerticesFBZ(self)
        xs, ys, linkList = reciprocalBackgroundMesh(self,vl,t)
        fbzRed = Polygon(vl, alpha=0.7, color = 'gray', label = "Super Red")
        ax.add_patch(fbzRed)

        ax.add_collection(linkList)#"Pintamos" los enlaces de la red de fondo calculado previamente
        ax.scatter(xs,ys, color='black',s=t)#"Pintamos" los Puntos de la red reciproca calculados previamente 
        ax.set(xlim=(-border,border), ylim=(-border,border))
        ax.legend(loc = 'upper right')
        print("...Terminado")
        if prnt:
            plt.savefig(('imagenes/SuperRed(RS)'),dpi=900, bbox_inches='tight')
        plt.show()
        return 1