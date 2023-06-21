from .Funciones import *

class system:
    def __init__(self, lor, name="", th=0.0):
        '''
        Inicializa el 'sistema' formado por la lista de Redes descrita en 'lol' dandole el nombre 'name'.
        lol  -> Lista de redes que conforman el sistema
        name -> Nombre del sistema
        '''
        self.redes = lor #Lista de Redes del Sistema
        if name=="":
            self.name = 'Sistema ['+','.join([l.name for l in lor])+']' #Nombre del Sistema
        self.resultados = []#Lista de Puntos de red comunes para todas las capas
        self.loMat = []  #Lista de Matrices de transformación sugeridas
        self.MaxNumM = 10 #Número maximo de Matrices sugeridas
        self.SuperRed = None  #Super-Red del Sistema
        self.MT = None
    
    def calculateTM(self, min_angle = 25, s=True):
        '''
        Genera una lista de matrices de trasformación optimas de menor tamaño a partir
        de la lista 'resultados' del sistema.
        '''
        #Ángulo interno minimo aceptado para la super red calculada
        if self.resultados == []:
            errmsg = "No hay puntos de red en común para las capas."
            errmsg +="\nEjecute la función searchLP previamente."
            errmsg +="\n*En caso de ya haberlo hecho aumente el rango de búsqueda (rangeOfSearch) o el error mínimo aceptado (epsilon)"
            print(errmsg)
            return -1
        if s : print("Generando Matrices a partir de {} puntos".format(len(self.resultados)))
        self.loMat = []
        if self.its_hexagonal_system():
            print("Sistema Hexagonal")
            for r in self.resultados:
                [m,n]=r[0]
                M = vTm((m,n),(-n,m-n))
                self.adjust(M)
                if s: print("Matrices guardadas en lista 'loMat'")
            return 0
        for i in range(len(self.resultados)):
            for j in range(i+1,len(self.resultados)):
                #print("Arreglo ({},{})".format(i,j))
                [a1,a2] = self.resultados[i][0]
                [b1,b2] = self.resultados[j][0]
                M = vTm((a1,a2),(b1,b2))
                if det(M)<0:
                    M = vTm((b1,b2),(a1,a2))
                a,b = transforma2v(self.redes[0].a,self.redes[0].b,M)
                if cAng(a,b) >= min_angle:
                    if 180-cAng(a,b) >= min_angle:
                        #print("Analizando:",M)
                        self.adjust(M)
                #else:
                    #print(" - -Red generada con ángulo menor a {}°".format(min_angle))
        if s : print("Matrices guardadas en lista 'loMat'")
        return 1
        
    def createSuperLattice(self,M):
        '''
        Crea la supercelda que representa al sistema.
        Los vectores primitivos de esta red se obtienen de transformar los vectores primitivos de la red en la
        primer capa con la matriz de transformación M dada.
        M -> Matríz de trasformación para los vectores primitivos
        '''
        sa, sb = transforma2v(self.redes[0].a,self.redes[0].b,M)
        self.SuperRed = superMesh(sa,sb,self.redes)
        self.MT = M
        return 1
    
    def searchLP(self, rangeOfSearch=15, epsilon=0.05, p=True):
        '''
        Busca las variables que desciben las transformaciones que debe tener cada Red del sistema para coincidir
        con una Super-Red común a todas ellas.
        rangeOfSearch -> Rango de busqueda de los posibles valores de transformaciòn de la red 1 (15 por Default)
        epsilon       -> Error máximo permitido para las transformaciones (0.1 por Default)
        '''
        lor = calculaPares(self.redes, max_val=rangeOfSearch, eps=epsilon)
        if p:print("Puntos de red en comun encontrados:",len(lor))
        self.resultados = lor
        if lor==[]:
            msg = '''
            No se encontraron soluciones en el rango de búsqueda rangeOfSearch={} con el error mínimo epsilon={}
            Aumente el rango de búsqueda (rangeOfSearch) o el error mínimo aceptado (epsilon) para tener resultados
            '''.format(rangeOfSearch,epsilon)
            if p:print(msg)
            return -1
        return 0
        
    def analyze(self, rangeOfAngleSearch=(0.0,180.0), rangeOfSearch=15, precision=2, maxErr=0.05):
        '''
        Analiza el sistema buscando los ángulos de rotaciòn para red2 en los que sin modificar las redes el error es mínimo.
        rangeOfAngleSearch -> Rango de búsqueda para el ángulo señalado por el par (ángulo inicial, ángulo final)
        rangeOfSearch      -> Rango de busqueda de los posibles valores de transformaciòn de la red 1
        precision          -> Precición de la busqueda dada por la cantidad de dìgitos despues del punto en que se busca
                              (1=decimos de grado,2=centecimos de grado,...)
        maxErr             -> Error màximo tolerable para ser aceptado
        '''
        if len(self.redes)==2:
            return analiza(self.redes[0],self.redes[1],roAng=rangeOfAngleSearch,erMax=maxErr,mor=rangeOfSearch,accuracy=precision)
        print("Este método sólo funcional en Sistemas binarios")
       
    def muestraSR(self, range_search=15, eps=0.05):
        '''
        Ejecuta las funciones 'searchLP' y 'calculateTM'
        Con la primer matriz en la lista loMat calculada de esto ejecuta 'createSuperLattice'.
        Muestra el resultado en pantalla con la función 'show'.
        range_search -> Rango de busqueda de los posibles valores de transformaciòn de la red 1 (15 por Default)
        eps          -> Error máximo permitido para las transformaciones (0.1 por Default)
        '''
        self.searchLP(rangeOfSearch=range_search, epsilon=eps, s=False)
        if self.resultados == []:
            print("No se encontraron puntos de red comunes en todas las capas, intente cambiando los valores 'range_search' o 'eps'")
            return -1
        self.calculateTM(s=False)
        print("Posibles Matrices de trasformación calculadas:{}\nOpción recomendada:".format(len(self.loMat)))
        if self.loMat == []:
            return -2
        self.createSuperLattice(self.loMat[0])
        self.show()
        return 0
    
    def analiza_Mat(self):
        '''
        Analiza las posibles Matrices de Transformación de menor tamaño para el sistema, guardando sus caracteristicas en una lista.
        Tambien propone la mejor opción.
        '''
        if self.resultados==[]:
            self.searchLP(p=False)
            self.calculateTM(s=False)
        m0 = [[1.0,0.0],[0.0,1.0]]
        dif = 100.0
        info = []
        mejor = 0
        k=0
        a,b = self.redes[0].get_pv()
        for T in self.loMat:
            sa, sb = transforma2v(a,b,T)
            err = 0.0
            for i in range(1, len(self.redes)):
                inf_c = []
                ai,bi = self.redes[i].get_pv()
                T_i = corresponding_points(self.redes[0],self.redes[i],T)
                sa_i, sb_i = transforma2v(ai,bi,T_i)
                V_i = vTm(ai,bi)
                V_i_Op = m2M(vTm(sa,sb),inv2x2(T_i))
                S_i = m2M(inv2x2(V_i),V_i_Op)
                d_a, d_b = long(sa_i)/long(sa), long(sb_i)/long(sb)
                t_a, t_b = cAng(sa,sa_i), cAng(sb,sb_i)
                e_a, e_b = dist(sa,sa_i)/long(sa), dist(sb,sb_i)/long(sb)
                dd = calc_dd(V_i,V_i_Op)
                err += dd
                #[Mat Ti, Mat Si, error a, error b, def a, def b, ang a, ang b, dd, No Atms]
                inf_c.append([T_i,S_i,e_a,e_b,d_a,d_b,t_a,t_b,dd,det(T_i)*self.redes[i].nOAtms()])
            info.append([T,(sa,sb),err,inf_c])
            if err<dif:
                dif=err
                mejor = k
            k+=1
        return info,mejor
    
    def muestra(self):
        '''
        Muestra una tabla con las caracteristicas de los mejores resultados para las trasnformaciones sugeridas para
        calcular la supercelda que describa el sistema.
        Regresa la matriz T sugerida.
        '''
        cont = 1
        lista, optimo = self.analiza_Mat()
        for e in lista:
            [T,(sa,sb),err,inf_c] = e
            print("\n**Opción {}. T <- Matriz loMat[{}] del sistema".format(cont,cont-1))
            table = PrettyTable(["Red","T","Deformación","G de Distorsión","delta--theta","#Átomos"])
            print("Tamaño de los vectores primitivos:|a|={:.4f},|b|={:.4f}\nÁngulo entre vectores:{:.2f}°".format(long(sa),long(sb),cAng(sa,sb)))
            totalAtms = self.redes[0].nOAtms()*det(T)
            table.add_row(["\n" + self.redes[0].name,
                           mtoStr(T),
                           mtoStr([[1.0,0.0],[0.0,1.0]]),
                           "{:.8f}".format(0.0),
                           "{:.3f}% -- {:.2}°\n{:.2f}% -- {:.2}°".format(0.0,0.0,0.0,0.0),
                           totalAtms])

            for i in range(len(inf_c)):
                [T_i,S_i,e_a,e_b,d_a,d_b,t_a,t_b,dd,nAtm] = inf_c[i]
                totalAtms+=nAtm
                table.add_row(["\n" + self.redes[i+1].name,
                               mtoStr(T_i),
                               mtoStr(S_i),
                               "{:.8f}".format(dd),
                               "{:.3f}% -- {:.2}°\n{:.3f}% -- {:.2}°".format((d_a-1)*100,t_a,(d_b-1)*100,t_b),
                               nAtm])

            cont+=1
            print(table)
            print("\t\tTotal de Átomos:{}\tGrado de Deformación del sistema:{:.10f}".format(totalAtms,err))
        print("***Se recomienda usar la matriz de transformacion loMat[{}]***".format(optimo))
        return self.loMat[optimo]
    
    def optimize_system(self, T):
        '''
        Genera la Supercelda del sistema con las capas optimizadas aplicando tención en ellas de tal manera
        que no generen error.
        Los vectores primitivos de la supercelda se calculan a partir de la matríz T
        '''
        s = self.clon()
        if len(s.redes)<2:
            return -1
        V_0 = vTm(s.redes[0].a,s.redes[0].b)
        V_s = m2M(V_0,T)
        deformaciones = []
        for i in range(1,len(s.redes)):
            V_i = vTm(s.redes[i].a,s.redes[i].b)
            T_i = corresponding_points(s.redes[0], s.redes[i], T)
            V_i_Op = m2M(V_s,inv2x2(T_i))
            S_i = m2M(inv2x2(V_i),V_i_Op)
            #Guarda la matríz de deformación
            deformaciones.append(S_i)
            #Actualiza la red en la capa con los nuevos vectores primitivos
            a, b = mTv(V_i_Op)
            s.redes[i].setNewVectors(a,b)
            s.redes[i].name = s.redes[i].name + "(Opt)"
        #Creando super Red para el sistema deformado
        s.createSuperLattice(T)
        self.SuperRed = s.SuperRed
        self.MT = s.MT
        print("***La supercelda calculada está optimizada")
        self.show()
        return s, deformaciones
    
    #----------------------------------------------------------------------------------------
    
    def its_hexagonal_system(self):
        '''
        Señala si el sistema está conformado sólo por redes hexagonales.
        '''
        err = 10**-6
        for r in self.redes:
            if abs(120.0-r.inAngle)>err:
                return False
        return True
    
    def adjust(self, M):
        '''
        Verifica si una matriz de trasformación M entra a la lista loMat y la coloca manteniendo el orden
        M -> Matriz a examinar
        '''
        if det(M)==0:
            #print("---Matriz invalida--")
            return 0
        b, i = self.goesHere(M,0)
        if b:
            if len(self.loMat)<self.MaxNumM:
                if i<len(self.loMat):
                    self.loMat.insert(i,M)
                else:
                    self.loMat.append(M)
            else:
                self.loMat.insert(i,M)
                del self.loMat[-1]
        #print("---Matriz rechazada--")
        return 1
    
    def goesHere(self, M, i):
        '''
        Determina si una Matriz de transformacion M entra o no en la posición i de la lista loMat
        M -> Matriz a examinar
        i -> posición que examinamos
        '''
        if i>=len(self.loMat):#Si el indice sobrepasa los validos
            if len(self.loMat)<self.MaxNumM:#Verifica si aaun hay espacio en loMat
                return True,i
            #print("\t--Matriz descartada")
            return False, -1 #La matriz M es peor que las de loMat
        dM = det(M)
        d = det(self.loMat[i])
        if dM < d:# M es mejor opción que la que está actualmente en i
            #print("\t--Nueva mejor en {}".format(i))
            return True, i
        if dM == d:
            #return False, -1
            a,b = transforma2v(self.redes[0].a,self.redes[0].b,M)
            c,d = transforma2v(self.redes[0].a,self.redes[0].b,self.loMat[i])
            if esRotacion(a,b,c,d):#Si el resultado es rotación del que ya tenemos lo descarta
                #print("\t--Matriz rotada existente")
                return False, -1
            bo, i = self.goesHere(M,i+1)#Verifica en la siguiente posición
            return bo, i
        bo, i = self.goesHere(M,i+1)#Verifica en la siguiente posición
        return bo, i
    
    def clon(self):
        '''
        Clona el sistema regresando una copia identica
        '''
        redes = [copy.copy(r) for r in self.redes]
        s = system(redes)
        s.name = self.name
        s.resultados = copy.copy(self.resultados)
        s.loMat = copy.copy(self.loMat)
        s.SuperRed = copy.copy(self.SuperRed)
        s.MT = self.MT
        return s
    
    def set_maximum_number_of_matrices(self, newMax):
        self.MaxNumM = newMax
        self.calculateTM(s=False)
        
    
    def show(self):
        if self.SuperRed is None:
            print('*'*84+"\n  Super Red no calculada aún, utilice la función 'createSuperLattice' para hacerlo  \n"+'*'*84)
        else:
            print("Matriz de trasformación:")
            pmat(self.MT)
            print(self.name,"\nCelda unitaria:")
            self.SuperRed.showme()
            print("Espacio Reciproco:")
            self.SuperRed.printReciprocalSpace()
    