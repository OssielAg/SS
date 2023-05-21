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
        self.theta = th  #Ángulo entre capas de un sistema binario
        self.resultados = []#Lista de Puntos de red comunes para todas las capas
        self.loMat = []  #Lista de Matrices de transformación sugeridas
        self.MaxNumM = 10 #Número maximo de Matrices sugeridas
        self.SuperRed = None  #Super-Red del Sistema
        self.MT = None
    
    def its_hexagonal_system(self):
        '''
        Señala si el sistema está conformado sólo por redes hexagonales.
        '''
        err = 10**-6
        for r in self.redes:
            if abs(120.0-r.inAngle)>err:
                return False
        return True
    
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
            return False, -1
            a,b = transforma2v(self.redes[0].a,self.redes[0].b,M)
            c,d = transforma2v(self.redes[0].a,self.redes[0].b,self.loMat[i])
            if esRotacion(a,b,c,d):#Si el resultado es rotación del que ya tenemos lo descarta
                #print("\t--Matriz rotada existente")
                return False, -1
            b, i = self.goesHere(M,i+1)#Verifica en la siguiente posición
            return b, i
        #if dM%d==0:#Si el resultado es un multiplo del que tenemos lo descarta
            #print("\t--Matriz crecida {} veces".format(dM/d))
            #return False, -1
        b, i = self.goesHere(M,i+1)#Verifica en la siguiente posición
        return b, i
    
    def adjust(self, M,i,j):
        '''
        Verifica si una matriz de trasformación M entra a la lista loMat y la coloca manteniendo el orden
        M -> Matriz a examinar
        '''
        if det(M)==0:
            #print("---Matriz invalida--")
            return 0
        b, i = self.goesHere(M,0)
        if b:
            print("--Matriz agregada por resultados {} y {}".format(i,j))
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
    
    def calculateTM(self, s=True):
        '''
        Genera una lista de matrices de trasformación optimas de menor tamaño a partir
        de la lista 'resultados' del sistema.
        '''
        min_angle = 25 #Ángulo interno minimo aceptado para la super red calculada
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
                self.adjust(M,0,0)
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
                        self.adjust(M,i,j)
                #else:
                    #print(" - -Red generada con ángulo menor a {}°".format(min_angle))
        if s : print("Matrices guardadas en lista 'loMat'")
        return 1
        
    def createSuperLattice(self,M):
        '''
        Crea la red que representa al sistema.
        Los vectores primitivos de esta red se obtienen de transformar los vectores primitivos de la red de la
        primer capa con la matriz de transformación M dada.
        M -> Matríz de trasformación para los vectores primitivos
        '''
        sa, sb = transforma2v(self.redes[0].a,self.redes[0].b,M)
        self.SuperRed = superMesh(sa,sb,self.redes)
        self.MT = M
        return 1
    
    def searchLP(self, rangeOfSearch=15, epsilon=0.1, s=True):
        '''
        Busca las variables que desciben las transformaciones que debe tener cada Red del sistema para coincidir
        con una Super-Red común a todas ellas.
        rangeOfSearch -> Rango de busqueda de los posibles valores de transformaciòn de la red 1 (15 por Default)
        epsilon       -> Error máximo permitido para las transformaciones (0.1 por Default)
        '''
        lor = calculaPares(self.redes, max_val=rangeOfSearch, eps=epsilon)
        print("Puntos de red en comun encontrados:",len(lor))
        self.resultados = lor
        if lor==[]:
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
    
    def muestraSR(self, range_search=15, eps=0.1):
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
    
    def clon(self):
        redes = [copy.copy(r) for r in self.redes]
        s = system(redes)
        s.name = self.name
        s.theta = self.theta
        s.resultados = copy.copy(self.resultados)
        s.loMat = copy.copy(self.loMat)
        s.SuperRed = copy.copy(self.SuperRed)
        s.MT = self.MT
        return s