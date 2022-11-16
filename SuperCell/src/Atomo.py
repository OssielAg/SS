from .Basicas import *
class Atomo:
    '''Clase que define un atomo en el sistema'''
    def __init__(self, pos, posZ=0.0, color='black', sig='C', lvl=1):
        '''Inicializa el objeto Atomo usando como entradas su posición (x,y) pos y su color caracteristico color'''
        (x,y) = pos
        self.pos = (x,y)
        self.posZ = posZ
        self.color = color
        self.sig = sig
        self.lvl = lvl
    
    def __str__(self):
        return self.sig + str(self.pos)

    def printAtomo(self, r, axs):
        '''Imprime un Atomo en un Patch en el axs señalado'''
        axs.add_patch(plt.Circle((self.pos), r , color=self.color))

    def setData(self, color='Silver', sig='C'):
        '''Cambia los datos del color y/o signo del Atomo'''
        self.color = color
        self.sig = sig
    
    def getPos(self):
        return self.pos
        
    def clasifica(self, lols):
        '''Clasifica al Atomo en una lista de atomos distinguidos por su signo "sig" '''
        if len(lols)==0:
            lols.append([self])
            return 1
        for l in lols:
            if len(l)==0:
                l.append(this)
                return 1
            else:
                if l[0].sig==self.sig:
                    l.append(self)
                    return 1
        lols.append([self])
        return 0