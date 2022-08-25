## Primer Reporte

Se desea simular la representación de mayas de átomos con patrones cuadrados o hexagonales con una periodicidad dada, para esto usamos el lenguaje Python auxiliado con Jupyter Notebook.

Tenemos varios métodos que usamos para cumplir nuestro propósito, para imprimir en pantalla (y en una imagen png si se requiere) usamos el método printPoints2(xList,yList,r), el cual imprime una lista de discos de radio r con centro en las posiciones dadas por las listas xList y yList, para generar una Maya de puntos con simetría cuadrada usamos el método rectangulares(nX,nY,pX,pY) donde nX y nY son el numero de veces que repetiremos en las direcciones 'x' y 'y', pX y pY la distancia entre ellos. De igual forma el método hexagonal(nX,nY,p) genera una Maya hexagonal que se repite nX veces en la dirección x y nY veces en sus dirección en ángulo con un periodo P.

Pese a que se cumplió el objetivo, se vieron varias deficiencias, la primera es en el tiempo que tarda en graficarse la maya, además de este también está el hecho de que las Mayas creadas sólo son rectangulares (sus ángulos internos son rectos) o hexagonales regulares, además no se pueden rotar. 