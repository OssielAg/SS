## Primer Reporte

Se desea simular la representación de mayas de atomos con patrones cuadrados o exagonales con una periodicidad dada, para esto usamos el lenguaje Python auxiliado con Jupyter Notebook.

Tenemos varios métodos que usamos para cumplir nuestro proposito, para imprimir en pantalla (y en una imagen png si se requiere) usamos el método printPoints2(xList,yList,r), el cual imprime una lista de discos de radio r con centro en las posiciones dadas por las listas xList y yList, para generar una Maya de puntos con simetría cuadrada usamos el método rectangles(nX,nY,pX,pY) donde nX y nY son el numero de veces que repetiremos en las direcciónes x y y y pX y pY la distancia entre ellos. De igual forma el método hexagons(nX,nY,p) genera una Maya exagonal que se repite nX veces en la dirección x y nY veces en sus dirección en ángulo con un periodo P.

Pese a que se cumplió el objetivo, se vieron varias deficiencias, la primera es en el tiempo que tarda en graficarse la maya, además de este tambien está el hecho de que las Mayas creadas sólo son rectangulares (sus ángulos internos son rectos) o exagonales regulares, además no se pueden rotar. 