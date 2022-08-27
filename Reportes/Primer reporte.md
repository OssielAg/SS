## Primer Reporte

El objetivo de la primer entrega es simular la representación de mallas de átomos con patrones cuadrados o hexagonales con una periodicidad dada con el fin de observar la nueva malla formada al colocar una sobre otra, para esto usamos el lenguaje Python auxiliado con Jupyter Notebook.

Se presentan varios métodos que usamos para cumplir nuestro propósito, para imprimir en pantalla (y en una imagen png si se requiere). Tenemos el método printPoints2(pList1,pList2,r), el cual dibuja en una sola imagen discos de radio r con centro en las posiciones dadas por las listas pList1 y pList2 respectivamente. Para generar una lista con la posición de puntos con simetría cuadrada usamos el método rectangles(nX,nY,pX,pY) donde nX y nY son el numero de veces que repetiremos en las direcciones 'x' y 'y', pX y pY la distancia entre ellos. De igual forma el método hexagonal(nX,nY,p) genera una lista de posiciones pero con estructura hexagonal que se repite nX veces en la dirección x y nY veces en sus dirección en ángulo con un periodo P.

Pese a cumplir el objetivo, se vieron varias deficiencias, la primera es en el tiempo que se toma en graficar la maya, otro problema está en el hecho de que las Mallas creadas se limitan a ser sólo rectangulares (sus ángulos internos son rectos) o hexagonales regulares, además no de no poderse rotar.
