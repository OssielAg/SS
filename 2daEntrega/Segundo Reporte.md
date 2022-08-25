## Segundo Reporte

Los objetivos para esta entrega son 2 principalmente, el primero es crear una nueva representación de una maya usando sus Vectores constructores u y v pudiendo hacer que en esta maya se puedan representar tanto estructuras cuadradas cómo hexagonales pudiendo también cambiar el color de sus átomos y el segundo objetivo es poder rotar en un ángulo dado la maya.

En el nuevo código tenemos 2 clases, Átomo que será la representación de los átomos contenidos en la maya y Maya, que es la representación de la maya en si. La construcción de un objeto Átomo requiere de entrada su posición (x,y) y opcionalmente el color con el que será dibujado, un Objeto Maya requiere por entrada 2 vectores constructores 'u' y 'v' y opcionalmente una variable 'Theta' que indicará el ángulo de rotación de la maya si es que esta está rotada, si no se especifica se le dará un valor 0 por default. En el código se cuenta con ejemplos de su funcionamiento.

En la impresión de la maya en pantalla (y en un documento de imagen si así se indica) hay un cambio con respecto a la anterior usando la función scatter de matpltlib en lugar de patch reduciendo el tiempo de dibujo pero generando problemas en la representación cuando el tamaño relativo de los Átomos con la Maya completa es muy pequeño.
