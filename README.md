# hopfield
Código elaborado para el ejercicio voluntario del Modelo de Hopfield de red neuronal en la asignatura Física Computacional

El archivo "version1.c" es el programa básico para recordar patrones a partir de un fichero de texto externo. La temperatura del sistema es fija. Devuelve un fichero 
con la evolución del sistema y otro con el solapamiento con cada patrón introducido.
El archivo "version2.c" permite incluir un cierto rango de temperaturas (por defecto, 100 valores de 0.0001 a 0.5 en escala log10) para calcular el
solapamiento en función de la temperatura. 
El archivo "version3.c" calcula el solapamiento de patrones aleatorios (el número de patrones se debe introducir) y devuelve el número de patrones
que superan un solapamiento de 0.75.
