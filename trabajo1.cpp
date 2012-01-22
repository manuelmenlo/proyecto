/*Este es un programa que busca simular el enfriamiento por evaporación
de un gas cuyos constituyentes obedecen la estadística de Bose - Einstein 
y además interactúan. Dicho gas está confinado por una trampa armónica bidimensional
Resulta ser que la condensación de Bose - Einstein es una transición de fase
y entonces, se puede demostrar, que son dos las ecuaciones que simulan la dinámica
del sistema: una ecuación se obedece cuando la temperatura es mayor a la crítica
y otra ecuación se obedece para temperaturas menores a la crítica.
Este programa intenta resolver numéricamente la ecuación para temperaturas
mayores a la crítica.
Estrictamente hablando, la dinámica queda simulada por un conjunto de
ecuaciones diferenciales acopladas, cada una simula el comportamiento temporal
del número de ocupación promedio en cada estado.*/

/*El método numérico a utilizar será el método de Euler.*/

/*A continuación se incluyen las librerías que se van a utilizar a lo largo del 
programa. En particular, "Extremos.h" es una librería que permite determinar
el máximo o el mínimo de un conjunto de números y "Aleatorios.h" es una librería
que permite calcular números aleatorios en un intervalo.*/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include "Extremos.h"
#include "Aleatorios.h"

using namespace std;

const int N = 16; /*El número de niveles energéticos.*/
const int E = 136; /*El número de estados posibles por los que pueden pasar las partículas.*/
const float g = 0.00007; /*La magnitud de la interacción.*/
const float dt = 0.0000005; /*El paso temporal.*/
const int iteraciones = 20; /*El tiempo total de la simulación.*/

int main () {
	
	ifstream entrada1;
	entrada1.open("interacciones.txt"); /*Archivo que contiene los términos de interacción, es decir, las integrales.*/
	ifstream entrada2;
	entrada2.open("edos_posibles.txt"); /*Archivo que contiene todos los estados posibles de las partículas.*/
	ofstream salida;
	salida.open("datos.txt"); /*Archivo que guardará el valor de los números de ocupación promedio.*/

	/*A continuación se específica cada una de los arreglos: TI será el arreglo que contendrá algunos
	de los términos de interacción, MTI será la matriz que incluirá todas las interacciones posibles. 
	Ahora bien, ya que se usará el método de Euler, será evidente que N0 es la condición inicial
	y N1 es la solución en el siguiente paso temporal. S es la función que se calcula
	en cada iteración temporal y que es multiplicada por el paso de tiempo dt; esta función
	S es una suma.*/

	/*Al tratarse de un problema en dos dimensiones, resulta ser que cada estado es en 
	realidad un vector y entonces las "x's" estarán en el arreglo I1 y las "y's" 
	estarán en el arreglo I2. En tanto, en la variable ENE se guardará la 
	energía correspondiente de cada estado.*/

	float *TI, *MTI, *N0, *N1, *S;
	int *I1, *I2, *ENE;

	/*Se reserva memoria en el host para cada uno de los arreglos.*/

	TI = new float[3876]; /*Son 3876 términos los que se encuentran en el archivos interacciones.txt.*/
	MTI = new float[N*N*N*N];
	N0 = new float[E];
	N1 = new float[E];
	S = new float[E];
	I1 = new int[E];
	I2 = new int[E];
	ENE = new int[E];

	/*A continuación se leerán los términos de interacción de la entrada1 y los
	valores se almacenarán en la variable TI.*/

	for (int i = 0; i < 3876; i++) {
		entrada1 >> TI[i];
	}
	entrada1.close();

	/*Estrictamente hablando no se han calculado todos los términos de interacción, sino más bien sólo algunas permutaciones que consisten en: 
	dados cuatro números, éstos se ordenan de menor a mayor y se determina el valor de la integral; ésta es una de las permutaciones. 
	Para clarificar lo anterior se tomará el siguiente ejemplo: dados los números (1,2,1,1), las permutaciones posibles son (1,1,1,2), (1,1,2,1) y (2,1,1,1). 
	Sin embargo la integral que se calculó fue la correspondiente a los índices (1,1,1,2).1
	Con el siguiente conjunto de for's se busca construir la matriz MTI que
	incluye todos los términos de interacción.*/

	int iter = 0;
	int min, max, norma;
	for(int i = 0; i < N; i++) {
		for(int j = i; j < N; j++) {
			for(int k = j; k < N; k++) {
				for(int l = k; l < N; l++) {
					min = MIN(i,j,k,l); max = MAX(i,j,k,l); norma = i*i+j*j+k*k+l*l;
					for(int ip = min; ip <= max; ip++) {
						for(int jp = min; jp <= max; jp++) {
							for(int kp = min; kp <= max; kp++) {
								for(int lp = min; lp <= max; lp++) {
									if((ip==i||ip==j||ip==k||ip==l)&&(jp==i||jp==j||jp==k||jp==l)
									&&(kp==i||kp==j||kp==k||kp==l)&&(lp==i||lp==j||lp==k||lp==l)
									&&(ip*ip+jp*jp+kp*kp+lp*lp==norma)) {
										MTI[ip*N*N*N+jp*N*N+kp*N+lp] = TI[iter]; /*MTI es un arreglo que ha sido mapeado en un arreglo lineal. La transformación es MTI[ip][jp][kp][lp] --> MTI[ip*N*N*N+jp*N*N+kp*N+lp].*/
									}
								}
							}
						}
					}
					iter++;
				}
			}
		}
	}

	/*A continuación se lee el archivo "edos_posibles.txt" la primer columna contiene
	los valores que irán en el arreglo I1, la segunda los correspondientes a la 
	variable I2 y la tercera columna del archivo son las energías que se 
	guardarán en el arreglo ENE.*/

	/*Puse este comentario al final para seguir usando github.*/ /*Cambié de lugar el comentario que puse al final.*/

	for (int i = 0; i < E; i++) {
		entrada2 >> I1[i] >> I2[i] >> ENE[i];
	}

	/*Lo que sigue ahora es establecer las condiciones iniciales en los
	números de ocupación promedio.*/

	float par0 = 0.0; /*par0 es el número de partículas inicial.*/
	float ene0 = 0.0; /*ene0 es la energía total del sistema inicial.*/
	N0[0] = 0.0;
	for (int i = 1; i < E; i++) {
		N0[i] = float(irand(20,65)); /*Se utilizan números aleatorios para asignar valores a los números de ocupación promedio.*/
		par0 += N0[i];
		ene0 += ENE[i] * N0[i];
	}
	cout << par0 << "\t" << ene0 << endl; /*Se imprimen en pantalla las variables par0 y ene0.*/

	/*A continuación se procederá a resolver E ecuaciones diferenciales acopladas.*/

	/*Las variables que se utilizan son las siguientes:
	a y b son dos sumas que conforman la función S, esto es S = g * (a + b). Cada una de las dos variables
	está sometida a una función delta de dirac, de ahí que en el siguiente conjunto de for's se implementen dos if's.
	Para explicar s1, s2, i1, i2, k1, k2, l1 y l2 se tomará en cuenta que al tratarse de un problema en dos dimensiones, 
	las integrales que representan la interacción deben ser escritas como sigue
	MTI(i,j,k,l) = MTI(i1,j1,k1,l1) * MTI(i2,j2,k2,l2) y entonces como cada cuarteto (i,j,k,l) corresponde a estados posibles (que son 136),
	las variables s1, s2, i1, i2, k1, k2, l1 y l2 pueden leerse de los arreglos I1 e I2. */

	float a, b;
	int s1, s2, i1, i2, k1, k2, l1, l2;
	float par1, ene1; /*par1 es el número de partículas y ene1 es la energía total del sistema en el siguiente paso temporal.*/
	float contemp = 0.0; /*Es la variable tiempo en la simulación.*/

	int contador = 1;

	for (int t = 0; t < iteraciones; t++) {
		contemp += dt; 
		for (int s = 0; s < E; s++) {
			a = 0.0;
			b = 0.0;
			s1 = I1[s], s2 = I2[s];
			for (int i = 0; i < E; i++) {
				for (int k = 0; k < E; k++) {
					for (int l = 0; l < E; l++) {
						i1 = I1[i]; i2 = I2[i];
						k1 = I1[k]; k2 = I2[k];
						l1 = I1[l]; l2 = I2[l];
						if (ENE[i]+ENE[s]-ENE[k]-ENE[l] == 0) { /*Primer función delta de dirac que garantiza conservación de energía en cada colisión..*/
							a += 2.0*MTI[i1*N*N*N+k1*N*N+l1*N+s1]*MTI[i2*N*N*N+k2*N*N+l2*N+s2]*MTI[i1*N*N*N+k1*N*N+l1*N+s1]*MTI[i2*N*N*N+k2*N*N+l2*N+s2]*(N0[s]*N0[i]*(1.0+N0[k]+N0[l])-N0[l]*N0[k]*(1.0+N0[s]+N0[i]));
						}
						if (ENE[s]-ENE[k] == 0) { /*Segunda función delta de dirac.*/
							b += 4.0*MTI[i1*N*N*N+i1*N*N+k1*N+s1]*MTI[i2*N*N*N+i2*N*N+k2*N+s2]*MTI[k1*N*N*N+l1*N*N+l1*N+s1]*MTI[k2*N*N*N+l2*N*N+l2*N+s2]*(N0[i]*N0[l]*(N0[s]-N0[k]));
						}
					}
				}
			}
			S[s] = g * (a + b);
		}
		/*El siguiente for es un paso de Euler para cada ecuación diferencial*/
		for (int j = 0; j < E; j++) {
			N1[j] = N0[j] - S[j] * dt;
			N0[j] = N1[j];
		}
		/*En este punto se guardan los valores de los números de ocupación 
		en el archivo "datos.txt" para algún estado de cada nivel energético. Recuérdese
		que existe degeneración.*/
		salida << contemp << "\t" << N1[0] << "\t" << N1[2] << "\t" << N1[4] << 
		"\t" << N1[7] << "\t" << N1[15] << "\t" << N1[20] << 
		"\t" << N1[27] << "\t" << N1[35] << "\t" << N1[44] << 
		"\t" << N1[54] << "\t" << N1[65] << "\t" << N1[77] << 
		"\t" << N1[90] << "\t" << N1[104] << "\t" << N1[119] << 
		"\t" << N1[135] << "\t" <<endl;
		par1 = 0.0; 
		ene1 = 0.0;
		for (int j = 0; j < E; j++) {
			par1 += N1[j];
			ene1 += ENE[j] * N1[j];
		}
		cout << par1 << "\t" << ene1 << endl; /*Aquí se imprimen en pantalla ene1 y par1 para verificar la conservación de la energía y del número de partículas.*/
		cout << contador << endl; /*Se imprime también un contador para saber el número de iteraciones que van.*/
		contador++;
	}

	/*Finalmente se limpia la memoria ocupada por cada uno de los arreglos definidos.*/

	delete[] TI;
	delete[] MTI;
	delete[] N0;
	delete[] N1;
	delete[] I1;
	delete[] I2;
	delete[] ENE;

	return 0;

}
