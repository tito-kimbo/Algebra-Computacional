# Entrega

En este documento especificamos los algoritmos hay en el repositorio

## Algoritmos básicos


* Algoritmo de euclides para R un D.E cualquiera 

   <b>Teoría</b>: página 20, Algoritmo 2.1.4 de [ICA] 
   
   <b>Algoritmo</b>: En la carpeta `algorithms`, en el fichero `divisibility.py`, es la función `gcd`


* Algoritmo de euclides extendido

   <b>Teoría</b>: página 60, Teorema 4.3 de [CINT]

   <b>Algoritmo</b>: En la carpeta `algorithms`, en el fichero `divisibility.py`, es la función `eea`


* Algoritmo para calcular el teorema chino del resto (i.e. calcular el inverso)

   En la carpeta `algorithms`, en el fichero `chinese_remainder.py`, es la función `chinese_remainder`


* Mcd en un D.F.U.

   <b>Teoría</b>: página 23, Algoritmo 2.1.10 de [ICA] 
   
   <b>Algoritmo</b>: En la carpeta `algorithms`, en el fichero `gcd_dfu.py`, es la función `gcd_dfu`


* Inverso de un elemento en un cuerpo finito.

   En la carpeta `algorithms`, en el fichero `chinese_remainder.py`, es la función `modinv`



* Test de irreducibilidad de un poliomio en Fq[x]

   En la carpeta `algorithms`, en el fichero `primality`, es la función `rabin_test`

* Logaritmo discreto en cuerpos Fq[x]/(f(x))

   En la carpeta `algorithms`, en el fichero `discrete_log.py`, es la función `discrete_log`

* Algoritmo de factorización de un polinomio en cuerpo finito parte 1, 2 y 3.

   En la carpeta `algorithms`, en el fichero `factorization`, functiones `squarefree_decomposition`, `distinct_degree_factorization`, `berlekamp_splitting`.


* Algoritmo de factorizacion de Berlekamp de polinomios sobre cuerpos finitos 

   En la carpeta `algorithms`, en el fichero `factorization`, funciones `berlekamp_factorization` y `berlekamp_cantor_zassenhaus`

* Algorimos de factorizacion en Z[x]

   En la carpeta `algorithms`, en el fichero `zx_factorization.py`, es la función `zx_factorization`

* Algorimo de primalidad de AKS

   En la carpeta `algorithms`, en el fichero `primality.py`, es la función `is_prime_aks`



## Algoritmos adicionales


* Algorimo de primalidad de Miller-Rabin

   En la carpeta `algorithms`, en el fichero `primality.py`, es la función `is_prime_miller_rabin`


* Algorithmo de extracción de raíces en cuerpos finitos de Adleman-Manders-Miller

   En la carpeta `algorithms`, en el fichero `root.py`, es la función `modsqrt`
   
* Algoritmo para calcular el mcd en los enteros empleando el método del binary gcd

   <b>Teoría</b>: página 57, Ejercicio 4.1 de [CINT]

   <b>Algoritmo</b>: En la carpeta `algorithms`, en el fichero `divisibility.py`, es la función `binary_gcd`
   
## Bibliografía
[ICA] Frühbis-Krüger, A., & Lossen, C. (s. f.). Introduction to Computer Algebra (Solving Systems of Polynomial Equations).

[CINT] Shoup, V. (s. f.). A Computational Introduction to Number Theory and Algebra.
