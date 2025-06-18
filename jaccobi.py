#!/usr/bin/env python3

import numpy as np

# Esta función toma como argumentos el tamaño lineal de la grilla cuadrada,
# V_p: Voltaje positivo
# V_n: Voltaje negativo
def jacobi_relaxation(L, M, V_p, V_n, tolerance):
    # Primero creamos los arreglos 2-dimensionales de la grilla
    # Vamos a necesitar dos según la regla de Jacobi
    # Note que usamos M+1, debido a que debemos contener la condición de frontera
    # phi contiene inicialmente los valores iniciales. Vamos a utilizar ceros.
    phi = np.zeros((M + 1, M + 1), dtype=float)
    # Ahora tenemos que colocar la condición inicial.
    # Recuerde accesos de listas en np.ndarray

    # --- Calculamos la reposición dependiendo el tamaño de M
    fil_start = int((2 * M) / L) # 2 cm desde arriba
    bar_len = int((6 * M) / L)   # 6 cm longitud de la barra
    fil_end = fil_start + bar_len

    col_plus = int((2 * M) / L)  # Voltaje positivo a 2 cm del borde izquierdo
    col_neg = col_plus + bar_len # Voltaje Negativo a 2 cm del borde derecho

    # Ahora tenemos que colocar la condición inicial.
    # Recuerde accesos de listas en np.ndarray
    phi[fil_start:fil_end, col_plus] = V_p  # Barra positiva
    phi[fil_start:fil_end, col_neg] = V_n    # Barra negativa

    # phiprime se necesita para la iteración
    phiprime = np.zeros((M + 1, M + 1), dtype=float)
    # Iteración de Jacobi
    delta = 1.0
    its = 0
    while delta > tolerance:
        # Calculamos la iteración
        its += 1
        for i in range(M + 1):
            for j in range(M + 1):
                # Condición de frontera
                if j == col_plus and fil_start <= i <= fil_end or j == col_neg and fil_start <= i <= fil_end:
                    phiprime[i, j] = phi[i,j]
                elif i == 0 or i == M or j == 0 or j == M:
                    phiprime[i, j] = phi[i, j]
                # Iteración principal
                else:
                # COMPLETE AQUÍ
                    phiprime[i,j] = 0.25 * (phi[i + 1, j] + phi[i - 1, j] + phi[i, j + 1] + phi[i, j - 1])
        # Calculamos la diferencia máxima con respecto a los valores anteriores
        delta = np.max(np.abs(phi - phiprime))
        # Ahora intercambiamos los arreglos para la nueva iteración
        # El nuevo phi es el phiprime
        temp = phi
        phi = phiprime
        # El nuevo phiprime es el phi viejo
        phiprime = temp

    return phi, its, delta


jacobi_vals, iterations, error = jacobi_relaxation(10, 100, 1.0, -1.0, 1e-5)
print(f"Convergencia alcanzada en {iterations} iteraciones con error {error}")
