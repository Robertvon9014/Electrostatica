#!/usr/bin/env python3

import numpy as np

def gauss_seidel(L, M, V_p, V_n, tolerance):
    """
    Resuelve el potencial eléctrico en una placa cuadrada usando el método de Gauss-seidel.

    Se modela la región cuadrada de tamaño `L x L`dividida en `M x M` puntos.
    Se tienen dos barras verticales que se colocan como condiciones de frontera internas:
    Una posee un potencial `V_p` y otra con `V_n`. La función itera hasta que el cambio
    máximo entre iteraciones sea menor que la `tolerance`.

    Args:
        L (int): Tamaño físico de la placa cuadrada (dado en cm).
        M (int): Número de divisiones de la grilla (grilla de (M+1) x (M+1)).
        V_p (float): Voltaje aplicado en la barra positiva.
        V_n (float): Voltaje aplicado en la barra negativa.
        tolerance (float): Tolerancia para el criterio de convergencia.

    Returns:
        phi (ndarray): Matriz de 2 dimensiones con potenciales verticales dados por dos barras.
        its (int): Número de iteraciones realizadas.
        delta (float): Error máximo alcanzado en la última iteración.

    Examples:
        >>> phi, its, error = jacobi_relaxation(10, 100, 1.0, -1.0, 1e-5)
        >>> print(f"Convergencia alcanzada en {itertions} iteraciones con error de {error:.2e}")
        Convergencia alcanzada en 1125 iteraciones con error de 9.97e-06
    """
    # Primero creamos el arreglo 2-dimensionales de la grilla
    # Note que usamos M+1, debido a que debemos contener la condición de frontera
    # phi contiene inicialmente los valores iniciales. Vamos a utilizar ceros.
    phi = np.zeros((M + 1, M + 1), dtype=float)

    # Validamos que M sea mayor que 10
    if M <= 10:
        raise ValueError("El tamaño de la grila (M) debe ser mayor a 10")

    # --- Calculamos la reposición dependiendo del valor de M
    fil_start = int(M * 0.2) # 2 cm desde arriba
    vol_len = int(M * 0.6)   # 6 cm longitud de la barra
    fil_end = fil_start + vol_len

    col_plus = int(M * 0.2)  # voltaje positivo a 2 cm del borde izquierdo
    col_neg = col_plus + vol_len # Voltaje negativo a 2 cm del borde derecho

    # Ahora tenemos que colocar la condición inicial.
    # Recuerde accesos de listas en np.ndarray
    phi[fil_start:fil_end, col_plus] = V_p
    phi[fil_start:fil_end, col_neg] = V_n
    # Vamos a guardar una copia para evaluar el error
    phi_copy = phi.copy()
    # Iteración de Gauss-Seidel
    delta = 1.0
    its = 0
    while delta > tolerance:
        # Calculamos la iteración
        its += 1
        for i in range(1, M):
            for j in range(1, M):
                # Condición de frontera
                # En este caso, en la frontera los valores no se modifican
                if fil_start <= i <= fil_end and (j == col_plus or j == col_neg):
                    continue
                # Iteración principal
                else:
                    # COMPLETE AQUÍ
                    phi[i,j] = 0.25 * (phi[i + 1, j] + phi[i - 1, j] + phi[i, j + 1] + phi[i, j - 1] )
        # Calculamos la diferencia máxima con respecto a los valores anteriores
        delta = np.max(np.abs(phi - phi_copy))
        # Copiamos los valores de la nueva iteración
        phi_copy = phi.copy()

    return phi, its, delta

gaussSeidel_vals, iterations, error = gauss_seidel(10, 100, 1.0, -1.0, 1e-5)
print(f"Convergencia alcanzada en {iterations} iteraciones con error de {error:.2e}")
print("Valores grilla phi[50][25] ", gaussSeidel_vals[50][25])
