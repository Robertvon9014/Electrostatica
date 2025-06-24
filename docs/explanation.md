# Explicación de los Métodos Numéricos

## Problema físico

Se resuelve la ecuación de Laplace bidimensional para el potencial eléctrico \(\phi(x,y)\):

\[
\frac{\partial^2 \phi}{\partial x^2} + \frac{\partial^2 \phi}{\partial y^2} = 0
\]

en una placa cuadrada de 10 cm × 10 cm, con dos barras internas que simulan las placas de un capacitor con potenciales \(V_p\) y \(V_n\).

---

## Métodos implementados

### 1. Método de relajación de Jacobi

Este método itera actualizando el potencial en cada punto de la malla usando el promedio de sus vecinos en la iteración anterior, manteniendo fijas las condiciones de frontera.

- Es sencillo pero convergencia lenta.
- Requiere dos matrices para mantener valores antiguos y nuevos.

---

### 2. Método de sobre-relajación de Jacobi (Jacobi modificado)

Se aplica un factor de sobre-relajación \(\omega\) para acelerar la convergencia, mezclando el valor nuevo calculado con el anterior:

\[
\phi^{new} = (1+\omega) \times \text{valor calculado} - \omega \times \phi^{old}
\]

- Acelera la convergencia respecto al Jacobi estándar.
- Requiere ajuste del parámetro \(\omega\) para optimizar rendimiento.

---

### 3. Método de Gauss-Seidel

En este método, la actualización del potencial en cada punto usa los valores más recientes disponibles (actualizados durante la misma iteración), mejorando la velocidad de convergencia frente a Jacobi.

- Solo se necesita una matriz.
- Mejor convergencia en menos iteraciones.

---

### 4. Método de Gauss-Seidel en C++

Implementación equivalente al método de Gauss-Seidel en Python, pero usando C++ con `std::vector` y `std::tuple`.

- Permite una posible extensión para paralelización.
- Mejora en eficiencia y manejo de memoria para casos grandes.

---

## Condiciones de frontera y configuración de la grilla

- La placa se discretiza en una malla de \((M+1) \times (M+1)\) puntos.
- Dos barras verticales se colocan a 2 cm de cada borde lateral, con longitudes de 6 cm, y potenciales fijos \(V_p\) y \(V_n\).
- El resto de bordes y posiciones iniciales son 0.

---

## Convergencia y criterio de parada

La iteración se detiene cuando la diferencia máxima absoluta entre iteraciones consecutivas es menor que una tolerancia dada.

---

## Comparación y resultados

- El método de Gauss-Seidel converge más rápido que Jacobi.
- La sobre-relajación mejora aún más la velocidad si se escoge un \(\omega\) adecuado (por ejemplo, 0.9).
- La implementación en C++ es más eficiente y es base para paralelización futura.


