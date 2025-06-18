"""
Este paquete proporciona funciones para resolver sistemas de ecuaciones lineales
utilizando métodos iterativos clásicos sin depender de un número máximo de iteraciones,
sino de una condición de convergencia basada en la tolerancia.

Módulos exportados por este paquete:

- `jaccobi`: Implementa el método de Jacobi tradicional para resolver sistemas Ax = b.
- `jaccobi_modified`: Variante del método de Jacobi que incorpora un factor de relajación (omega).
- `gauss_seidel`: Implementa el método de Gauss-Seidel, que actualiza los valores en cada iteración.

Cada módulo trabaja con matrices NumPy y detiene la iteración cuando se cumple
una condición de tolerancia establecida por el usuario.
"""
