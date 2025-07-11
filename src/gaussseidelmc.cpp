#include <iostream>
#include <cmath>
#include <vector>
#include <tuple>
#include <sys/time.h>
#include <omp.h>
#include <iomanip>

/**
 * @brief Devuelve el tiempo actual en segundos desde Epoch.
 *
 * Esta función se usa para calcular el tiempo de ejecución de secciones del código.
 *
 * @return double Tiempo en segundos como número decimal de alta precisión.
 */

double seconds(){
  struct timeval tmp;
  double sec;
  gettimeofday( &tmp, (struct timezone *)0 );
  sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;

  return sec;
}

/**
 * @brief Resuelve el potencial eléctrico en una placa cuadrada con método Gauss-Seidel usando red-black ordering y OpenMP.
 *
 * Se impone una grilla de tamaño `M x M` sobre la placa, con dos barras verticales:
 * una con voltaje positivo `V_p` y otra con voltaje negativo `V_n`. 
 * El algoritmo aplica el método de Gauss-Seidel con ordenamiento red-black en paralelo.
 *
 * @param M Número de divisiones de la grilla (debe ser mayor a 10)
 * @param V_p Voltaje positivo aplicado
 * @param V_n Voltaje negativo aplicado
 * @param tolerance Tolerancia para el criterio de convergencia
 * @param num_procs Referencia donde se guardará el número de procesos/hilos usados
 * @return std::tuple<int, double> Número de iteraciones realizadas y error final
 */

std::tuple<int, double> gaussseidel(int M, double V_p, double V_n, double tolerance, int &num_procs){
  // Creamos un arreglo de 2-dimensiones usando std::vector
  using Matrix = std::vector<std::vector<double>>;
  Matrix phi(M + 1, std::vector<double>(M + 1, 0.0));

  // ---- Validar que M sea mayor a 10
  if (M <= 10){
    throw std::runtime_error("El tamaño de la grilla (M) debe ser mayor a 10");
  }

  // ---- Calcular la reposición dependiendo del valor de M
  // ---- Posiciones expresadas en índices de grilla
  // ---- usamos static_cast<int> para convertir de double a int de forma segura
  int fil_start = static_cast<int>(M * 0.2);        // 2 cm desde arriba
  int bar_len = static_cast<int>(M * 0.6);          // longitud de la barra 6 cm
  int fil_end = fil_start + bar_len;

  int col_plus = static_cast<int>(M * 0.2);         // Voltaje positivo a 2 cm del borde izquierdo
  int col_neg = col_plus + bar_len;                 // Voltaje negativo a 2 cm del borde derecho

  // Colocamos las posiciones iniciales
  for (int i = fil_start; i < fil_end; ++i){
    phi[i][col_plus] = V_p;
    phi[i][col_neg] = V_n;
  }

  // Guardamos una copia
  Matrix phi_copy = phi;

  // Iteración de Gauss-seidel
  double delta = 1.0;
  int its = 0;
  while (delta > tolerance){
    its += 1;
    #pragma omp parallel 
    {
      #pragma omp single
      {
        num_procs = omp_get_num_threads();
      }
      /* --- Implementar el método de red-black ordering
       * Este método consiste en colorear todo el dominio 
       * de la grilla como si fuera un tablero de ajedrez:
       * en celdas rojas y negras
       * R N R N 
       * N R N R 
       * R N R N 
       * N R N R
       * La clave de este método es que cada paso iremos 
       * actualizando las celdas rojas primero (usando las
       * negras) y luego las celdas negras (usando las rojas
       * recién actualizadas). */

      // Actualizamos las celdas rojas
      #pragma omp for collapse(2)
      for (int i = 1; i < M; ++i){
        for (int j = 1; j < M; ++j){
          if ((i + j) % 2 == 0){
            if ((fil_start <= i && i <= fil_end) && (j == col_plus || j == col_neg)){
              continue;
            }
            else{
              phi[i][j] = 0.25 * (phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j - 1]);
            }
          }
        }
      }

      // Actualizamos las celdas negras
      #pragma omp for collapse(2)
      for (int i = 1; i < M; ++i){
        for (int j = 1; j < M; ++j){
          if ((i + j) % 2 == 1){
            if ((fil_start <= i && i <= fil_end) && (j == col_plus || j == col_neg)){
              continue;
            }
            else{
              phi[i][j] = 0.25 * (phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j -1]);
            }
          }
        }
      }
    }

    // Calculamos la diferencia máxima con respecto a los valores anteriores
    delta = 0.0;
    #pragma omp parallel reduction(max:delta)
    {
      #pragma omp for
      for (int i = 0; i <= M; ++i){
        for (int j = 0; j <= M; ++j){
          double diferencia = std::abs(phi[i][j] - phi_copy[i][j]);
          if (diferencia > delta){
            delta = diferencia;
          }
        }
      }
    }

    // Copiamos los valores de la nueva iteración

    phi_copy = phi;
  }

  //std::cout << "phi[50][25] = " << phi[50][25] << std::endl;
  return std::make_tuple(its, delta);
}

/**
 * @brief Función principal que ejecuta la simulación de Gauss-Seidel con OpenMP.
 *
 * Mide el tiempo de ejecución, llama a `gaussseidel` con M=500 y muestra el resultado.
 *
 * @return int Código de salida del programa (0 si termina correctamente)
 */

int main(){
  int iteraciones;
  double error;

  int num_procs;

  std:: cout.precision(10); // configuramos la salida de decimales con una precision de 10 números
  // iniciamos el temporizador
  double time_1 = seconds();

  std::tie(iteraciones, error) = gaussseidel(500, 1.0, -1.0, 1e-5, num_procs);
  std::cout << "Convergencia alcanzada en " << iteraciones << " iteraciones con error " << error << std::endl;


  // Finalizamos el temporizador
  double time_2 = seconds();

  std::cout << "Numero de procesos: " << num_procs << std::endl;
  std::cout << std::fixed << std::setprecision(4);
  std::cout << "Tiempo de ejecución: " << time_2 - time_1 << " segundos" << std::endl;

  return 0;
}
