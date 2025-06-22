#include <iostream>
#include <cmath>
#include <vector>
#include <tuple>

std::tuple<int, double> gaussseidel(int L, int M, double V_p, double V_n, double tolerance){
  // creamos un arreglo de 2-dimensiones usando std::vector
  using Matrix = std::vector<std::vector<double>>;
  Matrix phi(M + 1, std::vector<double>(M + 1, 0.0));

  // ---- Calcular la reposición dependiendo del valor de M
  // ---- Posiciones expresadas en índices de grilla
  int fil_start = (2 * M) / L; // 2 cm desde arriba
  int bar_len = (6 * M) / L;   // longitud de la barra 6 cm
  int fil_end = fil_start + bar_len;

  int col_plus = (2 * M) / L;  // Voltaje positivo a 2 cm del borde izquierdo
  int col_neg = col_plus + bar_len; // voltaje negativo a 2 cm del borde derecho 
  
  // Colocamos las condiciones iniciales 
  for (int i = fil_start; i < fil_end; ++i){
    phi[i][col_plus] = V_p;
    phi[i][col_neg] = V_n;
  }

  // Guardamos una copia
  Matrix phi_copy = phi;

  // Iteración Gauss-Seidel
  double delta = 1.0;
  int its = 0;
  while (delta > tolerance){
    // calculamos la iteración
    its += 1;
    for (int i = 1; i < M; ++i){
      for (int j = 1; j < M; ++j){
        // Condicion de frontera
        // Valores de frontera no se modifican
        if ((fil_start <= i && i <= fil_end) && (j == col_plus || j == col_neg)){
          continue;
        }
        else{
          phi[i][j] = 0.25 * (phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j -1]);
        }
      }
    }

    // Calculamos la diferencia máxima con respecto a los valores anteriores
    delta = 0.0;
    for (int i = 0; i <= M; ++i){
      for (int j = 0; j <= M; ++j){
        double diferencia = std::abs(phi[i][j] - phi_copy[i][j]);
        if (diferencia > delta){
          delta = diferencia;
        }
      }
    }

    // copiamos los valores de la nueva iteracion
    phi_copy = phi;
  }

  return std::make_tuple(its, delta);
}

int main(){
  int iteraciones;
  double error;
  std::tie(iteraciones, error) = gaussseidel(10, 100, 1.0, -1.0, 1e-5);

  std::cout << " Convergencia alcanzada en " << iteraciones << " iteraciones con error " << error << std::endl;
  return 0;
}
