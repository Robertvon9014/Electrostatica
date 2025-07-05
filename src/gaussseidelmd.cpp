#include <iostream>
#include <mpi.h>
#include <cmath>
#include <vector>
#include <tuple>
#include <sys/time.h>

// Temporizador preciso
double seconds(){
  struct timeval tmp;
  gettimeofday(&tmp, nullptr);
  return tmp.tv_sec + tmp.tv_usec * 1e-6;
}

// Función principal paralela usando memoria distribuida
std::tuple<int, double> gaussseidel_mpi(int L, int M, double V_p, double V_n, double tolerance){
  using Matrix = std::vector<std::vector<double>>;
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (M <= 10) throw std::runtime_error("M > 10");

  // Determinar cuántas filas maneja cada proceso
  int local_rows = M / size;
  int rem = M % size;
  int start = rank * (local_rows) + std::min(rank, rem) + 1;
  local_rows += (rank < rem ? 1 : 0);
  int end = start + local_rows - 1; // inclusive

  // Creamos matriz local con halos (2 filas extra)
  int cols = M + 1;
  Matrix phi(local_rows + 2, std::vector<double>(cols, 0.0));
  Matrix phi_copy = phi;

  // Calculamos posición de barras globales
  // ---- usamos static_cast<int> para convertir de double a int de forma segura
  int fil_start = static_cast<int>(M * 0.2);        // 2 cm desde arriba
  int bar_len = static_cast<int>(M * 0.6);          // longitud de la barra 6 cm
  int fil_end = fil_start + bar_len;

  int col_plus = static_cast<int>(M * 0.2);         // Voltaje positivo a 2 cm del borde izquierdo
  int col_neg = col_plus + bar_len;                 // Voltaje negativo a 2 cm del borde derecho

  // Inicialización: cada fila (global i = start..end)
  for (int gi = start; gi <= end; ++gi){
    if (gi >= fil_start && gi < fil_end){
      phi[gi - start + 1][col_plus] = V_p;
      phi[gi - start + 1][col_neg]  = V_n;
    }
  }

  double delta = 1.0;
  int its = 0;

  while (delta > tolerance){
    its += 1;

    // -- Intercambio halo superior/inferior
    MPI_Request active_reqs[4];
    int num_reqs = 0;

    if(rank > 0){
      MPI_Irecv(&phi[0][0], cols, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &active_reqs[num_reqs++]);
      MPI_Isend(&phi[1][0], cols, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &active_reqs[num_reqs++]);
    }
    if(rank < (size - 1)){
      MPI_Irecv(&phi[local_rows + 1][0], cols, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &active_reqs[num_reqs++]);
      MPI_Isend(&phi[local_rows][0], cols, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &active_reqs[num_reqs++]);
    }

    MPI_Waitall(num_reqs, active_reqs, MPI_STATUSES_IGNORE);

    // --- Actualizamos celdas rojas
    for (int i = 1; i <= local_rows; ++i){
      for (int j = 1; j < M; ++j){
        int gi = start + i -1; // gi = global i
        if ((gi + j) % 2 == 0){
          if ((fil_start <= gi && gi <= fil_end) && ( j == col_plus || j == col_neg)){
            continue;
          }
          else{
            phi[i][j] = 0.25 * (phi[i + 1][j] + phi[i - 1][j] + phi[i][j  + 1] + phi[i][j - 1]);
          }
        }
      }
    }

    // --- Actualizamos las celdas negras
    for (int i = 1; i <= local_rows; ++i){
      for (int j = 1; j < M; ++j){
        int gi = start + i - 1;
        if ((gi + j) % 2 == 1){
          if ((fil_start <= gi && gi <= fil_end) && (j == col_plus || j == col_neg)){
            continue;
          }
          else{
            phi[i][j] = 0.25 * (phi[i + 1][j] + phi[i - 1 ][j] + phi[i][j + 1] + phi[i][j - 1]);
          }
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD); // Sincronizamos los colores

    // --- Cálculo del delta local
    delta = 0.0;
    for (int i = 1; i <= local_rows; ++i){
      for (int j = 1; j <= M; ++j){
        double diferencia = std::abs(phi[i][j] - phi_copy[i][j]);
        if (diferencia > delta){
          delta = diferencia;
        }
      }
    }

    // --- Reducción global del delta
    MPI_Allreduce(MPI_IN_PLACE, &delta, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // --- Copiamos los valores de la nueva iteración
    phi_copy = phi;
  }

  // --- Imprimir el proceso que posee la fila 50
  int gi_target = 50;
  if (start <= gi_target && gi_target <= end){
    int li = gi_target - start + 1; 
    std::cout << "phi[50][25] = " << phi[li][25] << std::endl;
  }

  return std::make_tuple(its, delta);
}

int main(int argc, char**argv){
  MPI_Init(&argc,&argv);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  double t0 = seconds();
  auto [its, delta] = gaussseidel_mpi(10, 100, 1.0, -1.0, 1e-5);
  double t1 = seconds();

  if(rank == 0){
    std::cout<<"Iteraciones: "<<its<<" error: "<<delta<<"\n";
    std::cout<<"Tiempo MPI: "<<(t1-t0)<<" s\n";
  }

  MPI_Finalize();
  return 0;
}
