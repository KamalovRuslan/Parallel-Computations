#include "Poisson"

Poisson::Poisson(){
      grid = NULL;
      grid_x = 0; grid_y = 0;
}

Poisson::Poisson(const int x, const int y){
      grid_x = x; grid_y = y;
      grid = new double[x * y];
}

Poisson::~Poisson(){
      if (grid != NULL){
            delete [] grid;
      }
}

double& Poisson::operator ()(int x, int y){
      return grid[x * grid_y + y];
}

int Poisson::size_x(){
      return grid_x;
}

int Poisson::size_y(){
      return grid_y;
}

double Solver::x(const int index, const int shift){
      return lx + (index + shift) * delta;
}

double Solver::y(const int index, const int shift){
      return ly + (index + shift) * delta;
}

double Solver::DiffScheme(const grid& const p, const int i, const int j){
      return ((2 * p(i, j) - p(i - 1, j) - p(i + 1, j))  +
              (2 * p(i, j) - p(i, j - 1) - p(i, j + 1))) / delta2;
}

double Solver::ScalarDot(const grid& const p, const grid& const q){
      double scalar_dot = 0;
      int size_x = p.size_x();
      int size_y = p.size_y();
      for (int i = 1; i < size_x - 1; i++) {
          for (int j = 1; j < size_y - 1; j++) {
              scalar_dot += delta2 * p(i, j) * q(i, j);
          }
      }
    return scalar_dot;
}

float Solver::ProcessDot(float var, const int rank, const int size){
      float *processes_sum;
      if (rank == 0){
          processes_sum = new float[size];
      }

      MPI_Gather(&var, 1, MPI_FLOAT, sums, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

      float sum = 0.0f;
      if (rank == 0) {
          for (int i = 0; i < size; i++)
              sum += processes_sum[i];
          delete [] processes_sum;
      }

      MPI_Bcast(&sum, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
      return sum;
}
