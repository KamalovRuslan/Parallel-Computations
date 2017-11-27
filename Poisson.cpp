#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

#include <cmath>
#include <limits>
#include <iomanip>

#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif


#include "Poisson.h"


using std::cout;
using std::fstream;
using std::endl;
using std::setw;

using std::ceil;
using std::max;
using std::abs;

using std::swap;
using std::memset;

// ==================================================================================================================================================
//                                                                                                                       DHP_PE_RA_FDM::DHP_PE_RA_FDM
// ==================================================================================================================================================
Process_params::Process_params ():
x_proc_number(0), y_proc_number(0),
x_cells_number(0), y_cells_number(0),
x_cell_pos(0), y_cell_pos(0),
up(false), down(false),
left(false), right(false){
}


// ==================================================================================================================================================
//                                                                                                                       DHP_PE_RA_FDM::DHP_PE_RA_FDM
// ==================================================================================================================================================
Process_params::Process_params(const Process_info& Process_info, const int grid_x,
  const int grid_y, const int x_proc_num, const int y_proc_num){

    x_proc_number = x_proc_num;
    y_proc_number = y_proc_num;

    int x_cells_per_proc = (grid_x + 1) / x_proc_number;
    int x_redundant_cells_num = (grid_x + 1) % x_proc_number;
    int x_normal_tasks_num = x_proc_number - x_redundant_cells_num;

    if (Process_info.rank % x_proc_number < x_normal_tasks_num) {
        x_cells_num = x_cells_per_proc;
        x_cell_pos = Process_info.rank % x_proc_number * x_cells_per_proc;
    } else {
        x_cells_num = x_cells_per_proc + 1;
        x_cell_pos = Process_info.rank % x_proc_number * x_cells_per_proc +
                    (Process_info.rank % x_proc_num - x_normal_tasks_num);
    }

    int y_cells_per_proc = (grid_size_y + 1) / y_proc_number;
    int y_redundant_cells_num = (grid_size_y +1) % y_proc_number;
    int y_normal_tasks_num = y_proc_num - y_redundant_cells_num;

    if (Process_info.rank / x_proc_num < y_normal_tasks_num) {
        y_cells_num = y_cells_per_proc;
        y_cell_pos = Process_info.rank / x_proc_num * y_cells_per_proc;
    } else {
        y_cells_num = y_cells_per_proc + 1;
        y_cell_pos = Process_info.rank / x_proc_number * y_cells_per_proc +
                    (Process_info.rank / x_proc_number - y_normal_tasks_num);
    }

    up = Process_info.rank < x_proc_number;
    down = Process_info.rank >= x_proc_number * (y_proc_number - 1);
    left = Process_info.rank % x_proc_number == 0;
    right = Process_info.rank % x_proc_number == x_proc_num - 1;
}


Poisson_solver::Poisson_solver(const double x1, const double y1,
               const double x2, const double y2,
               const int grid_x_, const int grid_y_,
               const double eps_, const int step_iterations_ = 1):
x_1(x1), y_1(y1), x_2(x2), y_2(y2),
dx((x2 - x1) / grid_x), dy((y2 - y1) / grid_y),
grid_x(grid_x_), grid_y(grid_y_),
eps(eps_), num_iterations(0), step_iterations(step_iterations_),
p(NULL), p_prev(NULL),
send_lr(NULL), recv_lr(NULL),
send_rl(NULL), recv_rl(NULL),
send_ud(NULL), recv_ud(NULL),
send_du(NULL), recv_du(NULL),
receive(NULL), send(NULL){
  send = new MPI_Request [4];
  receive = new MPI_Request [4];
}

Poisson_solver::~Poisson_solver(){
  p != NULL ? delete [] p; p = NULL :
  p_prev != NULL ? delete [] p_prev; p_prev = NULL :

  send_lr != NULL ? delete [] send_lr; send_lr = NULL :
  recv_lr != NULL ? delete [] recv_lr; recv_lr = NULL :

  send_rl != NULL ? delete [] send_rl; send_rl = NULL :
  recv_rl != NULL ? delete [] recv_rl; recv_rl = NULL :

  send_ud != NULL ? delete [] send_ud; send_ud = NULL :
  recv_ud != NULL ? delete [] recv_ud; recv_ud = NULL :

  send_du != NULL ? delete [] send_du; send_du = NULL :
  recv_du != NULL ? delete [] recv_du; recv_du = NULL :

  send != NULL ? delete [] send; send = NULL :
  receive != NULL ? delete [] receive; receive = NULL :
}

MPI_Comm Poisson_solver::Build_MPI_communicator(const Process_info& process_in,
                                const int x_proc_num, const int y_proc_num){
  MPI_Comm comm;
  int color = 1;
  int id = process_in.rank
  if (id < x_proc_num * y_proc_num){
      MPI_Comm_split(MPI_COMM_WORLD, color, id, &comm);
  } else {
      MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, id, &comm);
  }
  return comm;
}

Poisson_solver::Solve(const Process_info& process_in,
           const int x_proc_num, const int y_proc_num){
  MPI_Comm comm = Build_MPI_communicator(process_in, x_proc_num, y_proc_num);
  Process_info = Precoess_info(comm)
  Process_params = Process_params(Process_info, grid_x, grid_y, x_proc_num, y_proc_num);

 int shape = ProcParams.x_cells_num * ProcParams.y_cells_num;

 if (p != NULL){ delete [] p; p = NULL; }
 if (p_prev != NULL){ delete [] p_prev; p_prev = NULL; }

 p = new double [shape]();
 p_prev = new double [shape]();
 double* r =  new double [shape]();
 double*  g = new double [shape]();
 double* delta_p = new double [shape]();
 double* delta_r = new double [shape]();;
 double* delta_g = new double [shape]();

 double dg_dot_g = 1;
 double dr_dot_g = 1;
 double r_dot_g = 1;
 double alpha = 0;
 double tau = 0;

 Init_p();

 num_iterations = 0;

 do{
   Compute_diff_scheme(delta_p, p_prev);
   Compute_r(r, delta_p);

 } while();

}
