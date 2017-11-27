#include <mpi.h>

struct Process_info{

  int rank, size;
  MPI_Comm communicator;

  public:

  Process_info(MPI_Comm comm = MPI_COMM_WORLD){
    communicator = comm;
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);
  }
};


struct Process_params{

  int x_proc_number, y_proc_number;
  int x_cells_number, x_cell_pos;
  int y_cells_number, y_cell_pos;

  bool up, down, left, right;

  Process_params();
  Process_params(const Process_info& Process_info, const int grid_x,
    const int grid_y, const int x_proc_num, const int y_proc_num);
};

class Poisson_solver{

public:
  Poisson_solver(const double x1, const double y1,
                 const double x2, const double y2,
                 const int grid_x_, const int grid_y_,
                 const double eps_, const int step_iterations_ = 1);
  virtual ~Poisson_solver();

  void Solve(const Process_info& process_in,
             const int x_proc_num, const int y_proc_num);

  double* get_solution () const { return p; }
  int get_iterations () const { return num_iterations; }
  Process_info get_process_info () const { return Process_info; }
  Process_params get_process_params () const { return Process_params; }

  const double x_1, y_1;
  const double x_2, y_2;

  const double dx, dy;

  const int grid_x, grid_y;
  const double eps;

  int num_iterations;
  int step_iterations;

protected:
  virtual double F(const double x, const double y) const = 0;
  virtual double phi(const double x, const double y) const = 0;
  virtual bool Stop(const double* const f1, const double* const f2);

  MPI_Comm Build_MPI_communicator(const Process_info& process_in,
                                  const int x_proc_num, const int y_proc_num) const;
private:

  void Compute_diff_scheme(double* const delta_f, const double* const f);
  double Scalar_dot(const double* const f1, const double* const f2);

  void Compute_r(double* const r, const double* const delta_p) const;
  void Compute_g(double* const g, const double* const r, const double alpha) const;
  void Compute_p(const double tau, const double* const g);

  Process_info process_info;
  Process_params process_params;
  int step_iterations;
  int num_iterations;

  double* p;
  double* p_prev;

  double* send_lr; double* recv_lr;
  double* send_rl; double* recv_rl;
  double* send_ud; double* recv_ud;
  double* send_du; double* recv_du;

  MPI_Request* receive;
  MPI_Request* send;
};
