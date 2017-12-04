#include <mpi.h>
#include <omp.h>

struct grid_info{
  int left, right, top, bottom;
};

struct buffer_set{
  double *start;
  int    left, right, top, bottom;
  int    d_left, d_right, d_top, d_bottom;
};

struct comm_data_t{
  MPI_Comm           grid_comm;
  struct grid_info   info;
  struct buffer_set  send_buffers, recv_buffers;
  int                rank;
  int                coord[2]
};

class PoissonSolver{
  public:

    comm_data_t GetCommData() const {return comm_data;}

    PoissonSolver(const double x1, const double x2,
                  const double y1, const double y2,
                  const int num_iteration, double epsilon);
    virtual ~PoissonSolver();

    void Solve(const comm_data_t& comm_data);

  protected:

    virtual double F(const double x, const double y) const = 0;
    virtual double phi(const double x, const double y) const = 0;

  private:
    comm_data_t comm_data;
    double hx, hy, hxhy, hx2, hy2;
    int iteration_counter, num_iteration;
    double epsilon;

    double* p;
    double* p_prev;

    double* send_message_lr;
    double* send_message_rl;
    double* send_message_td;
    double* send_message_bu;
    double* recv_message_lr;
    double* recv_message_rl;
    double* recv_message_td;
    double* recv_message_bu;
    MPI_Request* recv;
    MPI_Request* send;

    void BuildComm(comm_data_t& comm_data);

    void ComputeDiffScheme(double* const delta_f, const double* const f);
    double ScalarDot(const double* const f1, const double* const f2);
    void Compute_r(double* const r, const double* const delta_p) const;
    void Compute_g(double* const g, const double* const r, const double alpha) const;
    void Compute_p(const double tau, const double* const g);

}
