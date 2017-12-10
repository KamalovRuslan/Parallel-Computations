#ifndef POISSON_H
#define POISSON_H


class Poisson{
  private:
      int grid_x, grid_y;
      double* grid;

  public:
      Poisson();
      Poisson(const int x, const int y);
      ~Poisson();
      double& operator ()(const int x, const int y);
	  Poisson& operator = (const Poisson& tmp);
	  double max() const;

      int size_x() const;
      int size_y() const;
};

class Solver{
  private:
      int dimension;
      double lx; double rx;
      double ly; double ry;
      double delta; double delta2;
	  double eps;

      double F(const double x, const double y) const;
      double phi(const double x, const double y) const;
      double x(const int index, const int shift) const;
      double y(const int index, const int shift) const;

      double DiffScheme(Poisson& p, int i, int j) const;
      double ScalarDot(Poisson& p, Poisson& q) const;
      float  ProcessDot(float var, const int rank, const int size) const;
      float  ProcessMax(float var, const int rank, const int size) const;
      void   ProcessConform(Poisson& p, const int rank, const int blocks_x,
                                     const int blocks_y);

  public:
	  Solver();
	  Solver(const int dimension,
	  		 const int lx, const int rx,
		 	 const int ly, const int ry,
		 	 const double eps,
		 	 const double delta);

	  void   Solve(int argc, char** argv);
};
#endif //POISSON_H
