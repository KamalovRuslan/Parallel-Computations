class Poisson{
  private:
      int grid_x, grid_y;
      double* grid;

  public:
      Poisson();
      Poisson(const int x, const int y);
      ~Poisson()
      double& operator ()(const int x, const int y);
	  Poisson& operator = (const Poisson& tmp);

      int size_x() const;
      int size_y() const;
};

class Solver{
  private:
      const int dimension;
      const double lx;
      const double rx;
      const double l_y;
      const double r_y;
      const double delta;
      const double eps;
      const double delta2;

      double F(const double x, const double y) const;
      double phi(const double x, const double y) const;
      double x(const int index, const int shift) const;
      double y(const int index, const int shift) const;

      double DiffScheme(const Poisson& const p, const int i, const int j) const;
      double ScalarDot(const Poisson& const p, const Poisson& const q) const;
      float  ProcessDot(float var, const int rank, const int size) const;
      float  ProcessMax(float var, const int rank, const int size) const;
      void   ProcessConform(Poisson& p, const int rank, const int blocks_x,
                                     const int blocks_y);

	  void   Solve(int argc, char** argv);
  public:
	  Solver();
	  Solver(const int dimension,
	  		 const int lx, const int rx,
		 	 const int ly, const int ry,
		 	 const double eps,
		 	 const double delta);
};
