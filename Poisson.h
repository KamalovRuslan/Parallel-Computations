class Poisson{
  private:
      int grid_x, grid_y;
      double* grid;

  public:
      Poisson();
      Poisson(const int x, const int y);
      ~Poisson()
      double& operator ()(const int x, const int y);

      int size_x() const;
      int size_y() const;
};

class Solver{
  private:
      const int dimension = 1000;
      const double lx     = 0.0; 
      const double rx     = 2.0;
      const double l_y    = 0.0;
      const double r_y    = 2.0;
      const double delta  = (rx - lx) / dim;
      const double eps    = 1e-4;
      const double delta2 = delta * delta;

      double F(const double x, const double y) const;
      double phi(const double x, const double y) const;
      double x(const int index, const int shift) const;
      double y(const int index, const int shift) const;

      double DiffScheme(const grid& const p, const int i, const int j) const;
      double ScalarDot(const grid& const p, const grid& const q) const;
      float  ProcessDot(float var, const int rank, const int size) const;
      void   ProcessConform(grid& p, const int rank, const int blocks_x,
                                     const int blocks_y);
};
