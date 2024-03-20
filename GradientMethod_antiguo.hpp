#ifndef GRADIENTMETHOD
#define GRADIENTMETHOD

#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include "muParserXInterface.hpp" 

class GradientMethod
{
public:
  GradientMethod(const std::function<double(const std::array<double,2> &)> & fun_,
               const std::vector< std::function<double(const std::array<double,2>)>> & dfun_,
               const double alpha_ = 0.2,
               const unsigned int n_max_it_ = 100,
               const double tol_x_ = std::numeric_limits<double>::epsilon() * 1000.0,
               const double tol_fun_   = std::numeric_limits<double>::epsilon() * 1000.0)
    : m_fun(fun_)
    , m_dfun(dfun_)
    , m_alpha(alpha_)
    , m_n_max_it(n_max_it_)
    , m_tol_fun(tol_fun_)
    , m_tol_x(tol_x_)
    , m_x()
    , m_df_dx()
    , m_iter(0)
  {}

  void minimize(const std::vector<double> &x0)
  {
    m_old_x = x0;
    int N = x0.size();
    std::vector<double> step(N);

    for (size_t m_iter = 0; m_iter < N; ++m_iter)
      {
        for(size_t i = 0; i < N; ++i) step[i] = m_alpha* m_dfun[i](m_old_x);
        std::transform(m_old_x.begin(), m_old_x.end(), step.begin(), m_x.begin(),
                    [=](double x, double y){ return x - y ; });
        m_res = 0;
        for(size_t i = 0; i < x0.size(); ++i) m_res += (m_x[i]-m_old_x[i])*(m_x[i]-m_old_x[i]);
        if(m_res < m_tol_x)
            break;
        if(abs(m_fun(m_old_x)-m_fun(m_x)) < m_tol_fun)
            break;
      }
  }

  std::vector<double>  get_result() const
  {
    return m_x;
  };

  double get_residual() const
  {
    return m_res;
  };

  unsigned int get_iter() const
  {
    return m_iter;
  };

private:

  // Variables to initialize the Newton method
  // function f
  const std::function<double(const std::array<double,2> &)> m_fun;

  // derivative of f
  const std::vector< std::function<double(std::array<double,2>)>> m_dfun;

  // number of maximum iteration
  const unsigned int m_n_max_it;

  // tolerance on the residual
  const double       m_tol_fun;

  // tolerance on the increment
  const double       m_tol_x;

  double       m_alpha;

  // Variables employed by the solver
  // guess on the solution
  std::vector<double> m_x;

  std::vector<double> m_old_x;

  // value of the derivative at point x
  std::vector<double> m_df_dx;

  // number of iterations
  unsigned int m_iter;

  double m_res;
};

#endif