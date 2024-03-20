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
               const std::function<double(const std::array<double,2> &)> & d0fun_,
               const std::function<double(const std::array<double,2> &)> & d1fun_,
               const double alpha_ = 0.2,
               const unsigned int n_max_it_ = 100,
               const double tol_x_ = std::numeric_limits<double>::epsilon() * 1000.0,
               const double tol_fun_   = std::numeric_limits<double>::epsilon() * 1000.0)
    : m_fun(fun_)
    , m_d0fun(d0fun_)
    , m_d1fun(d1fun_)
    , m_alpha(alpha_)
    , m_n_max_it(n_max_it_)
    , m_tol_fun(tol_fun_)
    , m_tol_x(tol_x_)
    , m_x()
    , m_df_dx()
    , m_iter(0)
  {}

  void minimize(const std::array<double,2> &x0)
  {
    m_old_x = x0;
    int N = x0.size();
    std::array<double,2> step;

    for (size_t m_iter = 0; m_iter < N; ++m_iter)
      {
        step[0] = m_alpha* m_d0fun(m_old_x);
        step[1] = m_alpha* m_d1fun(m_old_x);
        std::transform(m_old_x.begin(), m_old_x.end(), step.begin(), m_x.begin(),
                    [=](double x, double y){ return x - y ; });
        m_res = 0;
        for(size_t i = 0; i < x0.size(); ++i) m_res += (m_x[i]-m_old_x[i])*(m_x[i]-m_old_x[i]);
        if(m_res < m_tol_x)
            break;
        if(abs(m_fun(m_old_x)-m_fun(m_x)) < m_tol_fun)
            break;
      };
  }

  std::array<double,2>  get_result() const
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
  const std::function<double(const std::array<double,2> &)> m_d0fun;
  const std::function<double(const std::array<double,2> &)> m_d1fun;

  // number of maximum iteration
  const unsigned int m_n_max_it;

  // tolerance on the residual
  const double       m_tol_fun;

  // tolerance on the increment
  const double       m_tol_x;

  double       m_alpha;

  // Variables employed by the solver
  // guess on the solution
  std::array<double,2> m_x;

  std::array<double,2> m_old_x;

  // value of the derivative at point x
  std::array<double,2> m_df_dx;

  // number of iterations
  unsigned int m_iter;

  double m_res;
};

#endif