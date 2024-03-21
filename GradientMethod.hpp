#ifndef GRADIENTMETHOD
#define GRADIENTMETHOD

#include <cmath>
#include <iostream>
#include <functional>
#include <iostream>
#include <limits>
#include "muParserXInterface.hpp" 

enum class alpha_strategies
{
  Exponencial_decay,
  Inverse_decay,
  Aproximate_line_search
};

class GradientMethod
{
public:
  GradientMethod(const std::function<double(const std::array<double,2> &)> & fun_,
               const std::function<double(const std::array<double,2> &)> & d0fun_,
               const std::function<double(const std::array<double,2> &)> & d1fun_,
               const double alpha_ = 0.2,
               const unsigned int n_max_it_ = 100,
               const double tol_x_ = std::numeric_limits<double>::epsilon() * 1000.0,
               const double tol_fun_   = std::numeric_limits<double>::epsilon() * 1000.0,
               const alpha_strategies strategy_ = alpha_strategies::Exponencial_decay)
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
    , m_strategy(strategy_)
  {}

  void minimize(const std::array<double,2> &x0)
  {
    int N = x0.size();
    m_x = x0;
    double const alpha0 = m_alpha;
    while(m_iter < m_n_max_it)
      {
        ++m_iter;
        m_old_x = m_x;
        if(m_strategy == alpha_strategies::Exponencial_decay) m_alpha = alpha0*pow(exp(-0.2),m_iter);
        if(m_strategy == alpha_strategies::Inverse_decay) m_alpha = alpha0/(1+0.2*m_iter);
        if(m_strategy == alpha_strategies::Aproximate_line_search){ 
          //we do alpha divided by 2 until Armijo rule is satisfyed
          m_alpha = alpha0;
          double const Armijo = 0.4 *alpha0*(m_d0fun(m_old_x)*m_d0fun(m_old_x)+m_d1fun(m_old_x)*m_d1fun(m_old_x));
          bool anti_loops = false;
          while(m_fun(m_old_x)-m_fun({m_old_x[0]-m_alpha*m_d0fun(m_old_x),m_old_x[1]-m_alpha*m_d1fun(m_old_x)}) < Armijo){
            m_alpha /= 2;
            if(m_alpha < 1e-5){ anti_loops = true; break;}
          } 
          if(anti_loops) break;
        }
        m_x[0] -= m_alpha*m_d0fun(m_old_x);
        m_x[1] -= m_alpha*m_d1fun(m_old_x);

        m_res = (m_x[0]-m_old_x[0])*(m_x[0]-m_old_x[0])+(m_x[1]-m_old_x[1])*(m_x[1]-m_old_x[1]);
        if(m_res < m_tol_x) break;
        if(abs(m_fun(m_old_x)-m_fun(m_x)) < m_tol_fun) break;     
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

  const alpha_strategies m_strategy;
};

#endif