#include "muParserXInterface.hpp" 
#include "GradientMethod.hpp"
#include "json.hpp"
#include <fstream>
#include <iostream>

using json = nlohmann::json;
#include <string>


int main(int argc, char **argv)
{
  std::ifstream ff("data.json");
  json data = json::parse(ff);

  //INITIALIZING VARIABLES FROM DATA.JSON
  std::string funString = data.value("fun","");
  std::string dx0_funString = data.value("dx[0]fun","");
  std::string dx1_funString = data.value("dx[1]fun","");
  double alpha = data.value("alpha", 0.2);
  const unsigned int max_it =  data.value("max_it", 50);
  const double tol_fun = data.value("tol_fun", 1e-6);
  const double tol_x = data.value("tol_x", 1e-6);
  std::array<double,2> x0{{0,0}}; 

  //PASSING FROM STRING TO FUNCTION WITH MUPARSERX
  using namespace MuParserInterface;
  muParserXInterface<2> f; 
  f.set_expression(funString);
  muParserXInterface<2> dfx0; 
  dfx0.set_expression(dx0_funString);
  muParserXInterface<2> dfx1; 
  dfx1.set_expression(dx1_funString);
  
  // GRADIENT METHOD (Exponential decay)
  std::cout<< "Calculating the minimum of the function with the gradient method (with exponential decay)" << std::endl;
  GradientMethod minimizer(f, dfx0, dfx1, alpha, max_it, tol_x, tol_fun, alpha_strategies::Exponencial_decay);
  minimizer.minimize(x0);   
  std::array<double,2> out = minimizer.get_result();
  std::cout<< "The minimum is found in the point: " << out[0] << ", " << out[1] << '.' << std::endl;
  std::cout<< "The value of the minimum is: " << f(out) << '.' << std::endl;
  std::cout<< "Calculated it in " << minimizer.get_iter() << " iterations." << std::endl << std::endl;

  // GRADIENT METHOD (Inverse decay)
  std::cout<< "Calculating the minimum of the function with the gradient method (with inverse decay)" << std::endl;
  GradientMethod minimizer2(f, dfx0, dfx1, alpha, max_it, tol_x, tol_fun, alpha_strategies::Inverse_decay);
  minimizer2.minimize(x0);   
  std::array<double,2> out2 = minimizer2.get_result();
  std::cout<< "The minimum is found in the point: " << out2[0] << ", " << out2[1] << '.' << std::endl;
  std::cout<< "The value of the minimum is: " << f(out2) << '.' << std::endl;
  std::cout<< "Calculated it in " << minimizer2.get_iter() << " iterations." << std::endl << std::endl;

  // GRADIENT METHOD (Aproximate line serach)
  std::cout<< "Calculating the minimum of the function with the gradient method (with aproximate_line_search)" << std::endl;
  GradientMethod minimizer3(f, dfx0, dfx1, alpha, max_it, tol_x, tol_fun, alpha_strategies::Aproximate_line_search);
  minimizer3.minimize(x0);   
  std::array<double,2> out3 = minimizer3.get_result();
  std::cout<< "The minimum is found in the point: " << out3[0] << ", " << out3[1] << '.' << std::endl;
  std::cout<< "The value of the minimum is: " << f(out3) << '.' << std::endl;
  std::cout<< "Calculated it in " << minimizer3.get_iter() << " iterations." << std::endl << std::endl;

  // GRADIENT METHOD (Exponential decay with momentum)
  std::cout<< "Calculating the minimum of the function with the gradient method (with exponential decay)" << std::endl;
  GradientMethod minimizer4(f, dfx0, dfx1, alpha, max_it, tol_x, tol_fun, alpha_strategies::Exponencial_decay);
  minimizer4.minimize_with_momentum(x0);   
  std::array<double,2> out4 = minimizer4.get_result();
  std::cout<< "The minimum is found in the point: " << out4[0] << ", " << out4[1] << '.' << std::endl;
  std::cout<< "The value of the minimum is: " << f(out4) << '.' << std::endl;
  std::cout<< "Calculated it in " << minimizer4.get_iter() << " iterations." << std::endl << std::endl;

  // GRADIENT METHOD (Inverse decay with momentum)
  std::cout<< "Calculating the minimum of the function with the gradient method (with inverse decay)" << std::endl;
  GradientMethod minimizer5(f, dfx0, dfx1, alpha, max_it, tol_x, tol_fun, alpha_strategies::Inverse_decay);
  minimizer5.minimize_with_momentum(x0);   
  std::array<double,2> out5 = minimizer2.get_result();
  std::cout<< "The minimum is found in the point: " << out5[0] << ", " << out5[1] << '.' << std::endl;
  std::cout<< "The value of the minimum is: " << f(out5) << '.' << std::endl;
  std::cout<< "Calculated it in " << minimizer5.get_iter() << " iterations." << std::endl << std::endl;

}
