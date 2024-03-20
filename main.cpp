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
  const unsigned int max_it =  data.value("max_it",100);
  const double tol_fun = data.value("tol_res", 1e-7);
  const double tol_x = data.value("tol_x", 1e-7);

  //PRINTING DATA TO PROVE IT WORKS
  std::cout<< "PRINTING THE DATA TAKEN WITH JSON: " << max_it << ' ' << tol_fun << ' ' << tol_x << std::endl;
  std::cout<< funString << std::endl << dx0_funString << std::endl<< std::endl;

  //PASSING FROM STRING TO FUNCTION WITH MUPARSERX
  using namespace MuParserInterface;
  muParserXInterface<2> f; 
  f.set_expression(funString);
  muParserXInterface<2> dfx0; 
  dfx0.set_expression(dx0_funString);
  muParserXInterface<2> dfx1; 
  dfx1.set_expression(dx1_funString);

  //TESTING IT WORKS
  std::cout << "TESTING THE MUPARSERX: " << dfx1({1,1}) << std::endl<< std::endl;

  //JOINING THE TWO FUNCTIONS IN ONE ONLY
  std::vector< std::function<double(std::array<double,2>)>> df(2);
  df[0] = dfx0; df[1] = dfx1;
  std::cout<< "TESTING THE JOIN OF THE FUNCTIONS: " << df[0]({1,1}) << std::endl << std::endl;

  
  // initialize minimizer
  
  GradientMethod minimizer(f, dfx0, dfx1, alpha, max_it, tol_x, tol_fun);
  std::array<double,2> x0{{0,0}};
  minimizer.minimize(x0);
  std::array<double,2> minimizer.get_result();

  // output the results
  //std::cout << "x =    " << solver.get_result() << std::endl;
  //std::cout << "r =    " << solver.get_residual() << std::endl;
  //std::cout << "iter = " << solver.get_iter() << std::endl;

}
