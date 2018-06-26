//Basic syntax test

#include <iostream>
#include "exprtk.hpp"

int main(){
  double x = 1;
  double sol;

  if (int(x==1 and 2==2)) std::cout << "true\n";
  sol = exp(x - 1);
  std::cout << sol << std::endl;
  return 0;
}
