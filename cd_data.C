// C++ Includes
#include <cmath>

// Mesh library includes
#include "libmesh/libmesh_common.h"
#include "libmesh/equation_systems.h"
#include "../operator/play.h"
#include "cd_data.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

extern int CASE;

Number initial_value_u (const Point & p,
                    const Parameters & parameters,
                    const std::string &,
                    const std::string &)
{
  long double pi = acos((long double)-1.0);
  switch (CASE)
  {
    case 1: return p(0)*(1-p(0))*p(1)*(1-p(1));
    case 2: return 2*sin(7*p(0));
    case 3: return sin(pi*p(0))*sin(pi*p(1))*sin(pi*p(2));
    default: return exp(p(0));
  }
}

Gradient initial_dvalue_u (const Point & p,
                    const Parameters & parameters,
                    const std::string &,
                    const std::string &)
{
  long double pi = acos((long double)-1.0);
  switch (CASE)
  {
    case 1: return Gradient((1-2*p(0))*p(1)*(1-p(1)), p(0)*(1-p(0))*(1-2*p(1)));
    case 2: return Gradient(14*cos(7*p(0)));
    case 3: return Gradient(pi*cos(pi*p(0))*sin(pi*p(1))*sin(pi*p(2)), pi*cos(pi*p(1))*sin(pi*p(0))*sin(pi*p(2)), pi*cos(pi*p(2))*sin(pi*p(1))*sin(pi*p(0)));
    default: return Gradient(exp(p(0)));
  }
}

Number initial_value_u (const Point & p,
                    const Parameters & parameters)
{
  long double pi = acos((long double)-1.0);
  switch (CASE)
  {
    case 1: return p(0)*(1-p(0))*p(1)*(1-p(1));;
    case 2: return 2*sin(7*p(0));
    case 3: return sin(pi*p(0))*sin(pi*p(1))*sin(pi*p(2));
    default: return exp(p(0));
  }
}

Gradient initial_dvalue_u (const Point & p,
                    const Parameters & parameters)
{
  long double pi = acos((long double)-1.0);
  switch (CASE)
  {
    case 1: return Gradient((1-2*p(0))*p(1)*(1-p(1)), p(0)*(1-p(0))*(1-2*p(1)));
    case 2: return Gradient(14*cos(7*p(0)));
    case 3: return Gradient(pi*cos(pi*p(0))*sin(pi*p(1))*sin(pi*p(2)), pi*cos(pi*p(1))*sin(pi*p(0))*sin(pi*p(2)), pi*cos(pi*p(2))*sin(pi*p(1))*sin(pi*p(0)));
    default: return Gradient(exp(p(0)));
  }
}


// Number boundary_value_u (const Point & p, const Real time)
// {
//     switch (CASE)
//     {
//        case 1: return (1-p(0))*(-time) + p(0)*(1-time);
//        default: return (1-p(0))*2*sin(7*time) + p(0)*2*cos(7*time)*sin((Real)7.0);
//     }
// }

