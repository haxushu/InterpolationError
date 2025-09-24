// C++ Includes
#include <cmath>

// Mesh library includes
#include "libmesh/libmesh_common.h"
#include "libmesh/equation_systems.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This function will initialize the system.
// Initialization functions are optional for systems.  They allow
// you to specify the initial values of the solution.  If an
// initialization function is not provided then the default (0)
// solution is provided.

void init_pj0_u (EquationSystems & es,
              const std::string & system_name);


Number initial_value_u (const Point & p,
                    const Parameters & parameters,
                    const std::string &,
                    const std::string &);

Number initial_value_u (const Point & p,
                    const Parameters & parameters);


Gradient initial_dvalue_u (const Point & p,
                    const Parameters & parameters,
                    const std::string &,
                    const std::string &);

Gradient initial_dvalue_u (const Point & p,
                    const Parameters & parameters);

// Number boundary_value_u (const Point & p, const Real time);
