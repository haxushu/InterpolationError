// The nonlinear solver and system we will be using
#include "libmesh/nonlinear_solver.h"
#include "libmesh/nonlinear_implicit_system.h"
// Bring in everything from the libMesh namespace
using namespace libMesh;

class L2ProjectionInit : public NonlinearImplicitSystem::ComputeResidual,
                                   public NonlinearImplicitSystem::ComputeJacobian
{
private:
  EquationSystems & es;
  std::string sys_name, var_name; 
  Number (* initial_value) (const Point & p,
                    const Parameters & parameters);
public:

  L2ProjectionInit (EquationSystems & es_in, std::string sys_name_in, std::string var_name_in, Number initial_value_in (const Point & p,
                    const Parameters & parameters)) :
    es(es_in), sys_name(sys_name_in), var_name(var_name_in), initial_value(initial_value_in)
  {}

  /**
   * Evaluate the Jacobian of the nonlinear system.
   */
  virtual void jacobian (const NumericVector<Number> & soln,
                         SparseMatrix<Number> & jacobian,
                         NonlinearImplicitSystem & /*sys*/);
  

  /**
   * Evaluate the residual of the nonlinear system.
   */
  virtual void residual (const NumericVector<Number> & soln,
                         NumericVector<Number> & residual,
                         NonlinearImplicitSystem & /*sys*/);
  


};
