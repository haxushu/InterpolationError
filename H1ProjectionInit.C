#include "H1ProjectionInit.h"
#include "../operator/play.h"


#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/mesh_generation.h"

// This example will solve a linear transient system,
// so we need to include the TransientLinearImplicitSystem definition.
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/vector_value.h"

// The definition of a geometric element
#include "libmesh/elem.h"

// To impose Dirichlet boundary conditions
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/analytic_function.h"
#include "libmesh/zero_function.h"
#include "libmesh/boundary_info.h"

  /**
   * Evaluate the Jacobian of the nonlinear system.
   */
void H1ProjectionInit::jacobian (const NumericVector<Number> & soln,
                        SparseMatrix<Number> & jacobian,
                        NonlinearImplicitSystem & /*sys*/)
{
        //  std::cout << "before jacobian" << std::endl;
  
  const MeshBase & mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  NonlinearImplicitSystem & system =
    es.get_system<NonlinearImplicitSystem>(sys_name);

  const unsigned int u_var = system.variable_number (var_name);

  const DofMap & dof_map = system.get_dof_map();

  FEType fe_type = dof_map.variable_type(u_var);
  std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);

  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

  DenseMatrix<Number> Ke;

  std::vector<dof_id_type> dof_indices;

  jacobian.zero();

  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      dof_map.dof_indices (elem, dof_indices);
      const unsigned int n_dofs = dof_indices.size();

      fe->reinit (elem);

      Ke.resize (n_dofs, n_dofs);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

          for (unsigned int dof_i=0; dof_i<n_dofs; dof_i++){

              for (unsigned int dof_j=0; dof_j<n_dofs; dof_j++){
                  
                  {
                      Ke(dof_i,dof_j) += JxW[qp] * dphi[dof_i][qp] * dphi[dof_j][qp];
                  }
              }
                         
          }
        }
      
      const BoundaryInfo boundary_info = mesh.get_boundary_info();
      // boundary_info.print_info();
      // auto s =  boundary_info.get_boundary_ids();
      // for (auto v : s)  out << v << " "; out << std::endl;
      // out << "mesh n_boundary_ids: " << boundary_info.n_boundary_ids() << std::endl;

      for (const auto &node : elem->node_ref_range()){

        // if (dim == 2){
        //   if (boundary_info.has_boundary_id(&node,3) || boundary_info.has_boundary_id(&node,1)) {
        //     dof_id_type local_dof_id = elem->local_node(node.id()); 
        //     for (unsigned int dof_j=0; dof_j<n_dofs; dof_j++) Ke(local_dof_id, dof_j) =0.0;
        //     Ke(local_dof_id,local_dof_id) = 1;
        //   }
        // }
        if (dim == 3 || dim == 1 || dim == 2){ 
          if (boundary_info.n_boundary_ids(&node)) {
            dof_id_type local_dof_id = elem->local_node(node.id()); 
            for (unsigned int dof_j=0; dof_j<n_dofs; dof_j++) Ke(local_dof_id, dof_j) =0.0;
            Ke(local_dof_id,local_dof_id) = 1;
          }
        }
      }
      
      jacobian.add_matrix (Ke, dof_indices);
    }

    // std::cout << "after jacobian" << std::endl;
    // jacobian.print_matlab();
}

/**
 * Evaluate the residual of the nonlinear system.
 */
void H1ProjectionInit::residual (const NumericVector<Number> & soln,
                        NumericVector<Number> & residual,
                        NonlinearImplicitSystem & /*sys*/)
{
              // std::cout << "before residual" << std::endl;

  const MeshBase & mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  NonlinearImplicitSystem & system =
    es.get_system<TransientNonlinearImplicitSystem>(sys_name);

  const unsigned int u_var = system.variable_number (var_name);

  const DofMap & dof_map = system.get_dof_map();

  FEType fe_type = dof_map.variable_type(u_var);
  std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, TENTH);
  fe->attach_quadrature_rule (&qrule);

  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

  const std::vector<Point> & q_point = fe->get_xyz();

  DenseVector<Number> Re;

  std::vector<dof_id_type> dof_indices;

  residual.zero();

  // std::cout << "before enumeration" << std::endl;

  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      dof_map.dof_indices (elem, dof_indices);
      const unsigned int n_dofs = dof_indices.size();

      fe->reinit (elem);

      Re.resize (n_dofs);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

          Gradient du_value = 0;
          
          for (unsigned int dof_j=0; dof_j<n_dofs; dof_j++){     
            du_value += dphi[dof_j][qp] * soln (dof_indices[dof_j]);
          }
            
          Gradient duxy = initial_dvalue(q_point[qp], es.parameters);

          for (unsigned int dof_i=0; dof_i<n_dofs; dof_i++){
            
            Re(dof_i) += JxW[qp] * dphi[dof_i][qp] * (du_value-duxy);       

          }

        }

      const BoundaryInfo boundary_info = mesh.get_boundary_info();

      for (const auto &node : elem->node_ref_range()){
        dof_id_type dof_id = node.dof_number(0, 0, 0);     // Can you understand the mapping rule ???
        // printf("rank %d node id %u dof_id %u\n", node.processor_id(), node.id(), dof_id);
        dof_id_type local_dof_id = elem->local_node(node.id()); 
        // if (dim == 2)
        //   if (boundary_info.has_boundary_id(&node,3) || boundary_info.has_boundary_id(&node,1)) {
        //     Re(local_dof_id) = soln(dof_id) - initial_value(node, es.parameters);
        //   }
        if (dim == 3 || dim == 1 || dim == 2)  
          if (boundary_info.n_boundary_ids(&node)) {
            Re(local_dof_id) = soln(dof_id) - initial_value(node, es.parameters);
          }
      }

      // dof_map.constrain_element_vector (Re, dof_indices);
      residual.add_vector (Re, dof_indices);
    }
    
    // std::cout << "after residual" << std::endl;
    // residual.print_global();

}


