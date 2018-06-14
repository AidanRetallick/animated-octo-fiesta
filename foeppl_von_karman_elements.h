// Header file for the Biharmonic Bell elements
#ifndef OOMPH_MYBIHARMONIC_ELEMENTS_HEADER
#define OOMPH_MYBIHARMONIC_ELEMENTS_HEADER


#include<sstream>

//OOMPH-LIB headers
#include "generic.h"
#include "MyBellShape.h"


namespace oomph
{
//=============================================================
/// A class for all subparametric elements that solve the 2D-
/// Biharmonic equations.
/// \f[
/// \frac{\partial^4 u}{\partial x_i^4} = f(x_j)
/// \f]
/// This contains the generic maths. Shape functions, geometric
/// mapping etc. must get implemented in derived class.
//=============================================================
template <unsigned DIM, unsigned NNODE_1D>
class FoepplVonKarmanEquations : public virtual TElement<DIM,NNODE_1D>
{

public:

 /// \short Function pointer to pressure function fct(x,f(x)) --
 /// x is a Vector!
 typedef void (*PressureFctPt)(const Vector<double>& x, double& f);


 /// \short Function pointer to in plane forcing function  fct(x,g(x)) --
 /// x is a Vector!
 typedef void (*InPlaneForcingFctPt)(const Vector<double>& x,
                                            Vector<double>& forcing);

 /// \short Function pointer to the Error Metric we are using
 ///  e.g could be that we are just interested in error on w etc.
 typedef void (*ErrorMetricFctPt)(const Vector<double>& x, const 
  Vector<double>& u, const Vector<double>& u_exact, double& error, double& norm);

 /// \shot Function pointer to the Error Metric we are using if we want multiple
 ///  errors.  e.g could be we want errors seperately on each displacment
 typedef void (*MultipleErrorMetricFctPt)(const Vector<double>& x, const 
  Vector<double>& u, const Vector<double>& u_exact, Vector<double>& error, 
  Vector<double>& norm);

 virtual void fix_out_of_plane_displacement_dof(const unsigned& dof_number, 
const unsigned& boundary_number, const PressureFctPt& w)=0;

 virtual void fix_in_plane_displacement_dof(const unsigned& dof_number,
const unsigned& boundary_number, const PressureFctPt& u)=0;

 virtual void precompute_association_matrix(DenseMatrix<double>& m)=0;
 virtual double n_basis_functions()=0;
 virtual double n_basic_basis_functions()=0;

 
 public:

 // Get pointer to association matrix 
 DenseMatrix<double> *get_association_matrix_pt()const {return Association_matrix_pt;};

 public: 

 /// Eta
 const double &eta() const {return *Eta_pt;}

 /// Pointer to eta
 double* &eta_pt() {return Eta_pt;}

 virtual void pin_all_deflection_dofs() const=0;

 /// Constructor (must initialise the Pressure_fct_pt to null)
 FoepplVonKarmanEquations() : Pressure_fct_pt(0),
   In_plane_forcing_fct_pt(0),Number_of_internal_dofs(0),
   Number_of_internal_dof_types(0), Association_matrix_pt(0),
Error_metric_fct_pt(0), Multiple_error_metric_fct_pt(0) {}

 /// Broken copy constructor
 FoepplVonKarmanEquations(const FoepplVonKarmanEquations& dummy)
  {
   BrokenCopy::broken_copy("FoepplVonKarmanEquations");
  }

 /// Broken assignment operator
 void operator=(const FoepplVonKarmanEquations&)
  {
   BrokenCopy::broken_assign("FoepplVonKarmanEquations");
  }


 
 /// \short Return the index at which the unknown value
 /// is stored.
 /// In derived multi-physics elements, this function should be overloaded
 /// to reflect the chosen storage scheme. Note that these equations require
 /// that the unknown is always stored at the same index at each node.
 virtual inline unsigned u_index_biharmonic() const {return this->required_nvalue(0);}

 /// Output with default number of plot points
 void output(std::ostream &outfile)
  {
   const unsigned n_plot=5;
   FoepplVonKarmanEquations::output(outfile,n_plot);
  }

 /// \short Output FE representation of soln: x,y,u or x,y,z,u at
 /// n_plot^DIM plot points
 void output(std::ostream &outfile, const unsigned &n_plot);

 /// C_style output with default number of plot points
 void output(FILE* file_pt)
  {
   const unsigned n_plot=5;
   FoepplVonKarmanEquations::output(file_pt,n_plot);
  }

 /// \short C-style output FE representation of soln: x,y,u or x,y,z,u at
 /// n_plot^DIM plot points
 void output(FILE* file_pt, const unsigned &n_plot);

 /// Output exact soln: x,y,u_exact or x,y,z,u_exact at n_plot^DIM plot points
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 FiniteElement::SteadyExactSolutionFctPt exact_soln_pt);

 /// \short Output exact soln: x,y,u_exact or x,y,z,u_exact at
 /// n_plot^DIM plot points (dummy time-dependent version to
 /// keep intel compiler happy)
 virtual void output_fct(std::ostream &outfile, const unsigned &n_plot,
                         const double& time,
                         FiniteElement::UnsteadyExactSolutionFctPt
                         exact_soln_pt)
  {
   throw OomphLibError(
    "There is no time-dependent output_fct() for these elements ",
    OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }


 /// Get error against and norm of exact solution
 void compute_error(std::ostream &outfile,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                    double& error, double& norm);

 /// Get error against and norm of exact solution
 void compute_error_in_deflection(std::ostream &outfile,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                    double& error, double& norm);

 /// Dummy, time dependent error checker
 void compute_error(std::ostream &outfile,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                    const double& time, double& error, double& norm)
  {
   throw OomphLibError(
    "There is no time-dependent compute_error() for these elements",
    OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

 // Fill in the stress tensor
 void get_sigma(DenseMatrix<double>& sigma, const DenseMatrix<double>& grad_u,
  const DenseMatrix<double>& grad_w )const
  {
   // Poisson ratio
   double nu(get_nu());
   // Truncated Green Lagrange strain tensor
   DenseMatrix<double> epsilon(2,2,0.0);
   for(unsigned alpha=0;alpha<2;++alpha)
    {
     for(unsigned beta=0;beta<2;++beta)
      {
       // Truncated Green Lagrange strain tensor
       epsilon(alpha,beta) += 0.5* grad_u(alpha,beta) + 0.5*grad_u(beta,alpha)
                           +(*Eta_pt)* 0.5*grad_w(0,alpha)*grad_w(0,beta);
      }
    }
   
   // Now construct the Stress
   for(unsigned alpha=0;alpha<2;++alpha)
    {
     for(unsigned beta=0;beta<2;++beta)
      {
       // The Laplacian term: Trace[ \epsilon ] I
       // \nu * \epsilon_{\alpha \beta} delta_{\gamma \gamma}
       sigma(alpha,alpha) += nu*epsilon(beta,beta)/(1-nu*nu);
      
       // The scalar transform term: \epsilon
       // (1-\nu) * \epsilon_{\alpha \beta}
       sigma(alpha,beta) += (1-nu)* epsilon(alpha,beta)/(1-nu*nu);
      }
    }
 
  }
 /// Access function: Pointer to pressure function
 PressureFctPt& pressure_fct_pt() {return Pressure_fct_pt;}

 /// Access function: Pointer to pressure function. Const version
 PressureFctPt pressure_fct_pt() const {return Pressure_fct_pt;}

 /// Access function: Pointer to error metric function
 ErrorMetricFctPt& error_metric_fct_pt() {return Error_metric_fct_pt;}

 /// Access function: Pointer to multiple error metric function
 MultipleErrorMetricFctPt& multiple_error_metric_fct_pt() 
  {return Multiple_error_metric_fct_pt;}

 /// Access function: Pointer to error metric function function
 ErrorMetricFctPt error_metric_fct_pt() const {return Error_metric_fct_pt;}

 /// Access function: Pointer to multiple error metric function
 MultipleErrorMetricFctPt multiple_error_metric_fct_pt() const 
  {return Multiple_error_metric_fct_pt;}

 /// Access function: Pointer to in plane forcing function
 InPlaneForcingFctPt& in_plane_forcing_fct_pt()
  {return In_plane_forcing_fct_pt;}

 /// Access function: Pointer to in plane forcing function. Const version
 InPlaneForcingFctPt in_plane_forcing_fct_pt() const
  {return In_plane_forcing_fct_pt;}

 ///Access function to the Poisson ratio.
 const double*& nu_pt() {return Nu_pt;}

 ///Access function to the Poisson ratio (const version)
 const double& get_nu() const {return *Nu_pt;}

 // Get the kth dof type at internal point l
 virtual double get_w_bubble_dof(const unsigned& l, const unsigned& k) const =0;

 // Get the kth equation at internal point l
 virtual int local_w_bubble_equation(const unsigned& l, const unsigned& k)
   const =0;

 /// Get pressure term at (Eulerian) position x. This function is
 /// virtual to allow overloading in multi-physics problems where
 /// the strength of the pressure function might be determined by
 /// another system of equations.
 inline virtual void get_pressure_biharmonic(const unsigned& ipt,
                                        const Vector<double>& x,
                                        double& pressure) const
  {
   //If no pressure function has been set, return zero
   if(Pressure_fct_pt==0)
    {
     pressure = 0.0;
    }
   else
    {
     // Get pressure strength
     (*Pressure_fct_pt)(x,pressure);
    }
  }

 /// Add the element's contribution to its residual vector (wrapper)
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   // Precompute the association matrix
   DenseMatrix<double> conversion_matrix (n_basis_functions(),
    n_basic_basis_functions(),0.0);
   this->precompute_association_matrix(conversion_matrix);
   this->Association_matrix_pt=&conversion_matrix;

   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   fill_in_generic_residual_contribution_biharmonic(
    residuals,GeneralisedElement::Dummy_matrix,0);

   Association_matrix_pt=0;
  }


 /// Add the element's contribution to its residual vector and
 /// element Jacobian matrix (wrapper)
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
   DenseMatrix<double> conversion_matrix (n_basis_functions(),
    n_basic_basis_functions(),0.0);
   this->precompute_association_matrix(conversion_matrix);
   this->Association_matrix_pt=&conversion_matrix;

   //Call the generic routine with the flag set to 1
   fill_in_generic_residual_contribution_biharmonic(residuals,jacobian,1);

   Association_matrix_pt=0;
  }

 /// Add the element's contribution to its residual vector and
 /// element Jacobian matrix (wrapper)
 void fill_in_contribution_to_jacobian_and_mass_matrix(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian,DenseMatrix<double> &mass_matrix)
  {
   //Call fill in Jacobian 
   fill_in_contribution_to_jacobian(residuals,jacobian);
   // There is no mass matrix: we will just want J w = 0

   // -- COPIED FROM DISPLACMENT FVK EQUATIONS --
   // Dummy diagonal (won't result in global unit matrix but
   // doesn't matter for zero eigenvalue/eigenvector
   unsigned ndof=mass_matrix.nrow();
   for (unsigned i=0;i<ndof;i++)
    {
     mass_matrix(i,i)+=1.0;
    }
  }


 /// \short Return FE representation of unknown values u(s)
 /// at local coordinate s
 inline Vector<double> interpolated_u_biharmonic(const Vector<double> &s, bool
output_stress_flag=false) const
  {
   //Find number of position dofs
   const unsigned n_position_type = this->nnodal_position_type();
   // Find the internal dofs
   const unsigned n_b_position_type = this->Number_of_internal_dof_types;
   //Find out how many nodes there are
   const unsigned n_w_node = 3;
  //Find out how many nodes there are
   const unsigned n_node = this->nnode();
   //Find out how many internal points there are
   const unsigned n_b_node = this->Number_of_internal_dofs;
   //Get the index at which the unknown is stored
   // const unsigned u_nodal_index = u_index_biharmonic();

   //Local c1-shape funtion
   Shape psi(n_w_node,n_position_type),test(n_w_node,n_position_type),
    psi_b(n_b_node,n_b_position_type),test_b(n_b_node,n_b_position_type);
   

   DShape dpsi_dxi(n_w_node,n_position_type,DIM),dtest_dxi(n_w_node,n_position_type,DIM),
    dpsi_b_dxi(n_b_node,n_b_position_type,DIM),dtest_b_dxi(n_b_node,n_b_position_type,DIM),
    d2psi_dxi2(n_w_node,n_position_type,3), d2test_dxi2(n_w_node,n_position_type,3),
    d2psi_b_dxi2(n_b_node,n_b_position_type,3), d2test_b_dxi2(n_b_node,n_b_position_type,3);
   
   // In--plane dofs
   Shape psi_u(n_node);
   DShape dpsi_u(n_node,DIM);
   Shape test_u(n_node);
   DShape dtest_u(n_node,DIM);

   //Initialise value of u
   Vector<double> interpolated_u(15,0.0);
   //Find values of c1-shape function
   d2shape_and_d2test_eulerian_biharmonic(s,psi,psi_b,dpsi_dxi,dpsi_b_dxi,
    d2psi_dxi2,d2psi_b_dxi2,test,test_b,dtest_dxi,dtest_b_dxi,d2test_dxi2,
    d2test_b_dxi2);

   // Get shape and test
   dshape_u_and_dtest_u_eulerian_biharmonic(s,psi_u,dpsi_u,test_u,dtest_u);
   //Interpolated unknown
   for(unsigned l=0;l<n_w_node;l++)
   {
    for(unsigned k=0;k<n_position_type;k++)
     {
      interpolated_u[0] += this->nodal_value(l,k+2)*psi(l,k);
      interpolated_u[1] += this->nodal_value(l,k+2)*dpsi_dxi(l,k,0);
      interpolated_u[2] += this->nodal_value(l,k+2)*dpsi_dxi(l,k,1);
      interpolated_u[3] += this->nodal_value(l,k+2)*d2psi_dxi2(l,k,0);
      interpolated_u[4] += this->nodal_value(l,k+2)*d2psi_dxi2(l,k,1);
      interpolated_u[5] += this->nodal_value(l,k+2)*d2psi_dxi2(l,k,2);
     }
   }

   // Bubble dofs
   for(unsigned l=0;l<Number_of_internal_dofs;l++)
   {
    for(unsigned k=0;k<Number_of_internal_dof_types;k++)
     {
      double u_value = get_w_bubble_dof(l,k);
      interpolated_u[0] += u_value * psi_b(l,k);
      interpolated_u[1] += u_value*dpsi_b_dxi(l,k,0);
      interpolated_u[2] += u_value*dpsi_b_dxi(l,k,1);
      interpolated_u[3] += u_value*d2psi_b_dxi2(l,k,0);
      interpolated_u[4] += u_value*d2psi_b_dxi2(l,k,1);
      interpolated_u[5] += u_value*d2psi_b_dxi2(l,k,2);
     }
   }
   // Now for the displacement 
   for(unsigned l=0; l< n_node;++l)
    {
     // Now for the two in--plane displacements
     interpolated_u[6] += this->nodal_value(l,0)*psi_u(l);
     interpolated_u[7] += this->nodal_value(l,1)*psi_u(l);
     // IF flag
     if(output_stress_flag)
      {
       for(unsigned i=0; i<2; ++i)
        {
        interpolated_u[8 +i] += this->nodal_value(l,0)*dpsi_u(l,i);
        interpolated_u[10+i] += this->nodal_value(l,1)*dpsi_u(l,i);
        }
      }
    }

   // IF flag
   if(output_stress_flag)
    {
     // For stress
     DenseMatrix<double> interpolated_dwdxi(1,DIM,0.0);
     DenseMatrix<double> interpolated_dudxi(2,DIM,0.0);
     interpolated_dwdxi(0,0)=interpolated_u[1];
     interpolated_dwdxi(0,1)=interpolated_u[2];
     interpolated_dudxi(0,0)=interpolated_u[8];
     interpolated_dudxi(0,1)=interpolated_u[9];
     interpolated_dudxi(1,0)=interpolated_u[10];
     interpolated_dudxi(1,1)=interpolated_u[11];
     // Get Stress
     DenseMatrix<double> sigma(2,2,0.0);
     get_sigma(sigma,interpolated_dudxi, interpolated_dwdxi);
     interpolated_u[12] = sigma(0,0);
     interpolated_u[13] = sigma(0,1);
     interpolated_u[14] = sigma(1,1);
    }

   return(interpolated_u);
  }

 /// \short Self-test: Return 0 for OK
 unsigned self_test();

 /// \short get the coordinate
 virtual void get_coordinate_x(const Vector<double>& s, Vector<double>& x) const=0;
 
// /// Wrapper
// void interpolated_x(const Vector<double>& s, Vector<double>& x) const
//  {
//   this->get_coordinate_x(s,x);
//  }

protected:
 virtual void shape_u(const Vector<double> &s,  Shape &psi) const=0;
 /// \short Shape/test functions and derivs w.r.t. to global coords at
 /// local coord. s; return  Jacobian of mapping
 virtual double d2shape_and_d2test_eulerian_biharmonic(const Vector<double> &s,
  Shape &psi, Shape& psi_b, DShape &dpsi_dx, DShape &dpsi_b_dx,
  DShape &d2psi_dx2,DShape& d2psi_b_dx2,
  Shape &test, Shape& test_b, DShape &dtest_dx, DShape &dtest_b_dx,
  DShape &d2test_dx2,DShape& d2test_b_dx2) const=0;

 /// \short Shape/test functions and derivs w.r.t. to global coords at
 /// local coord. s; return  Jacobian of mapping
 virtual double dshape_and_dtest_eulerian_biharmonic(const Vector<double> &s,
  Shape &psi, Shape& psi_b, DShape &dpsi_dx, DShape &dpsi_b_dx,
  Shape &test, Shape& test_b, DShape &dtest_dx, DShape &dtest_b_dx) const=0;

 /// \short Shape/test functions at local coordinate s
 virtual void shape_and_test_biharmonic(const Vector<double> &s,
  Shape &psi, Shape& psi_b, Shape &test, Shape& test_b) const=0;
 
 /// \short in--plane Shape/test functions at local coordinate s
 virtual double dshape_u_and_dtest_u_eulerian_biharmonic(const Vector<double> &s,
  Shape &psi,DShape &dpsidx, Shape &test,DShape &dtestdx) const=0;

// /// \short Shape/test functions and derivs w.r.t. to global coords at
// /// local coord. s; return  Jacobian of mapping
// virtual double d2shape_and_d2test_eulerian_at_knot_biharmonic(const
//  unsigned ipt, Shape &psi, Shape& psi_b, DShape &dpsi_dx, DShape &dpsi_b_dx,
//  DShape &d2psi_dx2,DShape& d2psi_b_dx2,
//  Shape &test, Shape& test_b, DShape &dtest_dx, DShape &dtest_b_dx,
//  DShape &d2test_dx2,DShape& d2test_b_dx2) const=0;
//
// /// \short Shape/test functions and derivs w.r.t. to global coords at
// /// local coord. s; return  Jacobian of mapping
// virtual double dshape_and_dtest_eulerian__at_knot_biharmonic(const
//  unsigned& ipt, Shape &psi, Shape& psi_b, DShape &dpsi_dx, DShape &dpsi_b_dx,
//  Shape &test, Shape& test_b, DShape &dtest_dx, DShape &dtest_b_dx) const=0;
//
// /// \short Shape/test functions at integral point ipt
// virtual void shape_and_test_biharmonic_at_knot(const unsigned& ipt,
//  Shape &psi, Shape& psi_b, Shape &test, Shape& test_b) const=0;

 /// \short Compute element residual Vector only (if flag=and/or element
 /// Jacobian matrix
 virtual void fill_in_generic_residual_contribution_biharmonic(
  Vector<double> &residuals, DenseMatrix<double> &jacobian,
  const unsigned& flag);

 /// Pointer to pressure function:
 PressureFctPt Pressure_fct_pt;

 /// Pointer to in plane forcing function (i.e. the shear force applied to
 /// the face)
 InPlaneForcingFctPt In_plane_forcing_fct_pt;

 /// Pointer to Poisson ratio, which this element cannot modify
 const double* Nu_pt;

 /// Pointer to global eta
 double *Eta_pt;

 /// \short unsigned that holds the internal 'bubble' dofs the element has -
 // zero for Bell Elements and 3 for C1 curved elements
 unsigned Number_of_internal_dofs;

 /// \short unsigned that holds the number of types of degree of freedom at each
 // internal point that the element has zero for Bell Elements and 1
 // for C1 curved elements
 unsigned Number_of_internal_dof_types;

 /// Default value for physical constants
 static double Default_Physical_Constant_Value;

 /// Pointer to error metric
 ErrorMetricFctPt Error_metric_fct_pt;

 /// Pointer to error metric when we want multiple errors
 MultipleErrorMetricFctPt Multiple_error_metric_fct_pt;


 protected:

 /// Pointer to precomputed matrix that associates shape functions to monomials
 DenseMatrix<double> *Association_matrix_pt;
};


} //end namespace oomph
#endif

