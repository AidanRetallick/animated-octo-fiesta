//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================

#ifndef OOMPH_BIHARMONIC_CURVED_BELL_ELEMENTS_HEADER
#define OOMPH_BIHARMONIC_CURVED_BELL_ELEMENTS_HEADER

#include "foeppl_von_karman_elements.h"
#include "MyBellShape.h"
#include "C1_curved_elements.h"
#include "my_geom_object.h"

namespace oomph
{
//===============================================================================
/// FoepplVonKarmanC1CurvedBellElement elements are a subparametric scheme
/// with  linear Lagrange interpolation for approximating the geometry and
/// the C1-functions for approximating variables.
//==============================================================================

template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
class FoepplVonKarmanC1CurvedBellElement : public virtual
 FoepplVonKarmanEquations<DIM,NNODE_1D>
{
public:
 /// \short Function pointer to pressure function fct(x,f(x)) --
 /// x is a Vector!
 typedef void (*DisplacementFctPt)(const Vector<double>& x, double& f);

 /// \short function to pin all deflection dofs
 void pin_all_deflection_dofs() const;

 /// \short function to pin particular out--of--plane displacement dofs
 void fix_out_of_plane_displacement_dof(const unsigned& dof_number, const unsigned& b,
const DisplacementFctPt& w);

 /// \short function to pin particular in--plane displacement dofs
 void fix_in_plane_displacement_dof(const unsigned& dof_number, const unsigned&
b, const DisplacementFctPt& u);

 /// \short Return boundary order
 unsigned boundary_order() {return BOUNDARY_ORDER;}
 
 /// \short get the location of the internal dofs
 inline void get_internal_dofs_location(const unsigned& s, Vector<double>& x) const;

 /// \short Function pointer to basis vectors function which sets  basis vectors
 /// b1 and b2 (which are in general functions of x)
 typedef void (*BasisVectorsFctPt) (const Vector<double>& x, Vector<double>& b1,
  Vector<double>& b2, DenseMatrix<double>& Db1, DenseMatrix<double>& Db2);

 /// \short enum to enumerate the possible edges that could be curved
 typedef typename MyC1CurvedElements::Edge Edge; 

 /// \short Get the pointer to the Curved shape class data member
 const MyC1CurvedElements::BernadouElementBasis<BOUNDARY_ORDER>* curved_shape_pt(){return &Curved_shape;};

 /// \short get the coordinate
 inline void get_coordinate_x(const Vector<double>& s, Vector<double>& x) const;

 /// \short get the coordinate i
 double interpolated_x (const Vector<double>s, const unsigned &i) const 
  { Vector<double> r(2); get_coordinate_x(s,r); return r[i]; }  

 /// \short Get the zeta coordinate
 inline void interpolated_zeta(const Vector<double> &s
  ,Vector<double> &zeta) const
 {
  /*
  // If there is a macro element use it
  if(this->Macro_elem_pt!=0) {this->get_x_from_macro_element(s,zeta);}
  */
  //Otherwise interpolate zeta_nodal using the shape functions
  /*else*/ {get_coordinate_x(s,zeta);}
 }

 // Upgrade an element to its curved counterpart
 inline void upgrade_to_curved_element(const Edge& curved_edge, const double& s_ubar,
  const double& s_obar,
  CurvilineGeomObject* parametric_edge);
 
 /// Lagrange interpolated shape for in--plane displacements
 void shape_u(const Vector<double> &s, Shape &psi) const;
  
 /// Lagrange interpolated d_shape for in--plane displacements
 void dshape_u_local(const Vector<double> &s,
                    Shape &psi, DShape &dpsids) const;

 /// Linear lagrange shape used for mapping  
 void Lshape(const Vector<double> &s, Shape &psi) const;
  
 /// Linear lagrange dshape used for mapping  
 void dLshape_local(const Vector<double> &s,
                    Shape &psi, DShape &dpsids) const;

 /// \short Get the Bell local to eulerian Jacobian determinant (based on 
 /// linear Lagrange shape)
 double J_eulerian1(const Vector<double> &s) const;

 /// \short Get the Bell local to eulerian Jacobian (based on linear 
 /// Lagrange shape)
 double local_to_eulerian_mapping2(const DShape &dpsids,
                                           DenseMatrix<double> &jacobian,
                                           DenseMatrix<double>
&inverse_jacobian) const;

  /// \short Get the bell shape (not basis) this is the linear Lagrange shape
  void my_interpolated_x(const Vector<double> &s, Vector<double> &x) const;

  ///\short Precompute the association matrix
  void precompute_association_matrix(DenseMatrix<double>& m)
   {
    // If the element has been upgraded
    if(Curved_edge ==MyC1CurvedElements::none)
     {} // Do nothing
    else 
     {Curved_shape.fill_in_full_association_matrix(m);}
   }; 
 
 /// Wrappers
 double n_basis_functions(){return Curved_shape.n_basis_functions();};
 double n_basic_basis_functions(){return Curved_shape.n_basic_basis_functions();};
protected:
 /// Get rotation matrices that change the degrees of freedom to the basis set
 /// by Rotated_basis_fct_pt
 inline void rotation_matrix_at_node(const unsigned& inode, DenseDoubleMatrix&
  rotation_matrix) const;
 
 /// \short Transform the shape functions so that they correspond to
 /// the new rotated dofs
 inline void rotate_shape(Shape& shape) const;

 /// \short Transform the shape functions and first derivatives so that they 
 /// correspond to the new rotated dofs
 inline void rotate_shape(Shape& shape, DShape& dshape) const;
  
 /// \short Transform the shape functions, first and second   derivatives so 
 /// that they correspond to the new rotated dofs
 inline void rotate_shape(Shape& shape, DShape& dshape,
   DShape& d2shape) const;
 
 /// \short Get the jth bubble dof at the lth internal point.
 inline double get_w_bubble_dof(const unsigned& l, const unsigned& j) const;

 /// \short Get the jth bubble dof at the lth internal point
 int local_w_bubble_equation(const unsigned& l, const unsigned& j) const;

public:
 ///\short  Constructor: Call constructors for C1CurvedBellElement and
 /// Biharmonic equations
 FoepplVonKarmanC1CurvedBellElement() :
  FoepplVonKarmanEquations<DIM,NNODE_1D>(), 
  Curved_edge(MyC1CurvedElements::none),
  Curved_shape(), Rotated_basis_fct_pt(0),  Nnodes_to_rotate(0)
  {
   this->set_nnodal_position_type(6);
   // Add the (zero) bubble dofs
   Bubble_w_internal_index = this->add_internal_data(new Data(0));
   #ifdef PARANOID
   Curved_edge_counter=0;
   #endif

  // Use the higher order integration scheme
  // delete this->integral_pt(); 
  TGauss<2,4>* new_integral_pt = new TGauss<2,4>;
  this->set_integration_scheme(new_integral_pt); 
  }

 ///\short  Destructor: clean up alloacations
 ~FoepplVonKarmanC1CurvedBellElement()
  {
   // May need a delete for this?
   //  Bubble_w_internal_index = this->add_internal_data(new Data(0));

   // Clean up allocation of integration scheme
   delete this->integral_pt(); 
  }

 /// Broken copy constructor
 FoepplVonKarmanC1CurvedBellElement(const
  FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>& dummy)
  {
   BrokenCopy::broken_copy("FoepplVonKarmanC1CurvedBellElement");
  }

 /// Broken assignment operator
 void operator=(const FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>&)
  {
   BrokenCopy::broken_assign("FoepplVonKarmanC1CurvedBellElement");
  }

 /// \short Set up the rotated degrees of freedom
 inline void set_up_rotated_dofs(const unsigned& nnodes_to_rotate ,
  const Vector<unsigned>& nodes_to_rotate, const BasisVectorsFctPt&
  basis_vectors_fct_pt);

 /// \short  Required  # of `values' (pinned or dofs)
 /// at node n
 inline unsigned required_nvalue(const unsigned &n) const
  {return Initial_Nvalue[n];}

 /// \short Output function:
 ///  x,y,u   or    x,y,z,u
 void output(std::ostream &outfile)
  {FoepplVonKarmanEquations<DIM,NNODE_1D>::output(outfile);}

 ///  \short Output function:
 ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
 void output(std::ostream &outfile, const unsigned &n_plot)
  {FoepplVonKarmanEquations<DIM,NNODE_1D>::output(outfile,n_plot);}

 ///  \short Output function:
 ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
 void output_interpolated_exact_soln(std::ostream &outfile,
  FiniteElement::SteadyExactSolutionFctPt exact_soln_pt, const unsigned& n_plot);

 /// \short C-style output function:
 ///  x,y,u   or    x,y,z,u
 void output(FILE* file_pt)
  {FoepplVonKarmanEquations<DIM,NNODE_1D>::output(file_pt);}

 ///  \short C-style output function:
 ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
 void output(FILE* file_pt, const unsigned &n_plot)
  {FoepplVonKarmanEquations<DIM,NNODE_1D>::output(file_pt,n_plot);}


 /// \short Output function for an exact solution:
 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {
   FoepplVonKarmanEquations<DIM,NNODE_1D>::output_fct(outfile,n_plot,
    exact_soln_pt);
  }

 /// \short Output function for a time-dependent exact solution.
 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
 /// (Calls the steady version)
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 const double& time,
                 FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
  {
   FoepplVonKarmanEquations<DIM,NNODE_1D>::output_fct(outfile,n_plot,time,
    exact_soln_pt);
  }


protected:

 /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
 void shape_and_test_foeppl_von_karman(const Vector<double> &s, Shape &psi,
  Shape& psi_b, Shape& test, Shape& test_b) const;

 /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
 inline double dshape_and_dtest_eulerian_foeppl_von_karman(const Vector<double> &s,
  Shape &psi, Shape &psi_b, DShape &dpsi_dx, DShape &dpsi_b_dx,
  Shape &test, Shape &test_b, DShape &dtest_dx, DShape &dtest_b_dx) const;

 /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
 inline double d2shape_and_d2test_eulerian_foeppl_von_karman(const Vector<double> &s,
  Shape &psi, Shape &psi_b, DShape &dpsi_dx, DShape &dpsi_b_dx,
  DShape &d2psi_dx2, DShape &d2psi_b_dx2,
  Shape &test, Shape &test_b, DShape &dtest_dx, DShape &dtest_b_dx,
  DShape &d2test_dx2, DShape &d2test_b_dx2) const;

// /// \short Shape, test functions & derivs. w.r.t. to global coords. at
// /// integration point ipt. Return Jacobian.
// inline double d2shape_and_d2test_eulerian_at_knot_foeppl_von_karman(const unsigned& ipt,
//                                                         Shape &psi,
//                                                         DShape &dpsidx,
//                                                         DShape &d2psidx,
//                                                         Shape &test,
//                                                         DShape &dtestdx,
//                                                         DShape &d2testdx)
//  const;
//
// inline double dshape_and_dtest_eulerian_at_knot_foeppl_von_karman(const unsigned &ipt,
//                                                         Shape &psi,
//                                                         DShape &dpsidx,
//                                                         Shape &test,
//                                                         DShape &dtestdx)
//  const;
 /// Get dshape and dtest for C0 (in--plane) displacements
 inline double dshape_u_and_dtest_u_eulerian_foeppl_von_karman(const Vector<double> &s, 
  Shape &psi, DShape &dpsidx, Shape &test, DShape &dtestdx) const;

// Private Data Members
private:

 #ifdef PARANOID
 /// Internal counter to check consistency
 unsigned Curved_edge_counter;
 #endif

 /// \short Static int that holds the number of variables at
 /// nodes: always the same
 static const unsigned Initial_Nvalue[];

 /// \short unsigned that holds the index the 'bubble' dofs of the element are
 /// stored
 unsigned Bubble_w_internal_index;

 ///  Which edge is curved, none by default
 Edge Curved_edge;

 /// Curved Shape function
 MyC1CurvedElements::BernadouElementBasis<BOUNDARY_ORDER> Curved_shape;

 /// A Pointer to the function that sets up the rotated basis at point x
 BasisVectorsFctPt Rotated_basis_fct_pt;

 /// Which nodes are we rotating
 Vector<unsigned> Nodes_to_rotate;

 /// Number of nodes to rotate
 unsigned Nnodes_to_rotate;


};




////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


//==============================================================================
/// Face geometry for the FoepplVonKarmanC1CurvedBellElement elements: The
/// spatial dimension of the face elements is one lower than that of the bulk
/// element but they have the same number of points along their 1D edges.
//==============================================================================
template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
class FaceGeometry<FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER> >:
 public virtual TElement<DIM-1,NNODE_1D>
{

  public:

 /// \short Constructor: Call the constructor for the
 /// appropriate lower-dimensional TElement
 FaceGeometry() : TElement<DIM-1,NNODE_1D>() {}

};



/////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//Inline functions:

//==============================================================================
/// Get the mapped position in the element. For straight sided elements this is
/// and affine mapping.
//==============================================================================
template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
void
FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::
get_internal_dofs_location(const unsigned& dof, Vector<double>& s) const
{
 // If the element has been upgraded
 if(Curved_edge ==MyC1CurvedElements::none)
  {
   throw OomphLibError(
    "There are no internal dofs for these elements as they have not been\
upgraded to curved elements.",
    OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }
 else 
  {Curved_shape.get_internal_dofs_location(dof,s);}
};

//==============================================================================
/// Get the mapped position in the element. For straight sided elements this is
/// and affine mapping.
//==============================================================================
template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
void FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::get_coordinate_x(const
 Vector<double>& s, Vector<double>& x) const
{
 // If the element has been upgraded
 if(Curved_edge ==MyC1CurvedElements::none)
  {this-> my_interpolated_x(s,x);}
 else 
  {Curved_shape.coordinate_x(s,x);}
};

//==============================================================================
/// Upgrade an element to its curved counterpart: this adds internal data to
/// elements and upgrades the shape class data member.
//==============================================================================
template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
inline void FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::upgrade_to_curved_element
 (const Edge& curved_edge, const double& s_ubar, const double& s_obar,
  CurvilineGeomObject* parametric_edge)
{
 #ifdef PARANOID
 // When upgrading add to count
 Curved_edge_counter+=1;
 // Check that we haven't upgraded this element already - we should introduce
 // a new function for this I think.
 if(Curved_edge_counter>1)
  {
   // SCREAM
   throw OomphLibError(
   "Cannot upgrade more than a single edge to be curved in C1 Curved Bell \
Elements.",OOMPH_CURRENT_FUNCTION,  OOMPH_EXCEPTION_LOCATION);
  }
 #endif
 using namespace MyC1CurvedElements;
 // Add the curved edge
 Curved_edge = curved_edge;

 // Set the integral pointer
 // HERE implement order 10 accuracy for faster 3rd order elements in clamped 
 // problems
 delete this->integral_pt();
 TGauss<2,13>* new_integral_pt = new TGauss<2,13>;
 this->set_integration_scheme(new_integral_pt); 

 // Set the number of internal dofs to 3
 this->Number_of_internal_dof_types = 1;
 this->Number_of_internal_dofs = Curved_shape.n_internal_dofs();
 Bubble_w_internal_index = this->add_internal_data(new Data(Curved_shape.n_internal_dofs()));

 // Set up the data of the element
 typename BernadouElementBasis<BOUNDARY_ORDER>::VertexList  vertices(3,Vector<double>(2,0.0));
 typename BernadouElementBasis<BOUNDARY_ORDER>::VertexList lvertices(3,Vector<double>(2,0.0));

 // Set up the local vertices
 lvertices[0][0]=1.0;
 lvertices[1][1]=1.0;

 // Now switch to upgrade
 // The shape functions are designed such that the curved edge is always edge
 // two. So this is where we set that up. This is temporary and not the final
 // solution we want
 switch(curved_edge)
  {
   // Throw an error if an edge is upgraded to none
   case none:
    throw OomphLibError( "Cannot upgrade edge 'none'. Curved elements must have\
one side defined by a parametric function.", OOMPH_CURRENT_FUNCTION,  
     OOMPH_EXCEPTION_LOCATION);
   break;
   case zero:
   // Everything cyclicly permutes
    for(unsigned i=0;i<2;++i)
     {
      vertices[2][i]=this->node_pt(0)->x(i);
      vertices[0][i]=this->node_pt(1)->x(i);
      vertices[1][i]=this->node_pt(2)->x(i);
     }
   break;
   case one:
   // Everything cyclicly permutes
    for(unsigned i=0;i<2;++i)
     {
      vertices[2][i]=this->node_pt(1)->x(i);
      vertices[0][i]=this->node_pt(2)->x(i);
      vertices[1][i]=this->node_pt(0)->x(i);
     }
   break;
   case two:
   // // Everything is just copied over
    for(unsigned i=0;i<2;++i)
     {
      vertices[2][i]=this->node_pt(2)->x(i);
      vertices[0][i]=this->node_pt(0)->x(i);
      vertices[1][i]=this->node_pt(1)->x(i);
     }
   break;
  }
 // Add the vertices to make the shape functions fully functional
 Curved_shape.upgrade_element(vertices, s_ubar,s_obar,curved_edge,*parametric_edge);
 // Now shift the nodes to be consistent with the new vertices
 unsigned n_node=this->nnode();
 // Only do the none--vertex nodes!
 for(unsigned l=0;l<n_node; ++l)
  {
   // Get local coordinate of node l
   Vector<double> s(2,0.0);
   Vector<double> x(2,0.0);
   this->local_coordinate_of_node(l,s);
   // Now use the curved mapping to reset the position
   Curved_shape.coordinate_x(s,x);
   // Now get the node
   Node* nod_pt=this->node_pt(l);
   // Reset its position
   nod_pt->x(0)=x[0];
   nod_pt->x(1)=x[1];
  }
}

//==============================================================================
/// Get the jth bubble dof at the lth internal point. Deliberately broken for
/// the case where there is no curved edge
//==============================================================================
template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
double FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::get_w_bubble_dof
 (const unsigned& l, const unsigned& j) const
{
  // Deliberately break this function for the below cases
  // If there is no curved edge then we cannot return anything meaningful
  if(Curved_edge==MyC1CurvedElements::none)
  {
  throw OomphLibError("There are no time-dependent internal 'bubble' dofs for \
this element.",OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  // Return dummy value 0.0
  return 0;
  }
  // For these elements we only have a single dof at each internal point
  else if(j!=0)
  {
  throw OomphLibError(
   "There is only a single degree of freedom at the internal points in this \
element.", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  // Return dummy value 0.0
  return 0;
  }
  // Now give the lth internal degree of freedom
  else
  {
   return this->internal_data_pt(Bubble_w_internal_index)->value(l);
  }
}

//==============================================================================
/// Get the jth bubble dof at the lth internal point. Deliberately broken for
/// case when there is no curved edge.
//==============================================================================
template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
int FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::local_w_bubble_equation(const
  unsigned& l, const unsigned& j) const
 {
  // Deliberately break this function for the below cases
  // If there is no curved edge then we cannot return anything meaningful
  if(Curved_edge==MyC1CurvedElements::none)
  {
  throw OomphLibError(
   "There are no time-dependent internal 'bubble' dofs for this element.",
    OOMPH_CURRENT_FUNCTION,  OOMPH_EXCEPTION_LOCATION);
   // Return dummy value -2
   return -2;
   }
   // For these elements we only have a single dof at each internal point
   else if(j!=0)
   {
   throw OomphLibError(
    "There is only a single degree of equation at the internal points in this \
element.", OOMPH_CURRENT_FUNCTION,OOMPH_EXCEPTION_LOCATION);
   // Return dummy value -2
   return -2;
   }
   // Now give the lth internal equation number
   else
   {
    return this->internal_local_eqn(Bubble_w_internal_index,l);
   }
  }

//==============================================================================
/// Set up the rotated degrees of freedom: includes a check for the number of
/// rotation nodes being greater than three.
//==============================================================================
template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
void FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::set_up_rotated_dofs(const unsigned&
  nnodes_to_rotate, const Vector<unsigned>& nodes_to_rotate, const
  BasisVectorsFctPt& basis_vectors_fct_pt)
{
 // Change the member Nnode_to_rotate
 Nnodes_to_rotate=nnodes_to_rotate;
 #ifdef PARANOID
  // Check that the number of nodes is smaller than 3
  if( nnodes_to_rotate > 3 )
   {
    throw OomphLibError(
     "There are only three nodes per element, so we cannot rotate more than\
 three ", OOMPH_CURRENT_FUNCTION,OOMPH_EXCEPTION_LOCATION);
   }
 #endif

 Nodes_to_rotate = nodes_to_rotate;

 // Point to the basis vectors function
 Rotated_basis_fct_pt = basis_vectors_fct_pt;
}

//==============================================================================
/// Rotate the shape functions according to
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//==============================================================================
template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
void FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::
 rotation_matrix_at_node (const unsigned& inode, DenseDoubleMatrix&
 rotation_matrix) const
{
 // Initialise x normal and tangent
 Vector<double> x(2,0.0);

 // Get the node pointer
 Node* nod_pt=this->node_pt(inode);

 // Get the position of the vertex
 x[0]=nod_pt->x(0);
 x[1]=nod_pt->x(1);

 // Initialise the two basis vectors
 Vector<Vector<double> > bi(2,Vector<double>(2,0.0));
 Vector<DenseMatrix<double> > Dbi(2,DenseMatrix<double>(2,2,0.0));

 // Now find the two new basis vectors
 // Get the normal - overload if not continuous
 (*Rotated_basis_fct_pt)(x,bi[0],bi[1], Dbi[0], Dbi[1]);

 // Rotation matrix, B
 DenseMatrix<double> b1(2,2,0.0),b22(3,3,0.0), b21(3,2,0.0);

 // Fill in the submatrices
 for(unsigned alpha=0; alpha<2; ++alpha)
  {
  for(unsigned beta=0; beta<2; ++beta)
   {
   // Fill in b1 - the Jacobian
   // Fill in the rotation of the first derivatives
   b1(alpha,beta) = bi[beta][alpha];
 
   // Avoid double counting the cross derivative
   if(alpha<=beta)
   {
    // Define row index
    const unsigned row = alpha + beta;
    for(unsigned gamma=0; gamma<2; ++gamma)
     {
     // Fill in b21 - the non affine part of the Jacobian derivative
     // Define column index
     unsigned col_b21 = gamma;
     // Fill in the non-affine part of the rotation of the second derivatives
     // if( beta>= alpha) ?
     b21(row,col_b21) += Dbi[gamma](alpha,beta);
     for(unsigned delta=0; delta<2; ++delta)
      {
      
      // Fill in b22 - the Affine part of the Jacobian derivative
      // Redefine column index for the next submatrix
      unsigned col_b22 = gamma + delta;
      // Fill in the affine part of the rotation of the second derivatives
      // if( beta>= alpha) ?
      b22(row,col_b22) += bi[gamma][alpha] * bi[delta][beta];
      }
     }
    }
   }
 }


 // Fill in the submatrices to the full (6x6) matrix - we need to right multiply
 // this matrix so we need the transpose of the Jacobian
 // W dof remains the same
 rotation_matrix(0,0)=1.0;
 // Fill in b1
 for(unsigned i=0; i<2; ++i){
  for(unsigned j=0; j<2; ++j)
   { rotation_matrix(1+j,1+i) = b1(i,j); }
  }
 // Fill in b21
 for(unsigned i=0; i<3; ++i){
  for(unsigned j=0; j<2; ++j)
   { rotation_matrix(1+j,3+i) = b21(i,j); }
  }
 // Fill in b22
 for(unsigned i=0; i<3; ++i){
  for(unsigned j=0; j<3; ++j)
   { rotation_matrix(3+j,3+i) = b22(i,j); }
  }
}


//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
 void FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::shape_and_test_foeppl_von_karman(
  const Vector<double> &s, Shape &psi, Shape& psi_b,  Shape &test, Shape& test_b
  ) const
{
 throw OomphLibError(
 "This still needs testing for curved elements.",
 "void FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::\
shape_and_test_foeppl_von_karman(...)", OOMPH_EXCEPTION_LOCATION); // HERE

 // Get dummy shape functions for the Bell call
 DShape dpsidx(3,6,2);
 DShape d2psidx(3,6,3);

 // Vertices
 Vector<Vector<double> > v(3,Vector<double>(2));
 for (unsigned inode=0;inode<3;++inode)
  {
   // Get the position vector
   Node* nod_pt=this->node_pt(inode);
   v[inode][0]=nod_pt->x(0);
   v[inode][1]=nod_pt->x(1);
  }

 // If the element has not been upgraded
 if(Curved_edge ==MyC1CurvedElements::none)
  { 
  // Get J
  this->J_eulerian1(s);
  MyShape::d2_basis_eulerian(s,v,psi,dpsidx,d2psidx);
  }
 else // i.e if has curved edge
  {Curved_shape.shape(s,psi,psi_b);}
 
 // Rotate the degrees of freedom
 rotate_shape(psi);
 
 // Galerkin
 // (Shallow) copy the basis functions
 test = psi;
 test_b = psi_b;

}


//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
 double FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::
 dshape_and_dtest_eulerian_foeppl_von_karman(const Vector<double> &s, Shape &psi,
 Shape& psi_b, DShape &dpsidx, DShape& dpsi_b_dx,  Shape &test, Shape& test_b,
 DShape &dtestdx,DShape &dtest_b_dx) const
{
 // Throw if called 
 throw OomphLibError(
 "This still needs testing for curved elements.",
 "void FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::\
dshape_and_dtest_foeppl_von_karman(...)",
  OOMPH_EXCEPTION_LOCATION);// HERE

 // Now set up dummy DShape so we can call Bell
 DShape d2psidx(3,6,3);
 double J =this->J_eulerian1(s);

 // Vertices
 Vector<Vector<double> > v(3,Vector<double>(2));
 for (unsigned inode=0;inode<3;++inode)
  {
   // Get the position vector
   Node* nod_pt=this->node_pt(inode);
   v[inode][0]=nod_pt->x(0);
   v[inode][1]=nod_pt->x(1);
  }

 // If the element has been upgraded
 if(Curved_edge ==MyC1CurvedElements::none)
  { 
  // Get J
  J=this->J_eulerian1(s);
  MyShape::d2_basis_eulerian(s,v,psi,dpsidx,d2psidx);
  }
 else // i.e if has curved edge
  {J=Curved_shape.d_shape_dx(s,psi,psi_b,dpsidx,dpsi_b_dx);}

 // Rotate the degrees of freedom
 rotate_shape(psi, dpsidx);
 // Galerkin
 // (Shallow) copy the basis functions
 test = psi;
 dtestdx= dpsidx;
 test_b = psi_b;
 dtest_b_dx= dpsi_b_dx;

 return J;
}

template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
 double FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::
  d2shape_and_d2test_eulerian_foeppl_von_karman(const Vector<double> &s,  Shape &psi,
  Shape &psi_b, DShape &dpsidx, DShape &dpsi_bdx,  DShape &d2psidx,
  DShape &d2psi_bdx,
  Shape &test, Shape &test_b, DShape &dtestdx, DShape &dtest_bdx,
  DShape &d2testdx, DShape &d2test_bdx) const
{
 //Call the geometrical shape functions and derivatives
 double J=this->J_eulerian1(s);
 // Vertices
 Vector<Vector<double> > v(3,Vector<double>(2));
 for (unsigned inode=0;inode<3;++inode)
  {
   // Get the position vector
   Node* nod_pt=this->node_pt(inode);
   v[inode][0]=nod_pt->x(0);
   v[inode][1]=nod_pt->x(1);
  }

 // If the element has been upgraded
 if(Curved_edge ==MyC1CurvedElements::none)
  { 
  // Get J
  J=this->J_eulerian1(s);
  MyShape::d2_basis_eulerian(s,v,psi,dpsidx,d2psidx);
  }
 else if(this->get_association_matrix_pt() != 0)// i.e if has curved edge and precomputed matrix
  {
   // Temporary
  // DenseMatrix<double> m=*(this->get_association_matrix_pt());
  J=Curved_shape.d2_shape_dx2(s,psi,psi_b,dpsidx,dpsi_bdx,d2psidx,d2psi_bdx,*(this->get_association_matrix_pt()));
  }
 else 
  {J=Curved_shape.d2_shape_dx2(s,psi,psi_b,dpsidx,dpsi_bdx,d2psidx,d2psi_bdx);}
 // Rotate the dofs
 rotate_shape(psi, dpsidx, d2psidx);
 
 // Galerkin
 //Set the test functions equal to the shape functions (this is a shallow copy)
 test = psi;
 dtestdx= dpsidx;
 d2testdx = d2psidx;
 test_b = psi_b;
 dtest_bdx= dpsi_bdx;
 d2test_bdx = d2psi_bdx;

 //Return the jacobian
 return J;
}

//======================================================================
/// Rotate the shape functions according to the specified basis on the 
/// boundary. This function does a DenseDoubleMatrix solve to determine 
/// new basis - which could be speeded up by caching the matrices higher
/// up and performing the LU decomposition only once
//======================================================================
template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
inline void FoepplVonKarmanC1CurvedBellElement<DIM, NNODE_1D,BOUNDARY_ORDER>::
 rotate_shape(Shape& psi) const
{
 // Loop over the nodes with rotated dofs
 for(unsigned i=0; i<Nnodes_to_rotate; ++i)
  {
   // Get the nodes
   unsigned inode = Nodes_to_rotate[i];

   // Construct the vectors to hold the shape functions
   Vector<double> psi_vector(6);

   // Get the rotation matrix
   DenseDoubleMatrix  rotation_matrix(6,6,0.0);
   this->rotation_matrix_at_node(inode,rotation_matrix);

   // Copy to the vectors
   for(unsigned l=0;l<6;++l)
    {
     // Copy to the vectors
     for(unsigned k=0;k<6;++k)
      {
      // Copy over shape functions
     // psi_vector[l]=psi(inode,l);
      psi_vector[l]+=psi(inode,k)*rotation_matrix(l,k); 
      }
    }

   // Copy back to shape the rotated vetcors
   for(unsigned l=0;l<6;++l)
    {
      // Copy over shape functions
      psi(inode,l)=psi_vector[l];
    }
  }
}

//======================================================================
/// Rotate the shape functions according to the specified basis on the 
/// boundary. This function does a DenseDoubleMatrix solve to determine 
/// new basis - which could be speeded up by caching the matrices higher
/// up and performing the LU decomposition only once
//======================================================================
template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
inline void FoepplVonKarmanC1CurvedBellElement<DIM, NNODE_1D,BOUNDARY_ORDER>::
 rotate_shape(Shape& psi, DShape& dpsidx) const
{
 // Loop over the nodes with rotated dofs
 for(unsigned i=0; i<Nnodes_to_rotate; ++i)
  {
   // Get the nodes
   unsigned inode = Nodes_to_rotate[i];

   // Construct the vectors to hold the shape functions
   Vector<double> psi_vector(6);
   Vector<Vector<double> > dpsi_vector_dxi(2,Vector<double>(6));

   // Get the rotation matrix
   DenseDoubleMatrix  rotation_matrix(6,6,0.0);
   this->rotation_matrix_at_node(inode,rotation_matrix);

   // Copy to the vectors
   for(unsigned l=0;l<6;++l)
    {
     // Copy to the vectors
     for(unsigned k=0;k<6;++k)
      {
      // Copy over shape functions
     // psi_vector[l]=psi(inode,l);
      psi_vector[l]+=psi(inode,k)*rotation_matrix(l,k); 
      // Copy over first derivatives
      for(unsigned i=0;i<2; ++i)
       {dpsi_vector_dxi[i][l]+=dpsidx(inode,k,i)*rotation_matrix(l,k);}
      }
    }

   // Copy back to shape the rotated vetcors
   for(unsigned l=0;l<6;++l)
    {
      // Copy over shape functions
      psi(inode,l)=psi_vector[l];
      // Copy over first derivatives
      for(unsigned i=0;i<2; ++i)
        {dpsidx(inode,l,i)=dpsi_vector_dxi[i][l];}
    }
  }
}

//======================================================================
/// Rotate the shape functions according to the specified basis on the 
/// boundary. This function does a DenseDoubleMatrix solve to determine 
/// new basis - which could be speeded up by caching the matrices higher
/// up and performing the LU decomposition only once
//======================================================================
template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
inline void FoepplVonKarmanC1CurvedBellElement<DIM, NNODE_1D, BOUNDARY_ORDER>::
 rotate_shape(Shape& psi, DShape& dpsidx, DShape& d2psidx) const
{
 // Loop over the nodes with rotated dofs
 for(unsigned i=0; i<Nnodes_to_rotate; ++i)
  {
   // Get the nodes
   unsigned inode = Nodes_to_rotate[i];

   // Construct the vectors to hold the shape functions
   Vector<double> psi_vector(6);
   Vector<Vector<double> > dpsi_vector_dxi(2,Vector<double>(6));
   Vector<Vector<double> > d2psi_vector_dxidxj(3,Vector<double>(6));

   // Get the rotation matrix
   DenseDoubleMatrix  rotation_matrix(6,6,0.0);
   this->rotation_matrix_at_node(inode,rotation_matrix);

   // Copy to the vectors
   for(unsigned l=0;l<6;++l)
    {
     // Copy to the vectors
     for(unsigned k=0;k<6;++k)
      {
      // Copy over shape functions
     // psi_vector[l]=psi(inode,l);
      psi_vector[l]+=psi(inode,k)*rotation_matrix(l,k); 
      // Copy over first derivatives
      for(unsigned i=0;i<2; ++i)
       {dpsi_vector_dxi[i][l]+=dpsidx(inode,k,i)*rotation_matrix(l,k);}
      for(unsigned i=0;i<3; ++i)
       {d2psi_vector_dxidxj[i][l]+=d2psidx(inode,k,i)*rotation_matrix(l,k);}
      }
    }

   // Copy back to shape the rotated vetcors
   for(unsigned l=0;l<6;++l)
    {
      // Copy over shape functions
      psi(inode,l)=psi_vector[l];
      // Copy over first derivatives
      for(unsigned i=0;i<2; ++i)
        {dpsidx(inode,l,i)=dpsi_vector_dxi[i][l];}
      // Copy over second derivatives
      for(unsigned i=0;i<3; ++i)
        {d2psidx(inode,l,i)=d2psi_vector_dxidxj[i][l];}
    }
  }
}

//============================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//=============================================================================
template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
 double FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::
  dshape_u_and_dtest_u_eulerian_foeppl_von_karman(const Vector<double> &s,Shape &psi,
  DShape &dpsidx,  Shape &test,  DShape &dtestdx) const
{
 // Initialise
 const unsigned n_node = this->nnode();
 DShape dpsids(n_node,DIM);
 dshape_u_local(s,psi,dpsids);
 double J;

 Shape shape(3);
 DShape dshape_ds(3,DIM);
 this->dLshape_local(s,shape, dshape_ds);

 // Get the inverse jacobian
 // If the element has been upgraded
 DenseMatrix<double> jacobian(DIM,DIM,0.0),inverse_jacobian(DIM,DIM,0.0);
 if(Curved_edge ==MyC1CurvedElements::none)
  { 
   // Get the Affine jacobian
   J = this->local_to_eulerian_mapping2(dshape_ds,jacobian,
   inverse_jacobian);
   // Now find the global derivatives
    for(unsigned l=0;l<n_node;++l)
    { 
      for(unsigned i=0; i<DIM; ++i)
       {
        // Initialise to zero
        dpsidx(l,i)=0.0;
        for(unsigned j=0; j<DIM; ++j)
         {
          // NB MyJ = oomphJ' consider fixing?
          // Convert to local coordinates
          dpsidx(l,i)+=inverse_jacobian(i,j)*dpsids(l,j);
         }
       }
    }
  }
 else // i.e if has curved edge
  {
   Curved_shape.get_jacobian(s,jacobian);
   J = MyShape::invert_two_by_two(jacobian,inverse_jacobian);
   // Now find the global derivatives
    for(unsigned l=0;l<n_node;++l)
    { 
      for(unsigned i=0; i<DIM; ++i)
       {
        // Initialise to zero
        dpsidx(l,i)=0.0;
        for(unsigned j=0; j<DIM; ++j)
         {
          // NB MyJ = oomphJ' consider fixing?
          // Convert to local coordinates
          dpsidx(l,i)+=dpsids(l,j)*inverse_jacobian(j,i);
         }
       }
    }
  }

 test = psi;
 dtestdx = dpsidx;

 return J;
}

template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
 void FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::pin_all_deflection_dofs() const
 {
  // CHECK HERE
  // Curved Bell elements only have deflection dofs at vertices
  for(unsigned n=0; n<3; ++n)
   {
    // Get node
    Node* nod_pt=this->node_pt(n);
    // Check if it is on the boundary
    for(unsigned i=0;i<6;++i)
     {
      // Pin and set the value
      nod_pt->pin(2+i);
      nod_pt->set_value(2+i,0.0);
     }
   }

  // Get number of internal dofs
  const unsigned n_b_node=this->Number_of_internal_dofs;
  // const unsigned n_b_position_type = this->Number_of_internal_dof_types;
  // Now fix internal dofs
  for(unsigned n=0; n<n_b_node; ++n)
   {
    // Get node
    // Pin and set the value
   // for(unsigned l=0;l<n_b_position_type;++l)
   //  {
      this->internal_data_pt(Bubble_w_internal_index)->pin(n);
   //  }
   }
 }


template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
 void FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::
 fix_out_of_plane_displacement_dof(const unsigned& dof_number, const unsigned&
b,const DisplacementFctPt& w)
 {
  // CHECK HERE
  // Bell elements only have deflection dofs at vertices
  for(unsigned n=0; n<3; ++n)
   {
    // Get node
    Node* nod_pt=this->node_pt(n);
    // Check if it is on the boundary
    bool is_boundary_node=nod_pt->is_on_boundary(b);
    if(is_boundary_node)
     {
      // Extract nodal coordinates from node:
      Vector<double> x(2);
      x[0]=nod_pt->x(0);
      x[1]=nod_pt->x(1);
      // Get value
      double value;
      w(x,value);
      // Pin and set the value
      nod_pt->pin(dof_number);
      nod_pt->set_value(dof_number,value);
     }
   }
 }

template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
 void FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::
 fix_in_plane_displacement_dof(const unsigned& dof_number, const unsigned& b, const
DisplacementFctPt& u)
 {
  // CHECK HERE
   const unsigned n_node = this->nnode();
  // Bell elements only have deflection dofs at vertices
  for(unsigned n=0; n<n_node; ++n)
   {
    // Get node
    Node* nod_pt=this->node_pt(n);
    // Check if it is on the boundary
    bool is_boundary_node=nod_pt->is_on_boundary(b);
    if(is_boundary_node)
     {
      // Extract nodal coordinates from node:
      // HERE shouldn't we get s of node and then use interpolated position
      // here to be consistent?
      // Since the element isn't necessarily isoparametric the nodes 'position'
      // is not necessarily correct?
      Vector<double> x(2),s(2);
      this->local_coordinate_of_node(n,s);
      get_coordinate_x(s, x); 
      /*
      x[0]=nod_pt->x(0);
      x[1]=nod_pt->x(1);
      */
      // Get value
      double value;
      u(x,value);
      // Pin and set the value
      nod_pt->pin(dof_number);
      nod_pt->set_value(dof_number,value);
     }
   }

 }
////==============================================================================
///// Define the shape functions and test functions and derivatives
///// w.r.t. global coordinates and return Jacobian of mapping.
/////
///// Galerkin: Test functions = shape functions
////==============================================================================
//template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
//double FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::
// dshape_and_dtest_eulerian_at_knot_foeppl_von_karman(
//  const unsigned &ipt,
//  Shape &psi,
//  DShape &dpsidx,
//  Shape &test,
//  DShape &dtestdx) const
//{
// const double J = this->dshape_and_dtest_eulerian_at_knot(ipt,psi,dpsidx,test
//  ,dtestdx);
//
// return J;
//}
//
//template<unsigned DIM, unsigned NNODE_1D>
//double FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::
// d2shape_and_d2test_eulerian_at_knot_foeppl_von_karman(
//  const unsigned &ipt,
//  Shape &psi,
//  DShape &dpsidx,
//  DShape &d2psidx,
//  Shape &test,
//  DShape &dtestdx,
//  DShape &d2testdx) const
//{
// //Call the geometrical shape functions and derivatives
// const double J = this->d2shape_and_d2test_eulerian_at_knot(ipt,psi,dpsidx,
//  d2psidx,test,dtestdx,d2testdx);
//
// //Return the jacobian
// return J;
//}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//=======================================================================
/// Shape function for specific TElement<2,2>
//=======================================================================
template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
  void FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::Lshape(const Vector<double> &s, Shape &psi) const
   {
    psi[0] = s[0];
    psi[1] = s[1];
    psi[2] = 1.0-s[0]-s[1];
   }
  

//=======================================================================
/// Derivatives of shape functions for specific TElement<2,2>
//=======================================================================
template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
  void FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::dLshape_local(const Vector<double> &s,
                    Shape &psi, DShape &dpsids) const
   {
    this->Lshape(s, psi);
    
    // Derivatives
    dpsids(0,0) = 1.0;
    dpsids(0,1) = 0.0;
    dpsids(1,0) = 0.0;
    dpsids(1,1) = 1.0;
    dpsids(2,0) = -1.0;
    dpsids(2,1) = -1.0;
   }
  


template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
double FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::local_to_eulerian_mapping2(const DShape &dpsids,
                                           DenseMatrix<double> &jacobian,
                                           DenseMatrix<double> &inverse_jacobian) const
  {
    //Assemble the jacobian
    //Locally cache the elemental dimension
   const unsigned el_dim = DIM;
   //The number of shape functions must be equal to the number
   //of nodes (by definition)
   const unsigned n_shape = 3;
   //The number of shape function types must be equal to the number
   //of nodal position types (by definition)
   const unsigned n_shape_type = 1;
  
   //Loop over the rows of the jacobian
   for(unsigned i=0;i<el_dim;i++)
    {
     //Loop over the columns of the jacobian
     for(unsigned j=0;j<el_dim;j++)
      {
       //Zero the entry
       jacobian(i,j) = 0.0;
       //Loop over the shape functions
       for(unsigned l=0;l<n_shape;l++)
        {
         for(unsigned k=0;k<n_shape_type;k++)
          {
           //Jacobian is dx_j/ds_i, which is represented by the sum
           //over the dpsi/ds_i of the nodal points X j
           //Call the Non-hanging version of positions
           //This is overloaded in refineable elements
           jacobian(i,j) += this->raw_nodal_position_gen(l,k,j)*dpsids(l,k,i);
          }
        }
      }
    }
  
   //Invert the jacobian (use the template-free interface)
   return this->invert_jacobian_mapping(jacobian,inverse_jacobian);
  }

 //========================================================================
 /// \short Calculate the Jacobian of the mapping between local and global
 /// coordinates at the position s
 //========================================================================
template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
 double FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::J_eulerian1(const Vector<double> &s) const
  {
   //Find the number of nodes and position types
   const unsigned n_node = 3;
   const unsigned n_position_type = 1;
   //Find the dimension of the node and element
   const unsigned n_dim_node = DIM;//nodal_dimension();
   const unsigned n_dim_element = DIM;
   
   
   //Set up dummy memory for the shape functions
   Shape psi(n_node,n_position_type);
   DShape dpsids(n_node,n_position_type,n_dim_element);
   //Get the shape functions and local derivatives
   this->dLshape_local(s,psi,dpsids);
   
   //Right calculate the base vectors
   DenseMatrix<double> interpolated_G(n_dim_element,n_dim_node);
   
   //Loop over the dimensions and compute the entries of the 
   //base vector matrix
   for(unsigned i=0;i<n_dim_element;i++)
    {
     for(unsigned j=0;j<n_dim_node;j++)
      {
       //Initialise the j-th component of the i-th base vector to zero
       interpolated_G(i,j) = 0.0;
       for(unsigned l=0;l<n_node;l++)
        {
         for(unsigned k=0;k<n_position_type;k++)
          {
           interpolated_G(i,j) += this->raw_nodal_position_gen(l,k,j)*dpsids(l,k,i);
          }
        }
      }
    }

   //Calculate the metric tensor of the element
   DenseMatrix<double> G(n_dim_element);
   for(unsigned i=0;i<n_dim_element;i++)
    {
     for(unsigned j=0;j<n_dim_element;j++)
      {
       //Initialise to zero
       G(i,j) = 0.0;
       for(unsigned k=0;k<n_dim_node;k++) 
        {
         G(i,j) += interpolated_G(i,k)*interpolated_G(j,k);
        }
      }
    }
   
   //Calculate the determinant of the metric tensor
   double det = 0.0;
   switch(n_dim_element)
    {
    case 0:
     throw OomphLibError("Cannot calculate J_eulerian() for point element\n",
                         "FiniteElement::J_eulerian()",
                         OOMPH_EXCEPTION_LOCATION);
     break;
    case 1:
     det = G(0,0);
     break;
    case 2:
     det = G(0,0)*G(1,1) - G(0,1)*G(1,0);
     break;
    case 3:
     det = G(0,0)*G(1,1)*G(2,2) + G(0,1)*G(1,2)*G(2,0) + G(0,2)*G(1,0)*G(2,1)
      - G(0,0)*G(1,2)*G(2,1) - G(0,1)*G(1,0)*G(2,2) - G(0,2)*G(1,1)*G(2,0);
     break;
    default:
     oomph_info << "More than 3 dimensions in J_eulerian()" << std::endl;
     break;
    }
   
#ifdef PARANOID
   this->check_jacobian(det);
#endif
   
   //Return the Jacobian (square-root of the determinant of the metric tensor)
   return sqrt(det);
   
  }

//=======================================================================
/// Return FE interpolated position x[] at local coordinate s as Vector
// (using linear Lagrange interpolation)
//=======================================================================
template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
  void FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::my_interpolated_x(const Vector<double> &s, Vector<double> &x)
   const
  {
   //Find the number of nodes
   const unsigned n_node = 3;
   //Find the number of positional types
   const unsigned n_position_type = 1;
   
   //Find the dimension stored in the node
   const unsigned nodal_dim = DIM;
   
   //Assign storage for the local shape function
   Shape psi(n_node);
   //Find the values of shape function
   Lshape(s,psi);
   
   //Loop over the dimensions
   for(unsigned i=0;i<nodal_dim;i++)
    {
     //Initilialise value of x[i] to zero
     x[i] = 0.0;
     //Loop over the local nodes
     for(unsigned l=0;l<n_node;l++) 
      {
       for(unsigned k=0;k<n_position_type;k++)
        {
         x[i] += this->nodal_position_gen(l,k,i)*psi[l];
        }
      }
    }
  }
}


#endif
