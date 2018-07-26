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
//Head file for the Biharmonic Bell elements
#ifndef OOMPH_BIHARMONIC_BELL_ELEMENTS_HEADER
#define OOMPH_BIHARMONIC_BELL_ELEMENTS_HEADER

#include "foeppl_von_karman_elements.h"
#include "MyBellShape.h"

namespace oomph
{

//======================================================================
/// FoepplVonKarmanBellElement elements are a subparametric scheme with
/// linear Lagrange interpolation for approximating the geometry and
/// the C1-functions for approximating variables.
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
class FoepplVonKarmanBellElement : public virtual FoepplVonKarmanEquations<DIM,NNODE_1D>
{
public:
 /// \short Function pointer to pressure function fct(x,f(x)) --
 /// x is a Vector!
 typedef void (*DisplacementFctPt)(const Vector<double>& x, double& f);

 /// \short function to pin all deflection dofs
 void pin_all_deflection_dofs() const;

 /// \short Function pointer to basis vectors function which sets  basis vectors
 /// b1 and b2 (which are in general functions of x)
 typedef void (*BasisVectorsFctPt) (const Vector<double>& x, Vector<double>& b1,
  Vector<double>& b2);

 /// \short function to pin particular out--of--plane displacement dofs
 void fix_out_of_plane_displacement_dof(const unsigned& dof_number, const
unsigned& b, const DisplacementFctPt& w);

 /// \short function to pin particular in--plane displacement dofs
 void fix_in_plane_displacement_dof(const unsigned& dof_number, const unsigned&b,
  const DisplacementFctPt& u);

 // Quick, lets rob some of WPs code! 
 void shape_u(const Vector<double> &s, Shape &psi) const;
 
 void Lshape(const Vector<double> &s, Shape &psi) const;
  
 void dshape_u_local(const Vector<double> &s,
                    Shape &psi, DShape &dpsids) const;

 void dLshape_local(const Vector<double> &s,
                    Shape &psi, DShape &dpsids) const;
 double J_eulerian1(const Vector<double> &s) const;

 double local_to_eulerian_mapping2(const DShape &dpsids,
                                           DenseMatrix<double> &jacobian,
                                           DenseMatrix<double>
&inverse_jacobian) const;

  void my_interpolated_x(const Vector<double> &s, Vector<double> &x) const;


  void precompute_association_matrix(DenseMatrix<double>& m){}; //Do nothing
  // Not used
  double n_basis_functions(){return 18;};
  double n_basic_basis_functions(){return 18;};
private:

 /// \short Static int that holds the number of variables at
 /// nodes: always the same
 static const unsigned Initial_Nvalue[];

protected:
 /// A Pointer to the function that sets up the rotated basis at point x
 BasisVectorsFctPt Rotated_basis_fct_pt;

 /// Which nodes are we rotating
 Vector<unsigned> Nodes_to_rotate;

 // Number of nodes to rotate
 unsigned Nnodes_to_rotate;
 /// Get rotation matrices that change the degrees of freedom to the basis set
 /// by Rotated_basis_fct_pt
 inline void rotation_matrix_at_node(const unsigned& inode, DenseMatrix<double>&
  rotation_matrix) const;

 double get_w_bubble_dof(const unsigned& l, const unsigned& j) const
  {
   throw OomphLibError(
    "There are no time-dependent internal 'bubble' dofs for these elements",
    OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
   // Return dummy value 0.0
   return 0;
  }
 /// Get the jth bubble dof at the lth internal point
 int local_w_bubble_equation(const unsigned& l, const unsigned& j) const
  {
   // If there is no curved edge then we cannot return anything meaningful
   throw OomphLibError(
    "There are no time-dependent internal 'bubble' dofs for this element.",
    OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
   // Return dummy value -2
   return -2;
  }

public:

 /// \short get the coordinate
 void get_coordinate_x(const Vector<double>& s, Vector<double>& x) const
  {
   // Just use WPs get_x
  this->my_interpolated_x(s,x);
  };

 ///\short  Constructor: Call constructors for BellElement and
 /// Biharmonic equations
 FoepplVonKarmanBellElement() :
  FoepplVonKarmanEquations<DIM,NNODE_1D>(), Rotated_basis_fct_pt(0),
  Nnodes_to_rotate(0)
  {
   this->set_nnodal_position_type(6);
   // Use the higher order integration scheme
   TGauss<2,5>* new_integral_pt = new TGauss<2,5>;
   this->set_integration_scheme(new_integral_pt); 
  }

 ///\short Destructor
 ~FoepplVonKarmanBellElement() 
  {
   // Use the higher order integration scheme
   delete this->integral_pt(); 
  }

 /// Broken copy constructor
 FoepplVonKarmanBellElement(const FoepplVonKarmanBellElement<DIM,NNODE_1D>& dummy)
  {
   BrokenCopy::broken_copy("FoepplVonKarmanBellElement");
  }

 /// Broken assignment operator
 void operator=(const FoepplVonKarmanBellElement<DIM,NNODE_1D>&)
  {
   BrokenCopy::broken_assign("FoepplVonKarmanBellElement");
  }

  /// Set up the rotated degrees of freedom
  inline void set_up_rotated_dofs(const unsigned& nnodes_to_rotate ,
   const Vector<unsigned>& nodes_to_rotate, const BasisVectorsFctPt&
   basis_vectors_fct_pt)
  {
   // Change the member Nnode_to_rotate
   Nnodes_to_rotate=nnodes_to_rotate;
   #ifdef PARANOID
    // Check that the number of nodes is smaller than 3
    if( nnodes_to_rotate > 3 )
     {
      throw OomphLibError(
       "There are only three nodes per element, so we cannot rotate more than\
 three ", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
     }
   #endif

   Nodes_to_rotate = nodes_to_rotate;

   // Point to the basis vectors function
   Rotated_basis_fct_pt = basis_vectors_fct_pt;
  }


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
  {FoepplVonKarmanEquations<DIM,NNODE_1D>::output_fct(outfile,n_plot,exact_soln_pt);}



 /// \short Output function for a time-dependent exact solution.
 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
 /// (Calls the steady version)
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 const double& time,
                 FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
  {FoepplVonKarmanEquations<DIM,NNODE_1D>::output_fct(outfile,n_plot,time,exact_soln_pt);}


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

 inline double dshape_u_and_dtest_u_eulerian_foeppl_von_karman(const Vector<double> &s, 
  Shape &psi, DShape &dpsidx, Shape &test, DShape &dtestdx) const;
};




////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


//=======================================================================
/// Face geometry for the FoepplVonKarmanBellElement elements: The spatial
/// dimension of the face elements is one lower than that of the
/// bulk element but they have the same number of points
/// along their 1D edges.
//=======================================================================
template<unsigned DIM, unsigned NNODE_1D>
class FaceGeometry<FoepplVonKarmanBellElement<DIM,NNODE_1D> >:
 public virtual TElement<DIM-1,NNODE_1D>
{

  public:

 /// \short Constructor: Call the constructor for the
 /// appropriate lower-dimensional TElement
 FaceGeometry() : TElement<DIM-1,NNODE_1D>() {}

};



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

//Inline functions:

//======================================================================
/// Rotate the shape functions according to
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned DIM, unsigned NNODE_1D>
void FoepplVonKarmanBellElement<DIM,NNODE_1D>::rotation_matrix_at_node
 (const unsigned& inode, DenseMatrix<double>& rotation_matrix) const
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

 // Now find the two new basis vectors
 // Get the normal - overload if not continuous
 (*Rotated_basis_fct_pt)(x,bi[0],bi[1]);
 // Rotate first derivatives
 DenseMatrix<double> b1(2,2,0.0),ib1(2,2,0.0), b2(3,3,0.0),ib2(3,3,0.0);

 // This function relies upon a particular ordering of the Hermite degrees of
 // freedom - It will no longer work if the shape functions are modified such
 // that the ordering of dofs chanages
 // No change to first dof
  rotation_matrix(0,0) = 1;

 // Rotate first derivative dofs, loop over directions and basis vectors
 b1(0,0) = bi[0][0];
 b1(0,1) = bi[1][0];
 b1(1,0) = bi[0][1];
 b1(1,1) = bi[1][1];

 // Get inverse ib1
 MyShape::invert_two_by_two(b1,ib1);

 // Rotate second derivative dofs, loop over directions and basis vectors
 b2(0,0) = bi[0][0]*bi[0][0];
 b2(1,0) = bi[0][1]*bi[0][0];
 b2(2,0) = bi[0][1]*bi[0][1];

 b2(0,1) = bi[1][0]*bi[1][0];
 b2(1,1) = bi[1][1]*bi[1][0];
 b2(2,1) = bi[1][1]*bi[1][1];

 b2(0,2) = bi[0][0]*bi[1][0];
 b2(1,2) = (bi[0][1]*bi[1][0]+bi[1][1]*bi[0][0]);
 b2(2,2) = bi[0][1]*bi[1][1];

 // Get inverse ib1
 MyShape::invert_three_by_three(b2,ib2);

 //Fill in the rotation matrix
 for(unsigned i=0;i<2;i++)
  {
   for(unsigned j=0;j<2;j++)
    {
     rotation_matrix(i+1,j+1)=ib1(i,j);
    }
  }

 //Fill in the rotation matrix
 for(unsigned i=0;i<3;i++)
  {
   for(unsigned j=0;j<3;j++)
    {
     rotation_matrix(i+3,j+3)=ib2(i,j);
    }
  }
}


//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned DIM, unsigned NNODE_1D>
 void FoepplVonKarmanBellElement<DIM,NNODE_1D>::shape_and_test_foeppl_von_karman(
  const Vector<double> &s, Shape &psi, Shape& psi_b,  Shape &test, Shape& test_b
  ) const
{
 //const double J = this->dbasis_eulerian(s,psi,dpsidx);
 DShape dpsidx(3,6,3);
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

 MyShape::d2_basis_eulerian(s,v,psi,dpsidx,d2psidx);
 // Loop over nodes with rotated dofs
 // Hierher?
 for(unsigned i=0; i<Nnodes_to_rotate; ++i)
  {
   // Get the nodes
   unsigned inode = Nodes_to_rotate[i];
   // Construct the new psi
   Vector<double> psi_new(6,0.0);
   DenseMatrix<double> rotation_matrix(6,6,0.0);
   this->rotation_matrix_at_node(inode,rotation_matrix);

   // Now loop over the shape functions at node i
   for(unsigned l=0;l<6;++l)
    {
     for(unsigned k=0;k<6;++k)
     {
      // Now rotate the shape functions
      psi_new[l]+=rotation_matrix(l,k)*psi(inode,k);
     }
    }

   // Copy over (the equals operator does a shallow copy ?)
   for(unsigned l=0;l<6;++l)
    { psi(inode,l)=psi_new[l]; }
  }

 // (Shallow?) copy the basis functions to test functions
 test = psi;

 // Leave the unused bubble shape functions empty
}

//============================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//=============================================================================
template<unsigned DIM, unsigned NNODE_1D>
 double FoepplVonKarmanBellElement<DIM,NNODE_1D>::
  dshape_u_and_dtest_u_eulerian_foeppl_von_karman(const Vector<double> &s,Shape &psi,
  DShape &dpsidx,  Shape &test,  DShape &dtestdx) const
{
 // Initialise
 const unsigned n_node = this->nnode();
 DShape dpsids(n_node,DIM);
 dshape_u_local(s,psi,dpsids);

 Shape shape(3);
 DShape dshape_ds(3,DIM);
 this->dLshape_local(s,shape, dshape_ds);
 // Get the inverse jacobian
 DenseMatrix<double> jacobian(DIM,DIM,0.0),inverse_jacobian(DIM,DIM,0.0);
 const double J = this->local_to_eulerian_mapping2(dshape_ds,jacobian,
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
       // Convert to local coordinates
       dpsidx(l,i)+=dpsids(l,j)*inverse_jacobian(i,j);
      }
    }
 }

 test = psi;
 dtestdx = dpsidx;

 return J;
}
//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned DIM, unsigned NNODE_1D>
 double FoepplVonKarmanBellElement<DIM,NNODE_1D>::
 dshape_and_dtest_eulerian_foeppl_von_karman(const Vector<double> &s, Shape &psi,
 Shape& psi_b, DShape &dpsidx, DShape& dpsi_b_dx,  Shape &test, Shape& test_b,
 DShape &dtestdx,DShape &dtest_b_dx) const
{
 //const double J = this->dbasis_eulerian(s,psi,dpsidx);
 DShape d2psidx(3,6,3);
 const double J =this->J_eulerian1(s);

 // Vertices
 Vector<Vector<double> > v(3,Vector<double>(2));
 for (unsigned inode=0;inode<3;++inode)
  {
   // Get the position vector
   Node* nod_pt=this->node_pt(inode);
   v[inode][0]=nod_pt->x(0);
   v[inode][1]=nod_pt->x(1);
  }

 MyShape::d2_basis_eulerian(s,v,psi,dpsidx,d2psidx);
 // Loop over nodes with rotated dofs
 for(unsigned i=0; i<Nnodes_to_rotate; ++i)
  {
   // Get the nodes
   unsigned inode = Nodes_to_rotate[i];
   // Construct the new psi
   Vector<double> psi_new(6,0.0);
   DenseMatrix<double> dpsi_new(6,2,0.0), rotation_matrix(6,6,0.0);
   this->rotation_matrix_at_node(inode,rotation_matrix);

   // Now loop over the shape functions at node i
   for(unsigned l=0;l<6;++l)
    {
     for(unsigned k=0;k<6;++k)
     {
      // Now rotate the shape functions
      psi_new[l]+=rotation_matrix(l,k)*psi(inode,k);
      dpsi_new(l,0)+=rotation_matrix(l,k)*dpsidx(inode,k,0);
      dpsi_new(l,1)+=rotation_matrix(l,k)*dpsidx(inode,k,1);
     }
    }

   // Copy over (the equals operator does a shallow copy ?)
   for(unsigned l=0;l<6;++l)
    {
     psi(inode,l)=psi_new[l];
     dpsidx(inode,l,0)=dpsi_new(l,0);
     dpsidx(inode,l,1)=dpsi_new(l,1);
    }
  }

 test = psi;
 dtestdx = dpsidx;
 // Leave the unused bubble shape functions empty

 return J;
}

template<unsigned DIM, unsigned NNODE_1D>
 double FoepplVonKarmanBellElement<DIM,NNODE_1D>::
  d2shape_and_d2test_eulerian_foeppl_von_karman(const Vector<double> &s,  Shape &psi,
  Shape &psi_b, DShape &dpsidx, DShape &dpsi_bdx,  DShape &d2psidx,
  DShape &d2psi_bdx,
  Shape &test, Shape &test_b, DShape &dtestdx, DShape &dtest_bdx,
  DShape &d2testdx, DShape &d2test_bdx) const
{
 //Call the geometrical shape functions and derivatives
 const double J =this->J_eulerian1(s);

 // Vertices
 Vector<Vector<double> > v(3,Vector<double>(2));
 for (unsigned inode=0;inode<3;++inode)
  {
   // Get the position vector
   Node* nod_pt=this->node_pt(inode);
   v[inode][0]=nod_pt->x(0);
   v[inode][1]=nod_pt->x(1);
  }

 MyShape::d2_basis_eulerian(s,v,psi,dpsidx,d2psidx);
 // Loop over nodes with rotated dofs
 // Hierher?
 for(unsigned i=0; i<Nnodes_to_rotate; ++i)
  {
   // Get the nodes
   unsigned inode = Nodes_to_rotate[i];
   // Construct the new psi
   Vector<double> psi_new(6,0.0);
   DenseMatrix<double> dpsi_new(6,2,0.0), d2psi_new(6,3,0.0),
     rotation_matrix(6,6,0.0);
   this->rotation_matrix_at_node(inode,rotation_matrix);

   // Now loop over the shape functions at node i
   for(unsigned l=0;l<6;++l)
    {
     for(unsigned k=0;k<6;++k)
     {
      // Now rotate the shape functions
      psi_new[l]+=rotation_matrix(l,k)*psi(inode,k);
      dpsi_new(l,0)+=rotation_matrix(l,k)*dpsidx(inode,k,0);
      dpsi_new(l,1)+=rotation_matrix(l,k)*dpsidx(inode,k,1);
      d2psi_new(l,0)+=rotation_matrix(l,k)*d2psidx(inode,k,0);
      d2psi_new(l,1)+=rotation_matrix(l,k)*d2psidx(inode,k,1);
      d2psi_new(l,2)+=rotation_matrix(l,k)*d2psidx(inode,k,2);
     }
    }

   // Copy over (the equals operator does a shallow copy ?)
/*   std::cout<<"w :"<<
      psi(inode,0)-psi_new[0]<<
      psi(inode,1)-psi_new[2]<<
      psi(inode,2)-psi_new[1]<<
      psi(inode,3)-psi_new[4]<<
      psi(inode,4)-psi_new[3]<<
      psi(inode,5)-psi_new[5]<<"\n";
   for(unsigned m=0;m<2;++m)
   std::cout<<"Dw: "<<
      dpsidx(inode,0,m)-dpsi_new(0,m)<<
      dpsidx(inode,1,m)-dpsi_new(2,m)<<
      dpsidx(inode,2,m)-dpsi_new(1,m)<<
      dpsidx(inode,3,m)-dpsi_new(4,m)<<
      dpsidx(inode,4,m)-dpsi_new(3,m)<<
      dpsidx(inode,5,m)-dpsi_new(5,m)<<"\n";
   for(unsigned m=0;m<3;++m)
   std::cout<<"D2w: "<<
      d2psidx(inode,0,m)-d2psi_new(0,m)<<
      d2psidx(inode,1,m)-d2psi_new(2,m)<<
      d2psidx(inode,2,m)-d2psi_new(1,m)<<
      d2psidx(inode,3,m)-d2psi_new(4,m)<<
      d2psidx(inode,4,m)-d2psi_new(3,m)<<
      d2psidx(inode,5,m)-d2psi_new(5,m)<<"\n";*/
   for(unsigned l=0;l<6;++l)
    {
     psi(inode,l)=psi_new[l];
     dpsidx(inode,l,0)=dpsi_new(l,0);
     dpsidx(inode,l,1)=dpsi_new(l,1);
     d2psidx(inode,l,0)=d2psi_new(l,0);
     d2psidx(inode,l,1)=d2psi_new(l,1);
     d2psidx(inode,l,2)=d2psi_new(l,2);
    }
  }

 //Set the test functions equal to the shape functions
 test = psi;
 dtestdx= dpsidx;
 d2testdx = d2psidx;

 //Return the jacobian
 return J;
}

template<unsigned DIM, unsigned NNODE_1D>
 void FoepplVonKarmanBellElement<DIM,NNODE_1D>::pin_all_deflection_dofs() const
 {
  // CHECK HERE
  // Bell elements only have deflection dofs at vertices
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
 }

template<unsigned DIM, unsigned NNODE_1D>
 void FoepplVonKarmanBellElement<DIM,NNODE_1D>::
 fix_out_of_plane_displacement_dof(const unsigned& dof_number, const unsigned&
boundary_number, const DisplacementFctPt& w)
 {
  // CHECK HERE
  // Bell elements only have deflection dofs at vertices
  for(unsigned n=0; n<3; ++n)
   {
    // Get node
    Node* nod_pt=this->node_pt(n);
    // Check if it is on the boundary
    bool is_boundary_node=nod_pt->is_on_boundary(boundary_number);
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

template<unsigned DIM, unsigned NNODE_1D>
 void FoepplVonKarmanBellElement<DIM,NNODE_1D>::
 fix_in_plane_displacement_dof(const unsigned& dof_number, const unsigned&
boundary_number, const DisplacementFctPt& u)
 {
  // CHECK HERE
   const unsigned n_node = this->nnode();
  // Bell elements only have deflection dofs at vertices
  for(unsigned n=0; n<n_node; ++n)
   {
    // Get node
    Node* nod_pt=this->node_pt(n);
    // Check if it is on the boundary
    bool is_boundary_node=nod_pt->is_on_boundary(boundary_number);
    if(is_boundary_node)
     {
      // Extract nodal coordinates from node:
      Vector<double> x(2);
      x[0]=nod_pt->x(0);
      x[1]=nod_pt->x(1);
      // Get value
      double value;
      u(x,value);
      // Pin and set the value
      nod_pt->pin(dof_number);
      nod_pt->set_value(dof_number,value);
     }
   }

 }
////======================================================================
///// Define the shape functions and test functions and derivatives
///// w.r.t. global coordinates and return Jacobian of mapping.
/////
///// Galerkin: Test functions = shape functions
////======================================================================
//template<unsigned DIM, unsigned NNODE_1D>
//double FoepplVonKarmanBellElement<DIM,NNODE_1D>::
// dshape_and_dtest_eulerian_at_knot_foeppl_von_karman(
//  const unsigned &ipt,
//  Shape &psi,
//  DShape &dpsidx,
//  Shape &test,
//  DShape &dtestdx) const
//{
// const double J = this->dshape_and_dtest_eulerian_at_knot(ipt,psi,dpsidx,test,dtestdx);
//
// return J;
//}


//template<unsigned DIM, unsigned NNODE_1D>
//double FoepplVonKarmanBellElement<DIM,NNODE_1D>::
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
// const double J = this->d2shape_and_d2test_eulerian_at_knot(ipt,psi,dpsidx,d2psidx,test,dtestdx,d2testdx);
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
template <unsigned DIM, unsigned NNODE_1D>
  void FoepplVonKarmanBellElement<DIM,NNODE_1D>::Lshape(const Vector<double> &s, Shape &psi) const
   {
    psi[0] = s[0];
    psi[1] = s[1];
    psi[2] = 1.0-s[0]-s[1];
   }
  

//=======================================================================
/// Derivatives of shape functions for specific TElement<2,2>
//=======================================================================
template <unsigned DIM, unsigned NNODE_1D>
  void FoepplVonKarmanBellElement<DIM,NNODE_1D>::dLshape_local(const Vector<double> &s,
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
  


template <unsigned DIM, unsigned NNODE_1D>
double FoepplVonKarmanBellElement<DIM,NNODE_1D>::local_to_eulerian_mapping2(const DShape &dpsids,
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
template <unsigned DIM, unsigned NNODE_1D>
 double FoepplVonKarmanBellElement<DIM,NNODE_1D>::J_eulerian1(const Vector<double> &s) const
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
template <unsigned DIM, unsigned NNODE_1D>
  void FoepplVonKarmanBellElement<DIM,NNODE_1D>::my_interpolated_x(const Vector<double> &s, Vector<double> &x)
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
