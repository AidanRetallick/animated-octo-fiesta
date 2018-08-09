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
// Non--inline functions for the foeppl_von_karman equations
#include "foeppl_von_karman.h"

namespace oomph
{
 /// \short Default to incompressible sheet Poisson ratio
 template <unsigned DIM, unsigned NNODE_1D>
 const double FoepplVonKarmanEquations<DIM,NNODE_1D>::Default_Nu_Value=0.5;
 
 /// \short Default to 'natural' FvK nondimensionalisation which is h 
 /// independent (NB NOT the same as in JFM paper)
 template <unsigned DIM, unsigned NNODE_1D>
 const double FoepplVonKarmanEquations<DIM,NNODE_1D>::Default_Eta_Value=1;


//======================================================================
/// Fill in generic residual and jacobian contribution 
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
void  FoepplVonKarmanEquations<DIM,NNODE_1D>::
fill_in_generic_residual_contribution_foeppl_von_karman(Vector<double> &residuals,
                                              DenseMatrix<double> &jacobian,
                                              const unsigned& flag)
{
 // CALL GET SHAPE ASSOCIATION MATRIX HERE
 //Find out how many nodes there are
 const unsigned n_node = this->nnode();
 const unsigned n_node_w = 3; // HERE
 //Find out how many bubble nodes there are
 const unsigned n_b_node = this->Number_of_internal_dofs;
 //Find out how many nodes positional dofs there are
 unsigned n_position_type = this->nnodal_position_type();
 // Find the internal dofs
 const unsigned n_b_position_type = this->Number_of_internal_dof_types;

 //Get the Poisson ratio of the plate
 const double nu=get_nu();

 //Set up memory for the shape and test functions
 Shape psi_u(n_node), test_u(n_node);
 DShape dpsi_udxi(n_node,DIM), dtest_udxi(n_node,DIM);

 //Local c1-shape funtion
 Shape psi(n_node_w,n_position_type),test(n_node_w,n_position_type),
  psi_b(n_b_node,n_b_position_type),test_b(n_b_node,n_b_position_type);
 
 DShape dpsi_dxi(n_node_w,n_position_type,DIM),dtest_dxi(n_node_w,n_position_type,DIM),
  dpsi_b_dxi(n_b_node,n_b_position_type,DIM),dtest_b_dxi(n_b_node,n_b_position_type,DIM),
  d2psi_dxi2(n_node_w,n_position_type,3), d2test_dxi2(n_node_w,n_position_type,3),
  d2psi_b_dxi2(n_b_node,n_b_position_type,3), d2test_b_dxi2(n_b_node,n_b_position_type,3);

 //Set the value of n_intpt
 const unsigned n_intpt = this->integral_pt()->nweight();
 //Integers to store the local equation and unknown numbers
 int local_eqn=0, local_unknown=0;

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the integral weight
   double w = this->integral_pt()->weight(ipt);

   //Calculate values of unknown
   Vector<double> interpolated_w(1,0.0);
   DenseMatrix<double> interpolated_dwdxi(1,DIM,0.0);
   DenseMatrix<double> interpolated_d2wdxi2(1,3,0.0);
   
   //Calculate values of unknown
   Vector<double> interpolated_u(2,0.0);
   DenseMatrix<double> interpolated_duidxj(2,DIM,0.0);

   //Allocate and initialise to zero
   Vector<double> interpolated_x(DIM,0.0);
   Vector<double> s(DIM);
   s[0] = this->integral_pt()->knot(ipt,0);
   s[1] = this->integral_pt()->knot(ipt,1);
   this->get_coordinate_x(s,interpolated_x);
   // CALL MODIFIED SHAPE (THAT TAKES ASSOCIATION MATRIX AS AN ARGUMENT) HERE
   // MAKE SURE THAT THE MULTIPLICATION IS EFFICIENT USING BLOCK STRUCTURE
   //Call the derivatives of the shape and test functions for the unknown
   double J = d2shape_and_d2test_eulerian_foeppl_von_karman(s,
    psi, psi_b, dpsi_dxi, dpsi_b_dxi, d2psi_dxi2, d2psi_b_dxi2,
    test, test_b, dtest_dxi, dtest_b_dxi, d2test_dxi2, d2test_b_dxi2);

   dshape_u_and_dtest_u_eulerian_foeppl_von_karman(s,psi_u,dpsi_udxi,test_u,dtest_udxi);
   //Premultiply the weights and the Jacobian
   double W = w*J;

   //Calculate function value and derivatives
   //-----------------------------------------
   // Loop over nodes
   for(unsigned l=0;l<n_node_w;l++)
    {
     for(unsigned k=0;k<n_position_type;k++)
      {
       //Get the nodal value of the unknown
       double w_value = this->raw_nodal_value(l,k+2);
       interpolated_w[0] += w_value*psi(l,k);
       // Loop over directions
       for(unsigned j=0;j<DIM;j++)
        {
         interpolated_dwdxi(0,j) += w_value*dpsi_dxi(l,k,j);
        }
       for(unsigned j=0;j<3;j++)
        {
          interpolated_d2wdxi2(0,j) += w_value*d2psi_dxi2(l,k,j);
        }
      }
    }

   // Loop over internal dofs
   for(unsigned l=0;l<Number_of_internal_dofs;l++)
   {
    for(unsigned k=0;k<n_b_position_type;k++)
     {
      //Get the nodal value of the unknown
      double u_value = get_w_bubble_dof(l,k);
      interpolated_w[0] += u_value*psi_b(l,k);
      // Loop over directions
      for(unsigned j=0;j<DIM;j++)
       {
        interpolated_dwdxi(0,j) += u_value*dpsi_b_dxi(l,k,j);
       }
      for(unsigned j=0;j<3;j++)
       {
        interpolated_d2wdxi2(0,j) += u_value*d2psi_b_dxi2(l,k,j);
       }
     }
    }

   // Loop over nodes
   for(unsigned l=0;l<n_node;l++)
    {
     for(unsigned alpha=0;alpha<2;alpha++)
      {
       double u_value = this->raw_nodal_value(l,alpha);
       interpolated_u[alpha] += u_value*psi_u(l);
  
        // Loop over directions
        for(unsigned beta=0; beta<DIM; beta++)
         {
          interpolated_duidxj(alpha,beta) += u_value*dpsi_udxi(l,beta);
         }
      }
     }

   // Get the stress
   DenseMatrix<double> sigma(2,2,0.0);
   get_sigma(sigma,interpolated_duidxj, interpolated_dwdxi); 
   
    // Check 
//    Vector<double> interpolated_w(6,0.0)pli:
//    interpolated_w=this-> interpolated_w_foeppl_von_karman(s);
//   interpolated_w[0]=interpolated_w[0];
//   interpolated_dwdxi(0,0)=interpolated_w[1];
//   interpolated_dwdxi(0,1)=interpolated_w[2];
//   interpolated_d2wdxi2(0,0)=interpolated_w[3];
//   interpolated_d2wdxi2(0,1)=interpolated_w[4];
//   interpolated_d2wdxi2(0,2)=interpolated_w[5];
   //Get pressure function
   //-------------------
   double pressure(0.0);
   Vector<double> pressure_gradient(2,0.0);
   get_pressure_foeppl_von_karman(ipt,interpolated_x,pressure);
   get_in_plane_forcing_foeppl_von_karman(ipt,interpolated_x,pressure_gradient);

   // Loop over the nodal test functions
   for(unsigned l=0;l<n_node_w;l++)
    {
    for(unsigned k=0;k<n_position_type;k++)
     {
     //Get the local equation
     local_eqn = this->nodal_local_eqn(l,k+2);
     //IF it's not a boundary condition
     if(local_eqn >= 0)
      {
      // Add body force/pressure term here
      residuals[local_eqn] -= pressure*test(l,k)*W;
      for(unsigned alpha=0;alpha<2;++alpha)
       {
       for(unsigned beta=0; beta<2;++beta)
        {
        // w_{,\alpha\beta} \delta \kappa_\alpha\beta
        residuals[local_eqn] += (1-nu)/(12.*(1-nu*nu))
         *interpolated_d2wdxi2(0,alpha+beta)*d2test_dxi2(l,k,alpha+beta)*W;
        // w_{,\alpha\alpha} \delta \kappa_\beta\beta
        residuals[local_eqn] += (nu)/(12.*(1-nu*nu))
         *interpolated_d2wdxi2(0,beta+beta)*d2test_dxi2(l,k,alpha+alpha)*W;
        // sigma_{\alpha\beta} w_\alpha \delta w_{,\beta}
        residuals[local_eqn] += eta()*sigma(alpha,beta)
         *interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)*W;
        }
       }
      // Calculate the jacobian
      //-----------------------
      if(flag)
      {
       //Loop over the in--plane unknowns again
       for(unsigned ll=0;ll<n_node;ll++)
        {
        // Loop over displacements
        for(unsigned gamma=0;gamma<2;++gamma)
         {
         local_unknown = this->nodal_local_eqn(ll,gamma);
         // If at a non-zero degree of freedom add in the entry
         if(local_unknown >= 0)
          {
          // Loop over displacements
          for(unsigned alpha=0;alpha<2;++alpha)
           {
           // Loop over dimensions 
           for(unsigned beta=0; beta<2;++beta)
            {
            // The Nonlinear Terms 
            if(alpha==beta)
             {
             jacobian(local_eqn,local_unknown) += eta()*nu*dpsi_udxi(ll,gamma)*
                   interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)/(1-nu*nu)*W;
             }
            if(alpha==gamma)
             {
             jacobian(local_eqn,local_unknown) += eta()*(1-nu)/2.*dpsi_udxi(ll,beta)
                   *interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)/(1-nu*nu)*W;
             }
            if(beta==gamma)
             {
             jacobian(local_eqn,local_unknown) += eta()*(1-nu)/2.*dpsi_udxi(ll,alpha)
                   *interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)/(1-nu*nu)*W;
             }
            }
           }
          }
         }
        }

       //Loop over the test functions again
       for(unsigned l2=0;l2<n_node_w;l2++)
        {
        // Loop over position dofs
        for(unsigned k2=0;k2<n_position_type;k2++)
         {
          local_unknown = this->nodal_local_eqn(l2,k2+2);
          // If at a non-zero degree of freedom add in the entry
          if(local_unknown >= 0)
           {
           // Loop over dimensions 
           for(unsigned alpha=0;alpha<2;++alpha)
            {
            for(unsigned beta=0; beta<2;++beta)
             {
             // Linear Terms
             // w_{,\alpha\beta} \delta \kappa_\alpha\beta
             jacobian(local_eqn,local_unknown) += (1-nu)/(12.*(1-nu*nu))
                *d2psi_dxi2(l2,k2,alpha+beta)*d2test_dxi2(l,k,alpha+beta)*W;
             // w_{,\alpha\alpha} \delta \kappa_\beta\beta
             jacobian(local_eqn,local_unknown) += nu/(12.*(1-nu*nu))
                *d2psi_dxi2(l2,k2,beta+beta)*d2test_dxi2(l,k,alpha+alpha)*W;
             // Nonlinear terms
             jacobian(local_eqn,local_unknown) += eta()*sigma(alpha,beta) 
                        *dpsi_dxi(l2,k2,alpha)*dtest_dxi(l,k,beta)*W;
             // Nonlinear terms
             jacobian(local_eqn,local_unknown) += eta()*nu/(1-nu*nu)*
               interpolated_dwdxi(0,alpha)*dpsi_dxi(l2,k2,alpha)
               *interpolated_dwdxi(0,beta)*dtest_dxi(l,k,beta)*W;
             jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
               interpolated_dwdxi(0,alpha)*dpsi_dxi(l2,k2,beta)
               *interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)*W;
             jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
               interpolated_dwdxi(0,beta)*dpsi_dxi(l2,k2,alpha)
               *interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)*W;
             }
            }
           }
         }
        }
       //Loop over the internal test functions
       for(unsigned l2=0;l2<Number_of_internal_dofs;l2++)
        {
        for(unsigned k2=0;k2<n_b_position_type;k2++)
         {
         local_unknown = local_w_bubble_equation(l2,k2);
         //If at a non-zero degree of freedom add in the entry
         if(local_unknown >= 0)
          {
          // Loop over dimensions 
          for(unsigned alpha=0;alpha<2;++alpha)
           {
           for(unsigned beta=0; beta<2;++beta)
            {
            // The Linear Terms 
            // w_{,\alpha\beta} \delta \kappa_\alpha\beta
            jacobian(local_eqn,local_unknown) += (1-nu)/(12.*(1-nu*nu))
               *d2psi_b_dxi2(l2,k2,alpha+beta)*d2test_dxi2(l,k,alpha+beta)*W;
            // w_{,\alpha\alpha} \delta \kappa_\beta\beta
            jacobian(local_eqn,local_unknown) += nu/(12.*(1-nu*nu))
               *d2psi_b_dxi2(l2,k2,beta+beta)*d2test_dxi2(l,k,alpha+alpha)*W;
            // Nonlinear terms
            jacobian(local_eqn,local_unknown) += eta()*sigma(alpha,beta) 
                       *dpsi_b_dxi(l2,k2,alpha)*dtest_dxi(l,k,beta)*W;
            // Nonlinear terms
            jacobian(local_eqn,local_unknown) += eta()*nu/(1-nu*nu)*
              interpolated_dwdxi(0,alpha)*dpsi_b_dxi(l2,k2,alpha)
              *interpolated_dwdxi(0,beta)*dtest_dxi(l,k,beta)*W;
            jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
              interpolated_dwdxi(0,alpha)*dpsi_b_dxi(l2,k2,beta)
              *interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)*W;
            jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
              interpolated_dwdxi(0,beta)*dpsi_b_dxi(l2,k2,alpha)
              *interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)*W;
            }
           }
          }
         }
        }
       } // End of flag
      }
     }
    }

    // Loop over the internal test functions
    for(unsigned l=0;l<Number_of_internal_dofs;l++)
    {
    for(unsigned k=0;k<n_b_position_type;k++)
     {
     //Get the local equation
     local_eqn = local_w_bubble_equation(l,k);
     //IF it's not a boundary condition
     if(local_eqn >= 0)
      {
      // Add body force/pressure term here
      residuals[local_eqn] -= pressure*test_b(l,k)*W;

      for(unsigned alpha=0;alpha<2;++alpha)
       {
       for(unsigned beta=0; beta<2;++beta)
        {
        // The Linear Terms 
        // w_{,\alpha\beta} \delta \kappa_\alpha\beta
        residuals[local_eqn] += (1-nu)/(12.*(1-nu*nu))
          *interpolated_d2wdxi2(0,alpha+beta)*d2test_b_dxi2(l,k,alpha+beta)*W;
        // w_{,\alpha\alpha} \delta \kappa_\beta\beta
        residuals[local_eqn] += (nu)/(12.*(1-nu*nu))
          *interpolated_d2wdxi2(0,beta+beta)*d2test_b_dxi2(l,k,alpha+alpha)*W;
        // sigma_{\alpha\beta} w_\alpha \delta w_{,\beta}
        residuals[local_eqn] += eta()*sigma(alpha,beta)
          *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)*W;
        }
       }
      // Calculate the jacobian
      //-----------------------
      if(flag)
       {
       //Loop over the in--plane unknowns again
       for(unsigned ll=0;ll<n_node;ll++)
        {
        // Loop over displacements
        for(unsigned gamma=0;gamma<2;++gamma)
         {
         local_unknown = this->nodal_local_eqn(ll,gamma);
         // If at a non-zero degree of freedom add in the entry
         if(local_unknown >= 0)
          {
          // Loop over displacements
          for(unsigned alpha=0;alpha<2;++alpha)
           {
           // Loop over dimensions 
           for(unsigned beta=0; beta<2;++beta)
            {
            // The Nonlinear Terms 
            if(alpha==beta)
             {
             jacobian(local_eqn,local_unknown) += eta()*nu*dpsi_udxi(ll,gamma)*
                   interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)/(1-nu*nu)*W;
             }
            if(alpha==gamma)
             {
             jacobian(local_eqn,local_unknown) += eta()*(1-nu)/2.*dpsi_udxi(ll,beta)
                   *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)/(1-nu*nu)*W;
             }
            if(beta==gamma)
             {
             jacobian(local_eqn,local_unknown) += eta()*(1-nu)/2.*dpsi_udxi(ll,alpha)
                   *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)/(1-nu*nu)*W;
             }
            }
           }
          }
         }

        }
       //Loop over the test functions again
       for(unsigned l2=0;l2<n_node_w;l2++)
        {
         // Loop over position dofs
         for(unsigned k2=0;k2<n_position_type;k2++)
          {
          local_unknown = this->nodal_local_eqn(l2,k2+2);
          //If at a non-zero degree of freedom add in the entry
          if(local_unknown >= 0)
           {
           // Loop over dimensions 
           for(unsigned alpha=0;alpha<2;++alpha)
            {
            for(unsigned beta=0; beta<2;++beta)
             {
             // The Linear Terms 
             // w_{,\alpha\beta} \delta \kappa_\alpha\beta
             jacobian(local_eqn,local_unknown) += (1-nu)/(12.*(1-nu*nu))
                *d2psi_dxi2(l2,k2,alpha+beta)*d2test_b_dxi2(l,k,alpha+beta)*W;
             // w_{,\alpha\alpha} \delta \kappa_\beta\beta
             jacobian(local_eqn,local_unknown) += nu/(12.*(1-nu*nu))
                *d2psi_dxi2(l2,k2,beta+beta)*d2test_b_dxi2(l,k,alpha+alpha)*W;
             // Nonlinear terms
             jacobian(local_eqn,local_unknown) += eta()*sigma(alpha,beta) 
                        *dpsi_dxi(l2,k2,alpha)*dtest_b_dxi(l,k,beta)*W;
             // Nonlinear terms
             jacobian(local_eqn,local_unknown) += eta()*nu/(1-nu*nu)*
               interpolated_dwdxi(0,alpha)*dpsi_dxi(l2,k2,alpha)
               *interpolated_dwdxi(0,beta)*dtest_b_dxi(l,k,beta)*W;
             jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
               interpolated_dwdxi(0,alpha)*dpsi_dxi(l2,k2,beta)
               *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)*W;
             jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
               interpolated_dwdxi(0,beta)*dpsi_dxi(l2,k2,alpha)
               *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)*W;
             }
            }
           }
          }
         }
       //Loop over the test functions again
       for(unsigned l2=0;l2<Number_of_internal_dofs;l2++)
        {
        // Loop over position dofs
        for(unsigned k2=0;k2<n_b_position_type;k2++)
         {
         local_unknown = local_w_bubble_equation(l2,k2);
         //If at a non-zero degree of freedom add in the entry
         if(local_unknown >= 0)
          {
          // Loop over dimensions 
          for(unsigned alpha=0;alpha<2;++alpha)
           {
           for(unsigned beta=0; beta<2;++beta)
            {
            // The Linear Terms 
            // w_{,\alpha\beta} \delta \kappa_\alpha\beta
            jacobian(local_eqn,local_unknown) += (1-nu)/(12.*(1-nu*nu))
               *d2psi_b_dxi2(l2,k2,alpha+beta)*d2test_b_dxi2(l,k,alpha+beta)*W;
            // w_{,\alpha\alpha} \delta \kappa_\beta\beta
            jacobian(local_eqn,local_unknown) += nu/(12.*(1-nu*nu))
               *d2psi_b_dxi2(l2,k2,beta+beta)*d2test_b_dxi2(l,k,alpha+alpha)*W;
            // Nonlinear terms
            jacobian(local_eqn,local_unknown) += eta()*sigma(alpha,beta) 
                       *dpsi_b_dxi(l2,k2,alpha)*dtest_b_dxi(l,k,beta)*W;
            // Nonlinear terms
            jacobian(local_eqn,local_unknown) += eta()*nu/(1-nu*nu)*
              interpolated_dwdxi(0,alpha)*dpsi_b_dxi(l2,k2,alpha)
              *interpolated_dwdxi(0,beta)*dtest_b_dxi(l,k,beta)*W;
            jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
              interpolated_dwdxi(0,alpha)*dpsi_b_dxi(l2,k2,beta)
              *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)*W;
            jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
              interpolated_dwdxi(0,beta)*dpsi_b_dxi(l2,k2,alpha)
              *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)*W;
            }
           }
          }
         }
        }
       } // End of flag
      }
     }
    } //end loop over nodes
  for(unsigned l=0;l<n_node;++l)
   {
    // Now loop over displacement equations
    for(unsigned alpha=0; alpha<2 ; ++alpha)
     {
     //Get the local equation
     local_eqn = this->nodal_local_eqn(l,alpha);
     //IF it's not a boundary condition
     if(local_eqn >= 0)
      {
      // Now add on the stress terms
      for(unsigned beta=0; beta<2;++beta)
        {
         // Variations in u_alpha
         residuals[local_eqn] += sigma(alpha,beta)*dtest_udxi(l,beta)*W;
        }
      // Add the forcing term
      residuals[local_eqn] += pressure_gradient[alpha]*psi_u(l)*W;

      // Now loop over Jacobian
      if(flag)
       {
       for(unsigned ll=0;ll<n_node;++ll)
        {
        // Now loop over displacement equations
        for(unsigned gamma=0; gamma<2 ; ++gamma)
         {
         //Get the local equation
         local_unknown = this->nodal_local_eqn(ll,gamma);
         // If the unknown is a degree of freedom 
         if(local_unknown>=0)
          {
          // Now add on the stress terms
          for(unsigned beta=0; beta<2;++beta)
           {
           // The Linear Terms 
           if(alpha==beta)
            {
            jacobian(local_eqn,local_unknown) += nu*dpsi_udxi(ll,gamma)
                  *dtest_udxi(l,beta)/(1-nu*nu)*W;
            }
           if(alpha==gamma)
            {
            jacobian(local_eqn,local_unknown) += (1-nu)/2.*dpsi_udxi(ll,beta)
                  *dtest_udxi(l,beta)/(1-nu*nu)*W;
            }
           if(beta==gamma)
            {
            jacobian(local_eqn,local_unknown) += (1-nu)/2.*dpsi_udxi(ll,alpha)
                  *dtest_udxi(l,beta)/(1-nu*nu)*W;
            }
           }// End loop beta
          }// End if local unknown
         }// End loop gamma dof
        }// End loop nodal dofs
       
       //Loop over the w nodal dofs
       for(unsigned l2=0;l2<n_node_w;l2++)
        {
        // Loop over position dofs
        for(unsigned k2=0;k2<n_position_type;k2++)
         {
         //If at a non-zero degree of freedom add in the entry
         //Get the local equation
         local_unknown = this->nodal_local_eqn(l2,k2+2);
         // If the unknown is a degree of freedom 
         if(local_unknown>=0)
          {
          // Now add on the stress terms
          for(unsigned beta=0; beta<2;++beta)
           {
           // The NonLinear Terms 
           jacobian(local_eqn,local_unknown) += eta()*nu*dpsi_dxi(l2,k2,beta)
                 *interpolated_dwdxi(0,beta)*dtest_udxi(l,alpha)/(1-nu*nu)*W;
           jacobian(local_eqn,local_unknown) +=eta()*(1-nu)/2.*dpsi_dxi(l2,k2,alpha)
                 *interpolated_dwdxi(0,beta)*dtest_udxi(l,beta)/(1-nu*nu)*W;
           jacobian(local_eqn,local_unknown) +=eta()*(1-nu)/2.*dpsi_dxi(l2,k2,beta)
                 *interpolated_dwdxi(0,alpha)*dtest_udxi(l,beta)/(1-nu*nu)*W;
           }// End loop beta
          }// End if local unknown
         }// End loop position type dof
        }// End loop w nodal dofs
       
       //Loop over the w internal dofs
       for(unsigned l2=0;l2<Number_of_internal_dofs;l2++)
        {
        // Loop over position dofs
        for(unsigned k2=0;k2<n_b_position_type;k2++)
         {
         //If at a non-zero degree of freedom add in the entry
         //Get the local equation
         local_unknown = local_w_bubble_equation(l2,k2);
         //If at a non-zero degree of freedom add in the entry
         if(local_unknown >= 0)
          {
          // Now add on the stress terms
          for(unsigned beta=0; beta<2;++beta)
           {
           // The NonLinear Terms 
           jacobian(local_eqn,local_unknown) += eta()*nu*dpsi_b_dxi(l2,k2,beta)
                 *interpolated_dwdxi(0,beta)*dtest_udxi(l,alpha)/(1-nu*nu)*W;
           jacobian(local_eqn,local_unknown) +=eta()*(1-nu)/2.*dpsi_b_dxi(l2,k2,alpha)
                 *interpolated_dwdxi(0,beta)*dtest_udxi(l,beta)/(1-nu*nu)*W;
           jacobian(local_eqn,local_unknown) +=eta()*(1-nu)/2.*dpsi_b_dxi(l2,k2,beta)
                 *interpolated_dwdxi(0,alpha)*dtest_udxi(l,beta)/(1-nu*nu)*W;
           } // End loop beta
          } // End if local_unknown
         } // End loop position type
        } // End loop internal dofs
       } // End if flag
      } // End if local eqn
     } // End loop over alpha
    } // End loop over nodes
  } // End of loop over integration points
}



//======================================================================
/// Self-test:  Return 0 for OK
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
unsigned  FoepplVonKarmanEquations<DIM,NNODE_1D>::self_test()
{

 bool passed=true;

 // Check lower-level stuff
 if (FiniteElement::self_test()!=0)
  {
   passed=false;
  }

 // Return verdict
 if (passed)
  {
   return 0;
  }
 else
  {
   return 1;
  }

}

//======================================================================
/// Output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
void  FoepplVonKarmanEquations<DIM,NNODE_1D>::output(std::ostream &outfile,
                                    const unsigned &nplot)
{

 //Vector of local coordinates
 Vector<double> s(DIM),x(DIM);

 // Tecplot header info
 outfile << this->tecplot_zone_string(nplot);

 // Loop over plot points
// Vector<double> u(this->required_nvalue(0),0.0);
 Vector<double> u(15,0.0);
 unsigned num_plot_points=this->nplot_points(nplot);
 Vector<double> r(3);

 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   // Get local coordinates of plot point
   this->get_s_plot(iplot,nplot,s);
   u = interpolated_u_foeppl_von_karman(s,true);

   // Get x position as Vector
   this->get_coordinate_x(s,x);

//  // // Get the interpolated x
//  unsigned n_node=this->nnode();
//  Vector<double> interpolated_x(2,0.0);
//  Shape psi_u(n_node), test_u(n_node);
//  // DShape dpsi_udxi(n_node,DIM), dtest_udxi(n_node,DIM);
//  // dshape_u_and_dtest_u_eulerian_foeppl_von_karman(s,psi_u,dpsi_udxi,test_u,dtest_udxi);
//
//   // Check the u geometry
//   // Loop over nodes
//   for(unsigned l=0;l<n_node;l++)
//    {
//     // Loop over directions
//     for(unsigned j=0;j<2;j++)
//      {
//       interpolated_x[j] += this->raw_nodal_position(l,j)*psi_u[l];
//      }
//    }
//
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << x[i] << " ";
    }

   // Loop for variables
  // for(unsigned j=0;j<this->required_nvalue(0);j++)
   for(unsigned j=0;j<15;j++)
    {
     outfile << u[j] << " " ;
    }

   outfile << std::endl;
  }

 // Write tecplot footer (e.g. FE connectivity lists)
 this->write_tecplot_zone_footer(outfile,nplot);
}


//======================================================================
/// C-style output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
void  FoepplVonKarmanEquations<DIM,NNODE_1D>::output(FILE* file_pt,
                                    const unsigned &nplot)
{
 //Vector of local coordinates
 Vector<double> s(DIM), x(DIM);;

 // Tecplot header info
 fprintf(file_pt,"%s",this->tecplot_zone_string(nplot).c_str());

 // Loop over plot points
 Vector<double> u(this->required_nvalue(0),0.0);
 unsigned num_plot_points=this->nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   // Get local coordinates of plot point
   this->get_s_plot(iplot,nplot,s);

   // Get x position as Vector
   this->get_coordinate_x(s,x);

   for(unsigned i=0;i<DIM;i++)
    {
     fprintf(file_pt,"%g ",x[i]);
    }
   u = interpolated_u_foeppl_von_karman(s);
   fprintf(file_pt,"%g \n",u[0]);//interpolated_w_poisson(s));
  }

 // Write tecplot footer (e.g. FE connectivity lists)
 this->write_tecplot_zone_footer(file_pt,nplot);
}



//======================================================================
 /// Output exact solution
 ///
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points.
 ///
 ///   x,y,u_exact    or    x,y,z,u_exact
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
void FoepplVonKarmanEquations<DIM,NNODE_1D>::output_fct(std::ostream &outfile,
                                       const unsigned &nplot,
                  FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
{
 //Vector of local coordinates
 Vector<double> s(DIM);

  // Vector for coordintes
  Vector<double> x(DIM);

 // Tecplot header info
 outfile << this->tecplot_zone_string(nplot);

 // Exact solution Vector (here a scalar)
 Vector<double> exact_soln(this->required_nvalue(0),0.0);

 // Loop over plot points
 unsigned num_plot_points=this->nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {

   // Get local coordinates of plot point
   this->get_s_plot(iplot,nplot,s);

   // Get x position as Vector
   this->get_coordinate_x(s,x);

   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);

   //Output x,y,...,u_exact
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << x[i] << " ";
    }
   // Loop over variables
   for(unsigned j=0;j<this->required_nvalue(0);j++)
   {
    outfile << exact_soln[j] << " ";
   }
   outfile <<  std::endl;
  }

 // Write tecplot footer (e.g. FE connectivity lists)
 this->write_tecplot_zone_footer(outfile,nplot);
}


//======================================================================
/// Validate against exact solution
///
/// Solution is provided via function pointer.
/// Plot error at a given number of plot points.
///
/// HERE THIS MAY BE SUPERFLUOUS NOW
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
void FoepplVonKarmanEquations<DIM,NNODE_1D>::compute_error_in_deflection(std::ostream &outfile,
                                          FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                                          double& error, double& norm)
{
 // Initialise
 error=0.0;
 norm=0.0;
 //Find out how many nodes there are
 // const unsigned n_node = this->nnode();
 const unsigned n_node_w = 3; // HERE
 //Find out how many bubble nodes there are
 const unsigned n_b_node = this->Number_of_internal_dofs;
 //Find out how many nodes positional dofs there are
 unsigned n_position_type = this->nnodal_position_type();
 // Find the internal dofs
 const unsigned n_b_position_type = this->Number_of_internal_dof_types;

 //Vector of local coordinates
 Vector<double> s(DIM);

 // Vector for coordintes
 Vector<double> x(DIM);

 //Set the value of n_intpt
 unsigned n_intpt = this->integral_pt()->nweight();

 // Exact solution Vector (here a scalar)
 Vector<double> exact_soln(this->required_nvalue(0),0.0);

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {

   //Assign values of s
   for(unsigned i=0;i<DIM;i++)
    {
     s[i] = this->integral_pt()->knot(ipt,i);
    }

   //Get the integral weight
   double w = this->integral_pt()->weight(ipt);

   // Get jacobian of mapping
   double J;
   {
   //Local c1-shape funtion
   Shape psi(n_node_w,n_position_type),test(n_node_w,n_position_type),
    psi_b(n_b_node,n_b_position_type),test_b(n_b_node,n_b_position_type);
   

   DShape dpsi_dxi(n_node_w,n_position_type,DIM),dtest_dxi(n_node_w,n_position_type,DIM),
    dpsi_b_dxi(n_b_node,n_b_position_type,DIM),dtest_b_dxi(n_b_node,n_b_position_type,DIM),
    d2psi_dxi2(n_node_w,n_position_type,3), d2test_dxi2(n_node_w,n_position_type,3),
    d2psi_b_dxi2(n_b_node,n_b_position_type,3), d2test_b_dxi2(n_b_node,n_b_position_type,3);

   J=this-> d2shape_and_d2test_eulerian_foeppl_von_karman(s,
    psi, psi_b, dpsi_dxi, dpsi_b_dxi, d2psi_dxi2, d2psi_b_dxi2,
    test, test_b, dtest_dxi, dtest_b_dxi, d2test_dxi2, d2test_b_dxi2);
   }
   //Premultiply the weights and the Jacobian
   double W = w*J;
   // double Wlin = w*Jlin;

   // Get x position as Vector
   this->get_coordinate_x(s,x);

   // Get FE function value
   Vector<double> u_fe(this->required_nvalue(0),0.0);
   u_fe = interpolated_u_foeppl_von_karman(s);

   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);

   //Output x,y,...,error
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << x[i] << " ";
    }
   for(unsigned ii=0;ii<this->required_nvalue(0);ii++)
    {
     outfile << exact_soln[ii] << " " << exact_soln[ii]-u_fe[ii] << " ";
    }
   outfile << std::endl;

   // Loop over variables
   double tmp1 = 0.0, tmp2 =0.0;
   for(unsigned ii=0;ii<1;ii++)
    {
     // Add to error and norm
     tmp1 = (exact_soln[ii]*exact_soln[ii]*W);
     tmp2 = ((exact_soln[ii]-u_fe[ii])*(exact_soln[ii]-u_fe[ii])*W);
     norm += tmp1;
     error += tmp2;
    }
  } //End of loop over integration pts
}
//HERE OVERLOAD COMPUTE_ERROR SO WE CAN GET A VECTOR OF NORMS


//======================================================================
 /// Validate against exact solution
 ///
 /// Solution is provided via function pointer.
 /// Plot error at a given number of plot points.
 ///
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
void FoepplVonKarmanEquations<DIM,NNODE_1D>::compute_error(std::ostream &outfile,
                                          FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                                          double& error, double& norm)
{
 // Initialise
 error=0.0;
 norm=0.0;

 //Find out how many nodes there are
 // const unsigned n_node = this->nnode();
 const unsigned n_node_w = 3; // HERE
 //Find out how many bubble nodes there are
 const unsigned n_b_node = this->Number_of_internal_dofs;
 //Find out how many nodes positional dofs there are
 unsigned n_position_type = this->nnodal_position_type();
 // Find the internal dofs
 const unsigned n_b_position_type = this->Number_of_internal_dof_types;

 //Vector of local coordinates
 Vector<double> s(DIM);

 // Vector for coordintes
 Vector<double> x(DIM);

 //Set the value of n_intpt
 unsigned n_intpt = this->integral_pt()->nweight();

 // Tecplot
 //outfile << "ZONE" << std::endl;

 // Exact solution Vector (here a scalar)
 Vector<double> exact_soln(this->required_nvalue(0),0.0);

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {

   //Assign values of s
   for(unsigned i=0;i<DIM;i++)
    {
     s[i] = this->integral_pt()->knot(ipt,i);
    }

   //Get the integral weight
   double w = this->integral_pt()->weight(ipt);

   // Get jacobian of mapping
   double J;
   // double Jlin = this->J_eulerian1(s);// Nope
   {
   //Local c1-shape funtion
   Shape psi(n_node_w,n_position_type),test(n_node_w,n_position_type),
    psi_b(n_b_node,n_b_position_type),test_b(n_b_node,n_b_position_type);
   

   DShape dpsi_dxi(n_node_w,n_position_type,DIM),dtest_dxi(n_node_w,n_position_type,DIM),
    dpsi_b_dxi(n_b_node,n_b_position_type,DIM),dtest_b_dxi(n_b_node,n_b_position_type,DIM),
    d2psi_dxi2(n_node_w,n_position_type,3), d2test_dxi2(n_node_w,n_position_type,3),
    d2psi_b_dxi2(n_b_node,n_b_position_type,3), d2test_b_dxi2(n_b_node,n_b_position_type,3);


   J=this-> d2shape_and_d2test_eulerian_foeppl_von_karman(s,
    psi, psi_b, dpsi_dxi, dpsi_b_dxi, d2psi_dxi2, d2psi_b_dxi2,
    test, test_b, dtest_dxi, dtest_b_dxi, d2test_dxi2, d2test_b_dxi2);
   }
   //Premultiply the weights and the Jacobian
   double W = w*J;
   // double Wlin = w*Jlin;

   // Get x position as Vector
   this->get_coordinate_x(s,x);

   // Get FE function value
   Vector<double> u_fe(this->required_nvalue(0),0.0);
   u_fe = interpolated_u_foeppl_von_karman(s);

   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);

   // Loop over variables
   double tmp1 = 0.0, tmp2 =0.0;
   if(Error_metric_fct_pt==0)
   {
   //Output x,y,...,error
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << x[i] << " ";
    }
   for(unsigned ii=0;ii<this->required_nvalue(0);ii++)
    {
     outfile << exact_soln[ii] << " " << exact_soln[ii]-u_fe[ii] << " ";
    }
   outfile << std::endl;

   for(unsigned ii=0;ii<1;ii++)
    {
     // Add to error and norm
     tmp1 = (exact_soln[ii]*exact_soln[ii]);
     tmp2 = ((exact_soln[ii]-u_fe[ii])*(exact_soln[ii]-u_fe[ii]));
    }
   }
   else
    {
     // Get the metric
     (*Error_metric_fct_pt)(x,u_fe,exact_soln,tmp2,tmp1);
    }
   norm += tmp1*W;
   error += tmp2*W;
  } //End of loop over integration pts
}


template class FoepplVonKarmanEquations<2,2>;

template class FoepplVonKarmanEquations<2,3>;

template class FoepplVonKarmanEquations<2,4>;

}
