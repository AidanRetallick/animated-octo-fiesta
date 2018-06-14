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
#include <fenv.h> 
//Generic routines
#include "generic.h" 

// The equations
#include "foeppl_von_karman.h"

// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;
using MathematicalConstants::Pi;

// Utility functions
namespace utility
{
// Template class to read in csv
template <class numeric_type>
void read_csv(ifstream& is, Vector<numeric_type>& data)
 {
  // Loop over lines
  while(!is.eof())
   {
    // Get line
    std::string line;
    getline(is,line,','); 
    // Convert to a stream
    std::stringstream ss;
    //  ss<<std::setprecision(20);
    ss<< line;
    // Put into numeric container
    numeric_type numeric_line;
    ss>>numeric_line;
    
    // Put it in
    data.push_back(numeric_line);
   } 
 }

}
namespace TestSoln
{
//Shape of the domain
double A = 1.0;
double B = 1.0;
// The coupling of the stretching energy
double eta = 1;
double p_mag = 1; 
double nu = 0.5;

// Function protoype
void get_approximate_w(const Vector<double>& x, Vector<double>& w);

//Dummy exact function - does nothing
void get_exact_w(const Vector<double>& xi, Vector<double>& w)
{}

//// Parametric function for boundary part 0
//Vector<double> parametric_edge_0(const double& s)
//{Vector<double> x(2,0.0);   x[0] =-std::sin(s);  x[1] =std::cos(s); 
// return x;};
//
//// Derivative of parametric function
//Vector<double> d_parametric_edge_0(const double& s)
//{Vector<double> dx(2,0.0); dx[0] =-std::cos(s);  dx[1] =-std::sin(s);
//  return dx;};
//
//// Parametric function for boundary part 1
//Vector<double> parametric_edge_1(const double& s)
//{Vector<double> x(2,0.0);   x[0] = std::sin(s);  x[1] =-std::cos(s); 
// return x;};
//
//// Derivative of parametric function
//Vector<double> d_parametric_edge_1(const double& s)
//{Vector<double> dx(2,0.0); dx[0] = std::cos(s);  dx[1] = std::sin(s);
//  return dx;};

// Parametric function for boundary part 0
void parametric_edge_0(const double& s, Vector<double>& x)
 { x[0] =-std::sin(s);  x[1] = std::cos(s);}
// Derivative of parametric function
void d_parametric_edge_0(const double& s, Vector<double>& dx)
 { dx[0] =-std::cos(s);  dx[1] =-std::sin(s);}
// Derivative of parametric function
void d2_parametric_edge_0(const double& s, Vector<double>& dx)
 { dx[0] = std::sin(s);  dx[1] =-std::cos(s);}

// Parametric function for boundary part 1
void parametric_edge_1(const double& s, Vector<double>& x)
{ x[0] = std::sin(s);  x[1] =-std::cos(s);}
// Derivative of parametric function
void  d_parametric_edge_1(const double& s, Vector<double>& dx)
{ dx[0] = std::cos(s);  dx[1] = std::sin(s);};
// Derivative of parametric function
void  d2_parametric_edge_1(const double& s, Vector<double>& dx)
{ dx[0] =-std::sin(s);  dx[1] = std::cos(s);};
// Get s from x
double get_s_0(const Vector<double>& x)
{
// The arc length (parametric parameter) for the upper semi circular arc
 return atan2(-x[0],x[1]);
}

// Get s from x
double get_s_1(const Vector<double>& x)
{
// The arc length (parametric parameter) for the lower semi circular arc
return atan2(x[0],-x[1]);
}

// Assigns the value of pressure depending on the position (x,y)
void get_pressure(const Vector<double>& x, double& pressure)
{
 pressure = p_mag; //constant pressure
}

// Assigns the value of pressure depending on the position (x,y)
void get_pressure_gradient(const Vector<double>& X, Vector<double>& grad)
{
 // Convenient definitions
 double x=X[0], y=X[1]; 
 // Return the (constant) pressure
 grad[0] = 0;//(x*(-1 + 16*x*x - 48*pow(x,4))*eta)/(-1 + nu*nu);
 grad[1] = 0;
}

// The normal and tangential directions.
void get_normal_and_tangent(const Vector<double>& x, Vector<double>& n, 
 Vector<double>& t, DenseMatrix<double>& Dn, DenseMatrix<double>& Dt)
{
 // Fill in the normal and derivatives of the normal
 n[0] = x[0]/sqrt(x[0]*x[0]+x[1]*x[1]);
 n[1] = x[1]/sqrt(x[0]*x[0]+x[1]*x[1]);

 // The derivatives of the x and y components
 //  Dn(0,0) = 1.0/sqrt(x[0]*x[0]+x[1]*x[1])-x[0]*x[0]*pow(x[0]*x[0]+x[1]*x[1],-1.5);
 //  Dn(0,1) =-x[0]*x[1]*pow(x[0]*x[0]+x[1]*x[1],-1.5);
 //  Dn(1,0) =-x[0]*x[1]*pow(x[0]*x[0]+x[1]*x[1],-1.5);
 //  Dn(1,1) = 1.0/sqrt(x[0]*x[0]+x[1]*x[1])-x[1]*x[1]*pow(x[0]*x[0]+x[1]*x[1],-1.5);

 Dn(0,0) = x[1]*x[1] * pow(x[0]*x[0]+x[1]*x[1],-1.5);
 Dn(1,0) =-x[1]*x[0] * pow(x[0]*x[0]+x[1]*x[1],-1.5);
 Dn(0,1) =-x[0]*x[1] * pow(x[0]*x[0]+x[1]*x[1],-1.5);
 Dn(1,1) = x[0]*x[0] * pow(x[0]*x[0]+x[1]*x[1],-1.5);

  // Fill in the tangent and derivatives of the tangent
 t[0] =-x[1]/sqrt(x[0]*x[0]+x[1]*x[1]);
 t[1] = x[0]/sqrt(x[0]*x[0]+x[1]*x[1]);

 Dt(0,0) =-Dn(1,0);
 Dt(1,0) = Dn(0,0); 
 Dt(0,1) =-Dn(1,1);
 Dt(1,1) = Dn(0,1);
}

// // Exact solution for constant pressure, circular domain and resting boundary
// // conditions
// void get_exact_w(const Vector<double>& x, Vector<double>& w)
// {
// //solution (r^4-br^2+c)/64 for w=0 and w'=0 or mrt=0
// w[0]= p_mag*pow(x[0]*x[0]+x[1]*x[1]-1,2)/64;
// w[1]= x[0]*p_mag*(x[0]*x[0]+x[1]*x[1]-1)/16;
// w[2]= x[1]*p_mag*(x[0]*x[0]+x[1]*x[1]-1)/16;
// w[3]= p_mag*(3*x[0]*x[0]+x[1]*x[1]-1)/16;
// w[4]= p_mag*x[0]*x[1]/8;
// w[5]= (x[0]*x[0]+3*x[1]*x[1]-1)/16;
// }

// Exact solution for constant pressure, circular domain and resting boundary
// conditions
void get_exact_w_radial(const Vector<double>& x, Vector<double>& w)
{
//solution (r^4-br^2+c)/64 for w=0 and w'=0 or mrt=0
w[0]= p_mag*pow(x[0]*x[0]+x[1]*x[1]-1,2)/64;
w[1]= p_mag*(x[0]*x[0]+x[1]*x[1]-1)*sqrt(x[0]*x[0]+x[1]*x[1])/16.;
w[2]= 0.0;
w[3]= p_mag*(3*x[0]*x[0]+3*x[1]*x[1]-1)/16.;
w[4]= p_mag*x[0]*x[1]/8;
w[5]= (x[0]*x[0]+3*x[1]*x[1]-1)/16;
}

void error_metric(const Vector<double>& x, const 
  Vector<double>& u, const Vector<double>& u_exact, double& error, double& norm)
{
 // We use the theta derivative of the out of plane deflection
 error = pow(-x[1]*u[1] + x[0]*u[2],2);
 norm  = pow( x[0]*u[1] + x[1]*u[2],2) ;
}
// Get the exact solution
void get_null_fct(const Vector<double>& X, double& exact_w)
{
 exact_w = 0.0;
}

}

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


//==start_of_problem_class============================================
/// Class definition
//====================================================================
template<class ELEMENT>
class UnstructuredFvKProblem : public virtual Problem
{

public:

/// Constructor
UnstructuredFvKProblem(double element_area = 0.1);

/// Destructor
~UnstructuredFvKProblem()
{
};

/// Update after solve (empty)
void actions_after_newton_solve()
{
}

/// Update the problem specs before solve: Re-apply boundary conditions
/// Empty as the boundary conditions stay fixed
void actions_before_newton_solve()
{
apply_boundary_conditions();
}

/// Doc the solution
void doc_solution(const std::string& comment="");

/// \short Overloaded version of the problem's access function to
/// the mesh. Recasts the pointer to the base Mesh object to
/// the actual mesh type.
TriangleMesh<ELEMENT>* mesh_pt()
{
return dynamic_cast<TriangleMesh<ELEMENT>*> (Problem::mesh_pt()); 
}

/// Doc info object for labeling output
DocInfo Doc_info;
private:

void actions_after_read_unstructured_meshes()
 {
 // Curved Edge upgrade
 upgrade_edge_elements_to_curve(0,Bulk_mesh_pt);
 upgrade_edge_elements_to_curve(1,Bulk_mesh_pt);
  
 // Rotate degrees of freedom
 rotate_edge_degrees_of_freedom(Bulk_mesh_pt);
  complete_problem_setup();
  apply_boundary_conditions();
 } 

/// Helper function to apply boundary conditions
void apply_boundary_conditions();

/// \short Helper function to (re-)set boundary condition
/// and complete the build of  all elements
void complete_problem_setup();

/// Trace file to document norm of solution
ofstream Trace_file;

// Keep track of boundary ids
enum
{
 Outer_boundary0 = 0,
 Outer_boundary1 = 1,
 Inner_boundary0 = 2,
 Inner_boundary1 = 3
};

double Element_area;

// The extra members required for flux type boundary conditions.
/// \short Number of "bulk" elements (We're attaching flux elements
/// to bulk mesh --> only the first Nkirchhoff_elements elements in
/// the mesh are bulk elements!)
// unsigned Nkirchhoff_elements;

/// \short Create bending moment elements on the b-th boundary of the
/// problems mesh 
void create_traction_elements(const unsigned &b, Mesh* const & bulk_mesh_py,
                            Mesh* const &surface_mesh_pt);

void upgrade_edge_elements_to_curve(const unsigned &b, Mesh* const & 
 bulk_mesh_pt);

void rotate_edge_degrees_of_freedom(Mesh* const &bulk_mesh_pt);


/// \short Delete traction elements and wipe the surface mesh
void delete_traction_elements(Mesh* const &surface_mesh_pt);

/// \short Set pointer to prescribed-flux function for all elements
/// in the surface mesh
void set_prescribed_traction_pt();

/// Pointer to "bulk" mesh
TriangleMesh<ELEMENT>* Bulk_mesh_pt;

/// Pointer to "surface" mesh
Mesh* Surface_mesh_pt;

private :
/// Finite difference bits
 unsigned N_fd=800;
 Vector<double> w_fd=Vector<double>(2*N_fd,0.0);

public :
void get_axisymmetric_solution(const double& r, double& w)
 {
  // Const reference to nu
  const double& nu = TestSoln::nu;
  // What point in the domain
  if(r<=1)
  {
   // Positions and octave index
   unsigned nlower = std::floor((N_fd -1)*r +1);
   double rlower = double(nlower-1)/double(N_fd -1);
   unsigned nupper = std::ceil((N_fd -1)*r +1);
   double rupper = double(nupper-1)/double(N_fd -1);
   // Convert to local coordinate s
   if(nlower != nupper)
    {
     // Get the local coordinate s
     const double s((r - rlower) / (rupper - rlower));
     // Linear interpolation, and scaling
     w = (w_fd[nlower-1]*s + w_fd[nupper-1]*(1-s))/sqrt(12.*(1-nu*nu));
    }
   // They are the same point
   else 
    {
     // NB octave uses indexing from one. Just compute and then scale
     w = w_fd[nlower-1]/sqrt(12.*(1-nu*nu));
    }
  }
   else
   {
    oomph_info <<"r greater than 1 detected: r="<<r <<"\n";
   }
 }

void fill_in_axisymmetric_solution()
{
 // Local reference to nu and converted pressure
 const double& nu = TestSoln::nu;
 const double p = TestSoln::p_mag*std::pow(12.*(1-nu*nu),1.5);
 // The command for octave
 char octave_command[1000];
 char csvname[100]="fd_octave_solve.reslt";
 const double p_step =2000;
 const unsigned max_iter =400;

 // Set up the command
 oomph_info<<"Running octave-cli for p_fd="<<p<<" or p="<<TestSoln::p_mag<<"\n";
 sprintf(octave_command,"octave-cli --eval\
 \"info=get_solution_at_pressure(%f,%f,\'%s\',%i,%i,%f);\
  printf('Returned with value %%i.\\n',info)\"",p,nu,csvname,N_fd,max_iter,p_step);
 oomph_info <<"Running: "<< octave_command <<std::endl;
 // Run it
 int sysret=system(octave_command);
 // Now check octave_out.txt 
 Vector<double> data;
 ifstream istream;
 istream.open(csvname);

 oomph_info<<"Reading in data from csv file.\n"; 
 utility::read_csv(istream,data);

 // Check the length
 if (data.size()!= 2*N_fd)
  {
   // Scream
   oomph_info << "Data read longer than expected. Something has gone wrong.\n"
              << "Length of data: "<< data.size()<<" expected length :"
              << 2*N_fd<<std::endl;
   exit(-1);
  }
 else
  {
   // Deep copy
   oomph_info << "Doing deep copy of data\n";
   w_fd=data; 
  }
 
 istream.close();

}

}; // end_of_problem_class


template<class ELEMENT>
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem(double element_area)
:
Element_area(element_area)
{
Vector<double> zeta(1);
Vector<double> posn(2);

//Outer boundary
//--------------

double A = 1.0;
double B = 1.0;
Ellipse* outer_boundary_ellipse_pt = new Ellipse(A, B);

TriangleMeshClosedCurve* outer_boundary_pt = 0;

Vector<TriangleMeshCurveSection*> outer_curvilinear_boundary_pt(2);

//First bit
double zeta_start = 0.0;
double zeta_end = MathematicalConstants::Pi;
unsigned nsegment = (int)(MathematicalConstants::Pi/sqrt(element_area));
outer_curvilinear_boundary_pt[0] =
new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
zeta_end, nsegment, Outer_boundary0);

//Second bit
zeta_start = MathematicalConstants::Pi;
zeta_end = 2.0*MathematicalConstants::Pi;
nsegment = (int)(MathematicalConstants::Pi/sqrt(element_area));
outer_curvilinear_boundary_pt[1] =
new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
zeta_end, nsegment, Outer_boundary1);

outer_boundary_pt =
new TriangleMeshClosedCurve(outer_curvilinear_boundary_pt);

// Internal open boundaries
// Total number of open curves in the domain
unsigned n_open_curves = 2;
// We want internal open curves
Vector<TriangleMeshOpenCurve *> inner_open_boundaries_pt(n_open_curves);

// Internal bit - this means we can have a boundary which is just the centre
// We start by creating the internal boundaries
// The boundary 2 is defined by its two vertices
// Open curve 1
 Vector<Vector<double> > vertices(2,Vector<double>(2,0.0));
 vertices[0][0] =-0.5;
 vertices[0][1] = 0.0;
 
 vertices[1][0] = 0.5;
 vertices[1][1] = 0.0;
 unsigned boundary_id = Inner_boundary0;

 TriangleMeshPolyLine *boundary2_pt =
   new TriangleMeshPolyLine(vertices, boundary_id);
// Open Curve 2
 vertices[0][0] = 0.0;
 vertices[0][1] =-0.5;
 
 vertices[1][0] = 0.0;
 vertices[1][1] = 0.5;
 boundary_id = Inner_boundary1;

 TriangleMeshPolyLine *boundary3_pt =
   new TriangleMeshPolyLine(vertices, boundary_id);

// Each internal open curve is defined by a vector of
// TriangleMeshCurveSection,
// on this example we only need one curve section for each internal boundary
 Vector<TriangleMeshCurveSection *> internal_curve_section1_pt(1);
 internal_curve_section1_pt[0] = boundary2_pt;

 Vector<TriangleMeshCurveSection *> internal_curve_section2_pt(1);
 internal_curve_section2_pt[0] = boundary3_pt;

// The open curve that define this boundary is composed of just one
// curve section
 inner_open_boundaries_pt[0] =
    new TriangleMeshOpenCurve(internal_curve_section1_pt);

 inner_open_boundaries_pt[1] =
    new TriangleMeshOpenCurve(internal_curve_section2_pt);

//Create the mesh
//---------------
//Create mesh parameters object
TriangleMeshParameters mesh_parameters(outer_boundary_pt);

mesh_parameters.element_area() = element_area;

// Specify the internal open boundaries
mesh_parameters.internal_open_curves_pt() = inner_open_boundaries_pt;

// Build an assign bulk mesh
Bulk_mesh_pt=new TriangleMesh<ELEMENT>(mesh_parameters);

// Create "surface mesh" that will contain only the prescribed-traction
// elements. The constructor creates the mesh without adding any nodes
// elements etc.
Surface_mesh_pt =  new Mesh;

//Add two submeshes to problem
add_sub_mesh(Bulk_mesh_pt);
add_sub_mesh(Surface_mesh_pt);

// Combine submeshes into a single Mesh
build_global_mesh();

// Curved Edge upgrade
upgrade_edge_elements_to_curve(0,Bulk_mesh_pt);
upgrade_edge_elements_to_curve(1,Bulk_mesh_pt);
 
// Rotate degrees of freedom
rotate_edge_degrees_of_freedom(Bulk_mesh_pt);

// Store number of bulk elements
complete_problem_setup();

char filename[100];
sprintf(filename, "RESLT/trace.dat");
Trace_file.open(filename);

oomph_info << "Number of equations: "
        << assign_eqn_numbers() << '\n';
}



//==start_of_complete======================================================
/// Set boundary condition exactly, and complete the build of 
/// all elements
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::complete_problem_setup()
{   

// Set the boundary conditions for problem: All nodes are
// free by default -- just pin the ones that have Dirichlet conditions
// here. 

// Pin the node that is at the centre in the domain!
// Get the num of nods on internal_boundary 0
unsigned num_int_nod=Bulk_mesh_pt->nboundary_node(2);
for (unsigned inod=0;inod<num_int_nod;inod++)
{
 // Get node point
 Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(2,inod);
 // If the node is on the other internal boundary too
 if( nod_pt->is_on_boundary(3))
 {
  oomph_info<<"Found centre point\n";

  oomph_info<<"x = ("<<nod_pt->x(0)<<" "<<nod_pt->x(1)<<") \n";
  // Pin it! It's the centre of the domain!
  nod_pt->pin(0);
  nod_pt->set_value(6,0.0);
  nod_pt->pin(1);
  nod_pt->set_value(7,0.0);
 }
}
//Just loop over outer boundary since inner boundary doesn't have boundary
//conditions
unsigned nbound = Outer_boundary1 + 1;
for(unsigned ibound=0;ibound<nbound;ibound++)
{
//  Node* first_nod_pt=Bulk_mesh_pt->node_pt(0);
 unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
for (unsigned inod=0;inod<num_nod;inod++)
{
// // Get first node and pin the inplane disp
 Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
// first_nod_pt->pin(6);
// first_nod_pt->set_value(6,0.0);
// first_nod_pt->pin(7);
// first_nod_pt->set_value(7,0.0);

 // Now for the rest
 Vector<double> x(2,0.0),w(6,0.0);
 x[0]=nod_pt->x(0);
 x[1]=nod_pt->x(1);

// TestSoln::get_exact_w_radial(x,w);
// // Pin unknown values (everything except for the second normal derivative)
// nod_pt->pin(2+0);
// nod_pt->set_value(0,0.0);
// nod_pt->pin(2+1);
// nod_pt->set_value(1,0.0);
// nod_pt->pin(2+2);
// nod_pt->set_value(2,0.0);
// nod_pt->pin(2+4);
// nod_pt->set_value(4,0.0);
// nod_pt->pin(2+5);
// nod_pt->set_value(5,0.0);
 }
} // end loop over boundaries 


// Complete the build of all elements so they are fully functional
unsigned n_element = Bulk_mesh_pt->nelement();
for(unsigned e=0;e<n_element;e++)
{
// Upcast from GeneralisedElement to the present element
ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

//Set the pressure function pointers and the physical constants
el_pt->pressure_fct_pt() = &TestSoln::get_pressure;
el_pt->pressure_fct_gradient_pt() = &TestSoln::get_pressure_gradient;   
el_pt->error_metric_fct_pt() = &TestSoln::error_metric;
el_pt->nu_pt() = &TestSoln::nu;
el_pt->eta_pt() = &TestSoln::eta;
//el_pt->pin_all_deflection_dofs();
}

// Set the boundary conditions
for(unsigned b=0;b<2;++b)
 {
 const unsigned nb_element = Bulk_mesh_pt->nboundary_element(b);
 for(unsigned e=0;e<nb_element;e++)
  {
   // Get pointer to bulk element adjacent to b
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->boundary_element_pt(b,e));
   el_pt->fix_out_of_plane_displacement_dof(2+0,b,TestSoln::get_null_fct);
   el_pt->fix_out_of_plane_displacement_dof(2+1,b,TestSoln::get_null_fct);
   el_pt->fix_out_of_plane_displacement_dof(2+2,b,TestSoln::get_null_fct);
   el_pt->fix_out_of_plane_displacement_dof(2+4,b,TestSoln::get_null_fct);
   el_pt->fix_out_of_plane_displacement_dof(2+5,b,TestSoln::get_null_fct);
  //  el_pt->fix_in_plane_displacement_dof(0,b,TestSoln::get_null_fct);
  //  el_pt->fix_in_plane_displacement_dof(1,b,TestSoln::get_null_fct);
  }
 }
// Loop over flux elements to pass pointer to prescribed traction function

/// Set pointer to prescribed traction function for traction elements
//set_prescribed_traction_pt();

// Re-apply Dirichlet boundary conditions (projection ignores
// boundary conditions!)
apply_boundary_conditions();
}

//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::apply_boundary_conditions()
{

// Loop over all boundary nodes
//Just loop over outer boundary conditions
unsigned nbound = Outer_boundary1 + 1;

for(unsigned ibound=0;ibound<nbound;ibound++)
{
unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
for (unsigned inod=0;inod<num_nod;inod++)
{
 // Get node
 Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
 
 // Extract nodal coordinates from node:
 Vector<double> x(2);
 x[0]=nod_pt->x(0);
 x[1]=nod_pt->x(1);
}
} 

} // end set bc

template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::
create_traction_elements(const unsigned &b, Mesh* const &bulk_mesh_pt, 
                             Mesh* const &surface_mesh_pt)
{
}// end create traction elements

template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::
upgrade_edge_elements_to_curve(const unsigned &b, Mesh* const &bulk_mesh_pt) 
{
 // How many bulk elements adjacent to boundary b
 unsigned n_element = bulk_mesh_pt-> nboundary_element(b);
 
 // These depend on the boundary we are on
 void (*parametric_edge_fct_pt)(const double& s, Vector<double>&x);
 void (*d_parametric_edge_fct_pt)(const double& s, Vector<double> &dx);
 double (*get_arc_position)(const Vector<double>& s);
 
// Define the functions for each part of the boundary
 switch (b)
  {
   // Upper boundary
   case 0:
    parametric_edge_fct_pt = &TestSoln::parametric_edge_0;
    d_parametric_edge_fct_pt = &TestSoln::d_parametric_edge_0;
    get_arc_position = &TestSoln::get_s_0;
   break;

   // Lower boundary
   case 1:
    parametric_edge_fct_pt = &TestSoln::parametric_edge_1;
    d_parametric_edge_fct_pt = &TestSoln::d_parametric_edge_1;
    get_arc_position = &TestSoln::get_s_1;
   break;

   default:
    throw OomphLibError(
     "I have encountered a boundary number that I wasn't expecting. This is very\
 peculiar.",
     "UnstructuredFvKProblem::upgrade_edge_elements_to_curve(...)",
     OOMPH_EXCEPTION_LOCATION);
   break;
  }
 
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0; e<n_element; e++)
  {
   // Get pointer to bulk element adjacent to b
   ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));
   
   // Loop over nodes
   unsigned nnode=3;//bulk_el_pt->nnode();
   unsigned index_of_interior_node=3;

   // The edge that is curved
   MyC1CurvedElements::BernadouElementBasis<3>::Edge edge;

   // Vertices positions
   Vector<Vector<double> > xn(3,Vector<double>(2,0.0));
 
   // Get vertices for debugging
   Vector<Vector<double> > verts(3,Vector<double>(2,0.0));
   // Loop nodes
   for(unsigned n=0;n<nnode;++n)
    {
     // If it is on boundary
     Node* nod_pt = bulk_el_pt->node_pt(n);
     verts[n][0]=nod_pt->x(0);
     verts[n][1]=nod_pt->x(1);
     if(nod_pt->is_on_boundary())
      {
       xn[n][0]=nod_pt->x(0);
       xn[n][1]=nod_pt->x(1);
      }
     // The edge is denoted by the index of the  opposite (interior) node
     else {index_of_interior_node = n;}
    }
   // Initialise s_ubar s_obar
   double s_ubar, s_obar;

   // s at the next (cyclic) node after interior
   s_ubar = (*get_arc_position)(xn[(index_of_interior_node+1) % 3]);
   // s at the previous (cyclic) node before interior
   s_obar = (*get_arc_position)(xn[(index_of_interior_node+2) % 3]);

   // Assign edge case
   switch(index_of_interior_node)
    {
     case 0: edge= MyC1CurvedElements::BernadouElementBasis<3>::zero; 
      break;
     case 1: edge= MyC1CurvedElements::BernadouElementBasis<3>::one; 
      break;
     case 2: edge= MyC1CurvedElements::BernadouElementBasis<3>::two; 
      break;
     // Should break it here HERE
     default: edge= MyC1CurvedElements::BernadouElementBasis<3>::none; 
      throw OomphLibError(
       "The edge number has been set to a value greater than two: either we have\
 quadrilateral elements or more likely the index_of_interior_node was never set\
 and remains at its default value.",
       "UnstructuredFvKProblem::upgrade_edge_elements_to_curve(...)",
       OOMPH_EXCEPTION_LOCATION);
      break;
     }
   if (s_ubar>s_obar)
    {std::cout<<s_ubar<<" "<<s_obar<<"\n";}
   // Check for inverted elements HERE

   // Upgrade it
    bulk_el_pt->upgrade_to_curved_element(edge,s_ubar,s_obar,
     *parametric_edge_fct_pt,*d_parametric_edge_fct_pt);
    
   // Get vertices for debugging
   Vector<Vector<double> > lverts(3,Vector<double>(2,0.0));
   lverts[0][0]=1.0;
   lverts[1][1]=1.0;
   Vector<Vector<double> > fkverts(3,Vector<double>(2,0.0));
   bulk_el_pt->get_coordinate_x(lverts[0],fkverts[0]);
   bulk_el_pt->get_coordinate_x(lverts[1],fkverts[1]);
   bulk_el_pt->get_coordinate_x(lverts[2],fkverts[2]);

  }
}// end upgrade elements

template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::
rotate_edge_degrees_of_freedom( Mesh* const &bulk_mesh_pt)
{
 // How many bulk elements
 unsigned n_element = bulk_mesh_pt-> nelement();
 
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0; e<n_element; e++)
  {
   // Get pointer to bulk element adjacent to b
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
 
   // Loop nodes  HERE
   unsigned nnode = 3;//nnode();
   unsigned nbnode=0 ;
   // Count the number of boundary nodes
   for (unsigned n=0; n<nnode;++n)
     {
      // Check it isn't on an internal boundary
      bool on_boundary_2=el_pt->node_pt(n)->is_on_boundary(2);
      bool on_boundary_3=el_pt->node_pt(n)->is_on_boundary(3);
      if(!(on_boundary_2 || on_boundary_3))
       {nbnode+=unsigned(el_pt->node_pt(n)->is_on_boundary());}
     }

   // Now if we have nodes on boundary 
   if(nbnode>0)
    {
     // Set up vector
     Vector<unsigned> bnode (nbnode,0);
     unsigned inode(0);

     // Fill in the bnode Vector
     for (unsigned n=0; n<nnode;++n)
      {
       // Check it isn't on an internal boundary
       bool on_boundary_2=el_pt->node_pt(n)->is_on_boundary(2);
       bool on_boundary_3=el_pt->node_pt(n)->is_on_boundary(3);
       if(!(on_boundary_2 || on_boundary_3))
       {
       // If it is on the boundary
       if(el_pt->node_pt(n)->is_on_boundary())
        {
         // Set up the Vector
         bnode[inode]=n;
         ++inode;
        }
       }
      }
    // Output that we have found element HERE
//    std::cout<<"Element "<<e<<" has "<<bnode<< " nodes on the boundary.\n";

    el_pt->set_up_rotated_dofs(nbnode,bnode,&TestSoln::get_normal_and_tangent);
   // Now rotate the nodes
   }
 }
}// end create traction elements

//==start_of_set_prescribed_traction_pt===================================
/// Set pointer to prescribed traction function for all elements in the 
/// surface mesh
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::set_prescribed_traction_pt()
{
}// end of set prescribed flux pt

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::doc_solution(const 
                                                    std::string& comment)
{ 
// Lets see how this goes
// fill_in_axisymmetric_solution();

ofstream some_file;
char filename[100];

// Number of plot points
unsigned npts = 2;

sprintf(filename,"RESLT/coarse_soln%i-%f.dat",Doc_info.number(),Element_area);
some_file.open(filename);
Bulk_mesh_pt->output(some_file,npts); 
some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" 
       << comment << "\"\n";
some_file.close();

// Number of plot points
npts = 5;

sprintf(filename,"RESLT/soln%i-%f.dat",Doc_info.number(),Element_area);
some_file.open(filename);
Bulk_mesh_pt->output(some_file,npts); 
some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" 
       << comment << "\"\n";
some_file.close();

//  Output exact solution
sprintf(filename,"%s/axisymmetric_soln%i-%f.dat","RESLT",Doc_info.number()
 ,Element_area);
some_file.open(filename);
Bulk_mesh_pt->output_fct(some_file,npts,TestSoln::get_exact_w); 
some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" 
       << comment << "\"\n";
some_file.close();

// Output boundaries
//------------------
sprintf(filename,"RESLT/boundaries%i-%f.dat",Doc_info.number(),Element_area);
some_file.open(filename);
Bulk_mesh_pt->output_boundaries(some_file);
some_file.close();

// Output regions
unsigned n_region = Bulk_mesh_pt->nregion();
if (n_region > 1)
{
for (unsigned r = 0; r < n_region; r++)
{
 //Attempt to output elements in different regions
 sprintf(filename,"RESLT/region%i%i-%f.dat",r,Doc_info.number(),
   Element_area);
 some_file.open(filename);
 unsigned nel = Bulk_mesh_pt->nregion_element(r);
 for (unsigned e = 0; e < nel; e++)
  {
   Bulk_mesh_pt->region_element_pt(r,e)->output(some_file,npts);
  }
 some_file.close();
}
}

// // Doc error and return of the square of the L2 error
// //---------------------------------------------------
// //double error,norm,dummy_error,zero_norm;
  double dummy_error,zero_norm;
 sprintf(filename,"RESLT/error%i-%f.dat",Doc_info.number(),Element_area);
 some_file.open(filename);
 
 Bulk_mesh_pt->compute_error(some_file,TestSoln::get_exact_w,
                        dummy_error,zero_norm);
 some_file.close();
 
 // Doc L2 error and norm of solution
 oomph_info << "Absolute norm of computed solution: " << sqrt(dummy_error) 
            << std::endl;
 
 oomph_info << "Norm of computed solution: " << sqrt(dummy_error/zero_norm)
            << std::endl;
 
 Trace_file << TestSoln::p_mag << " " << "\n ";

// Doc error and return of the square of the L2 error
//---------------------------------------------------
sprintf(filename,"RESLT/L2-norm%i-%f.dat",
        Doc_info.number(),
        Element_area);
some_file.open(filename);

some_file<<"### L2 Norm\n";
some_file<<"##  Format: err^2 norm^2 log(err/norm) \n";
// Print error in prescribed format
some_file<< dummy_error <<" "<< zero_norm <<" ";

// Only divide by norm if its nonzero
some_file<<0.5*(log10(fabs(dummy_error))-log10(zero_norm))<<"\n";
some_file.close();

// Increment the doc_info number
Doc_info.number()++;

} // end of doc

//============start_of_delete_flux_elements==============================
/// Delete Poisson Flux Elements and wipe the surface mesh
//=======================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>
::delete_traction_elements(Mesh* const &surface_mesh_pt)
{
// How many surface elements are in the surface mesh
unsigned n_element = surface_mesh_pt->nelement();

// Loop over the surface elements
for(unsigned e=0;e<n_element;e++)
{
// Kill surface element
delete surface_mesh_pt->element_pt(e);
}

// Wipe the mesh
surface_mesh_pt->flush_element_and_node_storage();

} // end of delete_flux_elements

// Namespace extension
namespace TestSoln{

UnstructuredFvKProblem<FoepplVonKarmanC1CurvedBellElement<2,4,3> >* problem_pt=0;

// Approximate axisymmetric solution: points to a member of problem.
void get_approximate_w(const Vector<double>& x, Vector<double>& w)
 {
 double r = sqrt(x[0]*x[0]+x[1]*x[1]);
 // if not null
 if(problem_pt!=0)
  {
  // Get the axisymmetric solution for the problem
  problem_pt->get_axisymmetric_solution(r,w[0]);
  }
 }
}

//=======start_of_main========================================
///Driver code for demo of inline triangle mesh generation
//============================================================
int main(int argc, char **argv)
{

 feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified

 // Validation?
 CommandLineArgs::specify_command_line_flag("--validation");

// The `restart' flag
 CommandLineArgs::specify_command_line_flag("--restart");
 
// Just do what you 'would' do without doing it
 CommandLineArgs::specify_command_line_flag("--dry-run");

 // Directory for solution
 string output_dir="RESLT";
 CommandLineArgs::specify_command_line_flag("--dir", &output_dir);

 // Poisson Ratio
 CommandLineArgs::specify_command_line_flag("--nu", &TestSoln::nu);
 
 // Applied Pressure
 CommandLineArgs::specify_command_line_flag("--p", &TestSoln::p_mag);

 // P_step
 double p_step=100;
 CommandLineArgs::specify_command_line_flag("--dp", &p_step);

 // Applied Pressure
 double p_max = 4000;
 CommandLineArgs::specify_command_line_flag("--p_max", &p_max);

 // Element Area
 double element_area=0.2;
 CommandLineArgs::specify_command_line_flag("--element_area", &element_area);

 // Parse command line
 CommandLineArgs::parse_and_assign(); 

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Problem instance
 UnstructuredFvKProblem<FoepplVonKarmanC1CurvedBellElement<2,4,3> >
   problem(element_area);
 // Set the Testsoln::problem_pt to be problem
 TestSoln::problem_pt = &problem;

 // If restart flag has been set
 const bool do_restart=CommandLineArgs::
   command_line_flag_has_been_set("--restart");
 const bool do_perturb=CommandLineArgs::
   command_line_flag_has_been_set("--p_pert");
 const bool dry_run=CommandLineArgs::
   command_line_flag_has_been_set("--dry-run");

 problem.max_residuals()=1e9;
 problem.max_newton_iterations()=20;
// problem.newton_solver_tolerance()=1e-9;

 // If we are restarting
 if(do_restart)
 {
  oomph_info<<"Restarting\n";
  // Read the restart file
  if(!dry_run)
  {
  std::ifstream filestream;
  filestream.open("fvk_problem_data.dump");
  problem.read(filestream);    
  filestream.close(); 

  // Check with a newton solve 
  problem.newton_solve(); 
 
  //Output solution
  problem.doc_solution();
  }
 }

// Now loop and solve
 while(TestSoln::p_mag<p_max)
  {
   oomph_info<<"Solving for p=" << TestSoln::p_mag<<"\n";
  if(!dry_run)
  {
   problem.newton_solve();

   // Document
   problem.doc_solution();
   }
   oomph_info << std::endl;
   oomph_info << "---------------------------------------------" << std::endl;
   oomph_info << " Pcos (" << TestSoln::p_mag << ")" << std::endl;
   oomph_info << "Current dp  (" << p_step << ")" << std::endl;
   oomph_info << "Poisson ratio (" << TestSoln::nu << ")" << std::endl;
   oomph_info << "Solution number (" <<problem.Doc_info.number()-1 << ")"
              << std::endl;
   oomph_info << "---------------------------------------------" << std::endl;
   oomph_info << std::endl;

  // Dump the data 
  char filename[100];
  std::ofstream filestream;
  filestream.precision(15);
  sprintf(filename,"%s/fvk_circle_data%i-%f.dump",
          problem.Doc_info.directory().c_str(),
          problem.Doc_info.number(),
          TestSoln::p_mag
         );
  oomph_info<<"\nDumping to: \""<<filename<<"\"\n";
  if(!dry_run)
  {
  filestream.open(filename);
  problem.dump(filestream);
  filestream.close();
  }
  // Increase the pressure 
  TestSoln::p_mag+=p_step;
  }

oomph_info<<"Exiting Normally\n";
} //End of main

