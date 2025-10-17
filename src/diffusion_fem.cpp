#include "inmost.h"
#include <stdio.h>


using namespace INMOST;
using namespace std;

// // Corresponds to tensor
// // [ 1  0 ]
// // [ 0 10 ]
// // rotated by M_PI/6
// const double Dxx = 3.25;
// const double Dyy = -0.433013;
// const double Dxy = 0.25;

const double dx = 1.0;
const double dy = 1.0;
const double dxy = 0.0;
const double pi = 3.1415926535898;
const double a = 1;

double C(double x, double y)
{
    return 1;
	return sin(a*x) * sin(a*y);
}

double source(double x, double y)
{
	return 0;
	return -a*a * (2.*dxy * cos(a*x)*cos(a*y) - (dx+dy) * sin(a*x)*sin(a*y));
}

// Class including everything needed
class Problem
{
private:
	/// Mesh
	Mesh &m;
	// =========== Tags =============
	/// Solution (concentration) tag defined on nodes
	Tag tagConc;
	/// Diffusion tensor tag (3 numbers per cell = triangle)
	Tag tagD;
	/// Boundary condition tag (1 real number per node)
	Tag tagBCval;
	/// Source tag (1 real number per node)
	Tag tagSource;
	/// Analytical solution tag (1 real number per node)
	Tag tagConcAn;
	/// Global index tag (1 integer number per node)
	Tag tagGlobInd;
	/// Numerical error in concentration (1 real number per node)
	Tag tagConcErr;

	// =========== Tag names ===========
	const string tagNameConc = "Concentration"; 
	const string tagNameD = "Diffusion_tensor";
	const string tagNameBCtype = "BC_type";
	const string tagNameBCval = "BC_value";
	const string tagNameSource = "Source";
	const string tagNameConcAn = "Concentration_analytical";
	const string tagNameGlobInd = "Global_Index";
	const string tagNameConcErr = "Error";

	// =========== Markers
	/// Marker for Dirichlet nodes
	MarkerType mrkDirNode;
	/// Number of Dir. nodes. Total number of unknowns = number of nodes - numDirNodes
	unsigned numDirNodes;
public:
	Problem(Mesh &m_);
	~Problem();

	// Creation of tags, markup of Dirichlet nodes
	void initProblem();

	// Assembly of the global linear system
	void assembleGlobalSystem(Sparse::Matrix &A, Sparse::Vector &rhs);

	// Assembly of local system of size 3x3
	void assembleLocalSystem(const Cell &c, rMatrix &A_loc, rMatrix &rhs_loc);

	// Solution process
	void run();

	// Basis piecewise-linear function
	double basis_func(const Cell& c, const Node& n, double x_, double y_);

	// Its gradient (= const on cell, so no (x,y))
	rMatrix basis_func_grad(const Cell& c, const Node& n);

	// Transform barycentric coordinates to (x,y)
	void coords_from_barycentric(double* node_x, double* node_y, double* eta, double* x, double* y);

	// Approximate in a cell a function defined by its nodal values in tag T  
	double approximate(const Cell& c, const Tag& T, double x, double y);

	// Integrate a function defined by its nodal values in tag T
	double integrate_over_triangle(const Cell& c, const Tag& T);

	// Get error C-norm
    double get_err_c_norm();

	// Get error L2-norm
    double get_err_L2_norm();
};

Problem::Problem(Mesh &m_) : m(m_)
{
}

Problem::~Problem()
{

}

void Problem::initProblem()
{
	// Init tags
	tagConc    = m.CreateTag(tagNameConc,    DATA_REAL,    NODE, NONE, 1); // 1 real number per node
	tagD       = m.CreateTag(tagNameD,       DATA_REAL,    CELL, NONE, 3); // 3 real numbers per cell
	tagBCval   = m.CreateTag(tagNameBCval,   DATA_REAL,    NODE, NODE, 1); // 1 real number per node, sparse on nodes
	tagSource  = m.CreateTag(tagNameSource,  DATA_REAL,    NODE, NONE, 1); // 1 real number per node
	tagConcAn  = m.CreateTag(tagNameConcAn,  DATA_REAL,    NODE, NONE, 1); // 1 real number per node
	tagGlobInd = m.CreateTag(tagNameGlobInd, DATA_INTEGER, NODE, NONE, 1); // 1 integer number per node
	tagConcErr = m.CreateTag(tagNameConcErr, DATA_REAL,    NODE, NONE, 1); // 1 real number per node

	// Cell loop
	// 1. Check that all cells are triangles
	// 2. Fill tensor
	for(Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++){
		Cell c = icell->getAsCell();
		ElementArray<Node> nodes = c.getNodes();
		if (nodes.size() != 3){
			cout<<"Cell is not a triangle!!!";
			exit(1);
		}
		
		c.RealArray(tagD)[0] = dx;
		c.RealArray(tagD)[1] = dy;
		c.RealArray(tagD)[2] = dxy;
	}

	// Node loop
	// 1. Set source and analytics
	// 2. Check boundary, set marker and BC tag
	// 3. Set global indices
	mrkDirNode = m.CreateMarker();
	numDirNodes = 0;
	int glob_ind = 0;
	for(Mesh::iteratorNode inode = m.BeginNode(); inode != m.EndNode(); inode++){
		Node n = inode->getAsNode();
		double xn[3];
		n.Centroid(xn);
		n.Real(tagConcAn) = C(xn[0], xn[1]);
		n.Real(tagSource) = source(xn[0], xn[1]);

		if(n.Boundary()){
			n.SetMarker(mrkDirNode);
			numDirNodes++;
			n.Real(tagBCval) = n.Real(tagConcAn);
		}
		else{
			n.Integer(tagGlobInd) = glob_ind;
			glob_ind++;
		}
		n.Real(tagConcErr) = 0.0;
	}
	printf("Number of Dirichlet nodes: %d\n", numDirNodes);
}

// Piecewise linear "pyramid" basis function, more precisely, its restricition to 1 triangle
// Since a triangle has 3 such function, a node is passed to determine specific function
double Problem::basis_func(const Cell &c, const Node &n, double x_, double y_)
{
    ElementArray<Node> nodes = c.getNodes();
	unsigned n_ind = 0;
	double x[3];
	double y[3];
	for(unsigned i = 0; i < 3; i++){
		if(n == nodes[i])
			n_ind = i;
		double coords[3];
		nodes[i].Centroid(coords);
		x[i] = coords[0];
		y[i] = coords[1];
	}
	
	if(n_ind == 0){
		return ((x_   - x[2])*(y[1] - y[2]) - (x[1] - x[2])*(y_   - y[2])) /
			   ((x[0] - x[2])*(y[1] - y[2]) - (x[1] - x[2])*(y[0] - y[2]));
	}
	else if(n_ind == 1){
		// implement yourself
	}
	else if(n_ind == 2){
		// implement yourself
	}
	else{
		printf("Unexpected n_ind = %d\n", n_ind);
		exit(1);
	}
}

// Gradient of the basis_func
// Since it's a vector, it's returned as INMOST::rMatrix (matrix of reals)
// Is constant in a cell
rMatrix Problem::basis_func_grad(const Cell &c, const Node &n)
{
    ElementArray<Node> nodes = c.getNodes();
	double x[3];
	double y[3];
	// gradient of the basis function
	rMatrix grad(2,1);
	unsigned n_ind = 0;
	for(unsigned i = 0; i < 3; i++){
		if(n == nodes[i])
			n_ind = i;
		double coords[3];
		nodes[i].Centroid(coords);
		x[i] = coords[0];
		y[i] = coords[1];
	}
	
	if(n_ind == 0){
		grad(0,0) = (y[1] - y[2]);
		grad(1,0) = - (x[1] - x[2]);
		grad /= ((x[0] - x[2])*(y[1] - y[2]) - (x[1] - x[2])*(y[0] - y[2]));
	}
	else if(n_ind == 1){
		// implement yourself
	}
	else if(n_ind == 2){
		// implement yourself
	}
	else{
		printf("Unexpected n_ind = %d\n", n_ind);
		exit(1);
	}
	return grad;
}

// Barycentric coordinates to normal ones
// For more info see page 39 in https://www.inm.ras.ru/wp-content/uploads/library/Monographies/yuv-kapyrin-svt-prak.pdf
void Problem::coords_from_barycentric(double* node_x, double* node_y, double* eta, double* x, double* y)
{
	*x = node_x[0] * eta[0] + node_x[1] * eta[1] + node_x[2] * eta[2];
	*y = node_y[0] * eta[0] + node_y[1] * eta[1] + node_y[2] * eta[2];
}

// Approximation by a linear combination of basis functions
// Coefficients are given by the nodal values of tag T
double Problem::approximate(const Cell &c, const Tag &T, double x, double y)
{
    ElementArray<Node> nodes = c.getNodes();
    double res = 0.0;
    for(unsigned i = 0; i < 3; i++){
        double xn[3];
        nodes[i].Barycenter(xn);
        res += nodes[i].Real(T) * basis_func(c, nodes[i], x, y);
    }
    return pow(C(x, y) - res, 2);
}

// Integration of a function over triangle
// The function is defined by its nodal values in tag T
// Problem::approximate is used to find values outside of nodes
// See pages 38-40 in https://www.inm.ras.ru/wp-content/uploads/library/Monographies/yuv-kapyrin-svt-prak.pdf
double Problem::integrate_over_triangle(const Cell &c, const Tag &T)
{
    double res = 0.0;
    double w3 = 0.205950504760887;
    double w6 = 0.063691414286223;
    double eta3[3] = {0.124949503233232, 0.437525248383384, 0.437525248383384};
    double eta6[3] = {0.797112651860071, 0.165409927389841, 0.037477420750088};

    ElementArray<Node> nodes = c.getNodes();
    if(nodes.size() != 3){
        printf("Cell is not a triangle, has %llu nodes!\n", nodes.size());
        exit(1);
    }
    // Coordinates of triangle nodes
    double node_x[3], node_y[3];
    // Set them
    for(unsigned i = 0; i < 3; i++){
        double c[3];
        nodes[i].Centroid(c);
        node_x[i] = c[0];
        node_y[i] = c[1];
    }

    // Add contribution from all combinations in eta3
    double x, y, val;
    double eta[3];
    eta[0] = eta3[0];
    eta[1] = eta3[1];
    eta[2] = eta3[2];
    coords_from_barycentric(node_x, node_y, eta, &x, &y);
    val = approximate(c, T, x, y);
    //printf("x = %e, y = %e, val = %e\n", x, y, val);
    res += w3 * val;
    eta[0] = eta3[1];
    eta[1] = eta3[2];
    eta[2] = eta3[0];
    coords_from_barycentric(node_x, node_y, eta, &x, &y);
    val = approximate(c, T, x, y);
    res += w3 * val;
    eta[0] = eta3[2];
    eta[1] = eta3[0];
    eta[2] = eta3[1];
    coords_from_barycentric(node_x, node_y, eta, &x, &y);
    val = approximate(c, T, x, y);
    res += w3 * val;


    // Add contribution from all combinations in eta6
    eta[0] = eta6[0];
    eta[1] = eta6[1];
    eta[2] = eta6[2];
    coords_from_barycentric(node_x, node_y, eta, &x, &y);
    val = approximate(c, T, x, y);
    res += w6 * val;
    eta[0] = eta6[0];
    eta[1] = eta6[2];
    eta[2] = eta6[1];
    coords_from_barycentric(node_x, node_y, eta, &x, &y);
    val = approximate(c, T, x, y);
    res += w6 * val;
    eta[0] = eta6[1];
    eta[1] = eta6[0];
    eta[2] = eta6[2];
    coords_from_barycentric(node_x, node_y, eta, &x, &y);
    val = approximate(c, T, x, y);
    res += w6 * val;
    eta[0] = eta6[1];
    eta[1] = eta6[2];
    eta[2] = eta6[0];
    coords_from_barycentric(node_x, node_y, eta, &x, &y);
    val = approximate(c, T, x, y);
    res += w6 * val;
    eta[0] = eta6[2];
    eta[1] = eta6[0];
    eta[2] = eta6[1];
    coords_from_barycentric(node_x, node_y, eta, &x, &y);
    val = approximate(c, T, x, y);
    res += w6 * val;
    eta[0] = eta6[2];
    eta[1] = eta6[1];
    eta[2] = eta6[0];
    coords_from_barycentric(node_x, node_y, eta, &x, &y);
    val = approximate(c, T, x, y);
    res += w6 * val;

    res *= c.Volume();
    return res;
}

// [ [phi1,phi1], [phi1,phi2]...   ]
// [   ]
// [   ]
// energetic scalar product: [u,v] = (D*grad(u), grad(v)) = (grad(v))^T * D * grad(u)
void Problem::assembleLocalSystem(const Cell &c, rMatrix &A_loc, rMatrix &rhs_loc)
{
     rMatrix D(2,2);
	 D(0,0) = c.RealArray(tagD)[0];
	 D(1,1) = c.RealArray(tagD)[1];
	 D(0,1) = D(1,0) = c.RealArray(tagD)[2];
	 
	 ElementArray<Node> nodes = c.getNodes();
	 
	 A_loc = rMatrix(3,3);
	 rhs_loc = rMatrix(3,1);
	 
	 for(unsigned i = 0; i < 3; i++){
         // fill matrix and rhs
	 }
	 
}

void Problem::assembleGlobalSystem(Sparse::Matrix &A, Sparse::Vector &rhs)
{
	// Cell loop
	// For each cell assemble local system
	// and incorporate it into global
	for(Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++){
		Cell c = icell->getAsCell();
		rMatrix A_loc, rhs_loc;
		assembleLocalSystem(c, A_loc, rhs_loc);
		// Now A_loc is 3x3, rhs_loc is 3x1
		//  
		//
		//
		//

		ElementArray<Node> nodes = c.getNodes();
		unsigned glob_ind[3] = {0,0,0};
		for(unsigned loc_ind = 0; loc_ind < 3; loc_ind++){
            // How to find global index?
            //glob_ind[loc_ind] = ???;
		}
		
		for(unsigned loc_ind = 0; loc_ind < 3; loc_ind++){
			// For node 'loc_ind' ADD (so use +=) its element in A_loc to the corresponding place in A
			
			// Check if this is a Dirichlet node
			if(nodes[loc_ind].GetMarker(mrkDirNode))
				continue;
			
			for(unsigned j = 0; j < 3; j++){
				if (nodes[j].GetMarker(mrkDirNode)) {
					///
				}
				else
					;//
				
			}
			// Don't forget the RHS
		}
	}
}

double Problem::get_err_c_norm() {
    double normC = 0.0;
    for(Mesh::iteratorNode inode = m.BeginNode(); inode != m.EndNode(); inode++){
        Node n = inode->getAsNode();
		// implement yourself
    }
    return normC;
}

double Problem::get_err_L2_norm() {
    double normL2 = 0.0;
    for(Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++){
        Cell c = icell->getAsCell();
        normL2 += integrate_over_triangle(c, tagConcErr);
    }
    normL2 = sqrt(normL2);
    return normL2;
}

void Problem::run()
{
	// Matrix size
	unsigned N = static_cast<unsigned>(m.NumberOfNodes()) - numDirNodes;
	// Global matrix called 'stiffness matrix'
	Sparse::Matrix A;
	// Solution vector
	Sparse::Vector sol;
	// Right-hand side vector
	Sparse::Vector rhs;

	A.SetInterval(0, N);
	sol.SetInterval(0, N);
	rhs.SetInterval(0, N);

	assembleGlobalSystem(A, rhs);

	A.Save("A.mtx");
	rhs.Save("rhs.mtx");

	string solver_name = "inner_mptiluc";
	Solver S(solver_name);

	S.SetMatrix(A);
	bool solved = S.Solve(rhs, sol);
	if(!solved){
		printf("Linear solver failed: %s\n", S.GetReason().c_str());
		printf("Number of iterations: %d\n", S.Iterations());
		printf("Residual:             %e\n", S.Residual());
		exit(1);
	}

	for(Mesh::iteratorNode inode = m.BeginNode(); inode != m.EndNode(); inode++){
		Node n = inode->getAsNode();
		if(n.GetMarker(mrkDirNode)){
			n.Real(tagConc) = n.Real(tagBCval);
			continue;
		}
		unsigned ind = static_cast<unsigned>(n.Integer(tagGlobInd));
		n.Real(tagConc) = sol[ind];
		// Implement yourself: setting error tag
	}
	m.Save("res.vtk");
}

int main(int argc, char ** argv)
{
	if( argc < 2 )
	{
		printf("Usage: %s mesh_file\n",argv[0]);
		return -1;
	}

	Mesh m;
	m.Load(argv[1]);
	Problem P(m);
	P.initProblem();
	P.run();

    cout << "|u - u_approx|_C = "  << P.get_err_c_norm() << endl;
    cout << "|u - u_approx|_L2 = " << P.get_err_L2_norm() << endl;
	printf("Success\n");
	return 0;
}
