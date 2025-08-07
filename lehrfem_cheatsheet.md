

```c
//create mesh pointer

std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space_p->Mesh();
//generating a mesh
const std::shared_ptr<lf::mesh::Mesh> mesh_ptr =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);


//from Triplet to Sparse
Eigen::SparseMatrixXd A; 
A.setFromTriplets(trip_.begin(), trip_.end()); 
Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Lower> solver; 



//getting the dofh from  finite element space
const lf::assemble::DofHandler &dofh{fe_space_p->LocGlobMap()};



// Obtain area of triangle
const double area = lf::geometry::Volume(*cell->Geometry());


//initilizing Matrix in COO format
lf::assemble::COOMatrix<double> M_COO(N_dofs, N_dofs);


//the mesh functions for the Reaction Diff

auto zero_mf = lf::mesh::utils::MeshFunctionConstant(0.0);
  auto one_mf = lf::mesh::utils::MeshFunctionConstant(1.0);
  
  
  
  //get Galerkin Matrices with alpha and gamma factor used
  
  lf::uscalfe::ReactionDiffusionElementMatrixProvider M_locmat_builder(
      fe_space_p, zero_mf, one_mf);
  lf::uscalfe::ReactionDiffusionElementMatrixProvider A_locmat_builder(
      fe_space_p, one_mf, zero_mf);
      
      
  //Galerkin Matrix for edge Mass matrix
  
  lf::mesh::utils::CodimMeshDataSet<bool> bd_flags{
      lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  lf::uscalfe::MassEdgeMatrixProvider B_locmat_builder(fe_space_p, one_mf,
                                                       bd_flags);
                                                       
                                                       
//assemble the Galerkin Matrix        
lf::assemble::AssembleMatrixLocally(0, dofh, dofh, M_locmat_builder, M_COO);


```
### integration on reference element
```c
//?????
//integration element, reference element, transformation function 


//lf::quad::QuadRule, lf::quad::make_QuadRule ! how to use 

//lf::geometryintegrationElement

//JacobianInverseGramian

//EvalReferenceShapeFunctions

//GradientsRegerenceShapeFUnctions

//ScalarReferenceFiniiteElement (EvalRegerenceShapeFunctions and GradientsRegerenceShapeFunctions

//lf::geometry::Geometry::Global




```
### general remarks
```c

//fixflaggedSolutionCompAlt

//lf::geometry::Geometry
//lf::geometry::Volume

//from local to global indices
Eigen::Vector3i dofhk = Mesh._elements.row(i) (where i is 0 till M (GSF)




// (Eigen::Vector2d() << 2, 2).finished()//tells eigen to actually make the vector// .replicate() //replicates in direction of vector vertically if transposed horizontally 


//get dofh from mesh
std::shared_ptr<fe::ScalarFESpace<double>> fe_space = std::make_shared<fe::FeSpaceLagrangeO1<double>>(mesh)); 
std::shared_ptr<assemble::DofHandler> dofh= fe_space->LocGlobMap()


//Fixing boundary for Dirichlet Data 

 // Obtain an array of boundary flags for edges (codim-1 entities !)
  lf::mesh::utils::CodimMeshDataSet<bool> boundary_edges{
      lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 1)};
 
  
  
 //fix the Entities on the boundary, with the Dirichlet data or a function g 
 //the structure is FixFlaggedSolutionComponents<double>(selector, A, phi)

auto selector_zero_bd = [&boundary_edges, &dof_handler](
          lf::assemble::glb_idx_t gdof_idx) -> std::pair<bool, double> {
        const lf::mesh::Entity &edge{dof_handler.Entity(gdof_idx)};
        return {boundary_edges(edge), 0.0};

auto selector_function_g_bd = [&boundary_edges, &dof_handler, &g](
          lf::assemble::glb_idx_t gdof_idx) -> std::pair<bool, double> {
        Eigen::VectorXd gt_vec = lf::fe::NodalProjection(*fe_space, g);
        const bool boundary = bd_flags(dofh.Entity(gdof_idx)); 
        if(boundary){
          return {true, gt_vec[gdof_idx]}
        }else{
          return {false, 0.0}; 
        }
        
        return {boundary_edges(edge), 0.0};


  lf::assemble::FixFlaggedSolutionComponents<double>(selector, A_COO, phi);
      
      
      
  // Set up Galerkin matrix in CRS format
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();
  // ... and solve the linear system of equations by Gaussian elimination
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  sol = solver.solve(phi);


//lf::mesh::utils::AllCodimMeshDataSet<bool> inflow_flag or boundary flags




 auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};
 
  // Assigning zero to the boundary values of phi
  for (const lf::mesh::Entity *vertex : mesh_p->Entities(2)) {
    if (bd_flags(*vertex)) {
      auto dof_idx = dofh.GlobalDofIndices(*vertex);
      LF_ASSERT_MSG(
          dofh.NumLocalDofs(*vertex) == 1,
          "Too many global indices were returned for a vertex entity!");
      phi(dof_idx[0]) = 0.0;
    }
  }
  return phi;



//setting certain position in A matrix
A.setZero([&selector](lf::assemble::gdof_idx_t i, lf::assemble::gdof_idx_t j){ return selector(i) || selecotr(j)};); 




//when splitting up Matrices like M = M_00 M_0BC
when use auto pred = [bd_flags, dofh](int i, int j){return bd_flags(dofh.Entity(i));};
M_00.setZero(pred)


//if uv is only valid on the boundary 
lf::uscalfe::MassEdgeMatrixProvider(fe_space_p,  one_mf, bd_flags); 



```
### computing error

```c

// either L2 error or H1 error

//L2 error, H1 error simply with GradFE
int quad_degree = 4; 

lf::fe::MeshFunctionGradeFE mf(fes_p, mu_vec);
lf::fe::MeshFunctionFE mf(fes_p, mu_vec) 


auto u1 = [&v_ex](Eigen::Vector2d x) -> double { return v_exact(x)[0]; };
const lf::mesh::utils::MeshFunctionGlobal mf_u1{u1}; 

lf::fe::IntegrateMeshFunction(*(fes_p->Mesh()), lf::mesh::utils::squaredNorm(mf - mf_u1), quad_degree); 


//computing error only over boundary term 
//for error over boundary terms, use the boundary selector
auto bd_flags = lf::mesh::utils::flagEntitiesOnBoundary(fes_p->Mesh(), 1);
auto bd_sel = [&bd_flags](const lf::mesh::Entity &edge) {
    return bd_flags(edge);
};
lf::fe::IntegrateMeshFunction(*(fes_p->Mesh()), lf::mesh::utils::squaredNorm(mf - mf_sol), 4, bd_sel, 1); 
//(multiplicative trace inequality for error rate and CSI) 





```
### Element Provider
```c

//LinearFElaplaceElementMatrix
//LinearFELocalLoadVector




auto fes_o2_ptr = std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(lev_mesh_p);


//DONT DO THIS THERE IS  MeshFunctionFE and IntegrateMeshFunctionFE
//we dont need to calculate the gradients manually
//get the H1 error (error of u exact and u_h) 
//some usefull code to get the H1 norm by using RefGradients and then pullback via JacInv.tranpose()

//define first reference coordinates
 Eigen::MatrixXd ref_coords(2, 6);
  // Example: midpoints of edges + barycenter + vertices
  ref_coords << 0.0, 1.0, 0.0, 0.5, 0.5, 1.0/3.0,
                0.0, 0.0, 1.0, 0.5, 0.5, 1.0/3.0;

 unsigned num_nodes = 6;
 std::vector<std::vector<Eigen::MatrixXd>> vec_physical_gradients;
 
 //iterating over cells (dofh vs. fes_p)
for(const lf::mesh::Entity& cell : 	fes_p->Entities(0)) 
OR
for(const lf::mesh::Entity* cell : dofh.Mesh()->Entities(0)) 


{
  auto geo = cell.Geometry(); // get geometry object
  Eigen::MatrixXd ref_coords = geo; // your N evaluation points
  auto jacobians = geo->Jacobian(ref_coords); // returns vector of 2x2 matrices

  auto ref_gradients = fe_o2.GradientsReferenceShapeFunctions();
  std::vector<Eigen::MatrixXd> physical_gradients(num_nodes); // each will be 2×6

  for (int i = 0; i < num_nodes; ++i) {
    Eigen::Matrix2d J_invT = jacobians[i].inverse().transpose();
    physical_gradients[i] = J_invT * ref_gradients[i]; // still 2×6
  }
  vec_physical_gradients.emplace_back(physical_gradients);


//Adding to COO Matrix
A.AddToEntry(n, tent_idx, area / 3.0)


//calculating gradients
G = J^{-T} * Gradients_on_reference_triangle

//TWO approaches
Eigen::Matrix<double, 2, 3> G =
  cell.Geometry()->JacobianInverseGramian(dpt).block(0, 0, 2, 2) *
  (Eigen::Matrix<double, 2, 3>(2,3) << -1, 1, 0,
                                       -1, 0, 1).finished(); 
//OR
auto endpoints = lf::geometry::Corners(*(cell.Geometry()));
  Eigen::Matrix<double, 3, 3> X;  // temporary matrix
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = endpoints.transpose();
  // This matrix contains $\cob{\grad \lambda_i}$ in its columns
  const auto G{X.inverse().block<2, 3>(1, 0)};

//Get subentities of a cell and iterate over  them 

    const std::span<const lf::mesh::Entity* const> nodes{cell->SubEntities(2)};
    // 0: cell itself, 1: segments, 2: points/nodes
    // Loop over nodes
    for (const lf::mesh::Entity* node : nodes)


Fixing for Dirichlet Data

//how to use bd_flags
auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fes_p->Mesh(),2); 
auto bdy_vertices_selector = [&bd_flags, &dofh](unsigned int idx) -> bool {
	return bd_flags(dofh.Entity(idx)));}; 

FixFlaggedSolutionComponents() or FixFlaggedSolutionCompAlt()

```
### Offset function and use of Selector
```c
//first (for problem a(.,.) =0;  U = u0+w => a(w, v) = - a(u0, v) <=> a(.,.) = -l(.)

//reacationDiffusionElementMatrixProvider
//AssembleMatrixLocally
//selector uses nodeflags(node) function
//Selector = [](lf::assemble::gdof_idx_t idx)->std::pair<bool, double>{}
//fixFlaggedSolutionComponents<double>(selectir, A, phi); 

```
### Non-linear elliptic BVP
```c

FunctionMFWrapper mf_coeff(mf_uh_prev, lambda_function); 
//simply represent the non-linear function/term
auto lamda_function = [](double xi)->double {return sinh(xi); }; 
FunctionMFWrapper mf_coeff_LHS(mf_grad_uh_prev, lambda_function);


```
### Stokes Problem
```c

auto edge_orientations = cell_p->RelativeOrientations(); 
lf::mesh::Orientation::positive
//internal odering depen on the orientation of the edge

//gradbarycoordinates(lf::mesh::Entity &entity){
 auto endpoints = lf::geometry::Corners(*(cell.Geometry()));
  Eigen::Matrix<double, 3, 3> X;  // temporary matrix
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = endpoints.transpose();
  // This matrix contains $\cob{\grad \lambda_i}$ in its columns
  const auto G{X.inverse().block<2, 3>(1, 0)};
  

// fixing boundary conditions (Dirichlet) (P2-P1)
auto bd_flags{//as usual}

//fixing for the nodes in x and y direction
 for (const lf::mesh::Entity *node : mesh_p->Entities(2)) {
    if (bd_flags(*node)) {
      // Indices of global shape functions sitting at node
      std::span<const lf::assemble::gdof_idx_t> dof_idx{
          dofh.InteriorGlobalDofIndices(*node)};
      LF_ASSERT_MSG(dof_idx.size() == 3, "Node must carry 3 dofs!");
      // Position of node
      const Eigen::Vector2d pos{Corners(*(node->Geometry())).col(0)};
      // Dirichlet data
      const Eigen::Vector2d g_val{g(pos)};
      // x-component of the velocity
      ess_dof_select[dof_idx[0]] = {true, g_val[0]};
      // y-component of the velocity
      ess_dof_select[dof_idx[1]] = {true, g_val[1]};
    }
  }
 for (const lf::mesh::Entity *edge : mesh_p->Entities(1)) {
    if (bd_flags(*edge)) {
      // Indices of global shape functions associated with the edge
      std::span<const lf::assemble::gdof_idx_t> dof_idx{
          dofh.InteriorGlobalDofIndices(*edge)};
      LF_ASSERT_MSG(dof_idx.size() == 2, "Edge must carry 2 dofs!");
      // Midpoint of edge
      const Eigen::MatrixXd endpoints{Corners(*(edge->Geometry()))};
      const Eigen::Vector2d pos{0.5 * (endpoints.col(0) + endpoints.col(1))};
      // Dirichlet data
      const Eigen::Vector2d g_val{g(pos)};
      // x-component of the velocity
      ess_dof_select[dof_idx[0]] = {true, g_val[0]};
      // y-component of the velocity
      ess_dof_select[dof_idx[1]] = {true, g_val[1]};
    }
  }
  
  //this is the core part, which help altering the boundary values
  // modify linear system of equations
  
  //required input
  std::vector<std::pair<bool, double>> ess_dof_select(n + 1, {false, 0.0});
  //we then adjust the the value in the vector based on some Dirichlet data  
  lf::assemble::FixFlaggedSolutionComponents<double>(
      [&ess_dof_select](lf::assemble::glb_idx_t dof_idx)
          -> std::pair<bool, double> { return ess_dof_select[dof_idx]; },
      A, phi);




//get convergence

   lf::assemble::UniformFEDofHandler dofh(lev_mesh_p,
                                           {{lf::base::RefEl::kPoint(), 3},
                                            {lf::base::RefEl::kSegment(), 2},
                                            {lf::base::RefEl::kTria(), 0},
                                            {lf::base::RefEl::kQuad(), 0}});

//dofh for U_x, U_y and p
 auto fes_o1_ptr =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(lev_mesh_p);
    auto fes_o2_ptr =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(lev_mesh_p);
    // Fetch dof handler for the components of the velocity
    const lf::assemble::DofHandler& dofh_u = fes_o2_ptr->LocGlobMap();
    //  Fetch dof handler for the pressure p
    const lf::assemble::DofHandler& dofh_p = fes_o1_ptr->LocGlobMap();

//solve problem
    Eigen::VectorXd res = StokesPipeFlow::solvePipeFlow(dofh, v_ex);
 
//calculate error for u_1 in L2 and H1, we assume coeff_vec_u1 already constains the solution mu
 const lf::fe::MeshFunctionFE mf_o2_u1(fes_o2_ptr, coeff_vec_u1);
 const lf::fe::MeshFunctionGradFE mf_o2_grad_u1(fes_o2_ptr, coeff_vec_u1);
//convert exact solution to mesh function
    auto u1 = [&v_ex](Eigen::Vector2d x) -> double { return v_ex(x)[0]; };
    const lf::mesh::utils::MeshFunctionGlobal mf_u1{u1};
//calculate error L2 and H1
L2err_u1 = std::sqrt(lf::fe::IntegrateMeshFunction(
        *lev_mesh_p, lf::mesh::utils::squaredNorm(mf_o2_u1 - mf_u1), 4));

H1err_u1 = std::sqrt(lf::fe::IntegrateMeshFunction(
        *lev_mesh_p, lf::mesh::utils::squaredNorm(mf_o2_grad_u1 - mf_grad_u1),
        4));



```

### FEEC

```c

//cell orientation
//indicates whereedge and cellorientation differ
auto cellOrientations = cell->RelativeOrientations(); 
//to be used like this for Indicence matrix cell-edge
auto ornt_It = cellOrientations.begin(); 
D(cell_idx, edge_idx) += lf::mesh::to_sign(*ornt_It); 


//Building EdgeVertexINdcidenceMatrix
lf::mesh::Mesh::Index(Entity & e) //index ranges from 0 till no. of entities of same co-dimension!
//example:
auto nodes = edge->SubEntities(1); //gives two nodes, start- and endpoint
lf::mesh::Mesh::size_type firstNodeIdx = mesh.Index(*nodes[0]) gives index of this node in comparison to all Nodes

//building CellEdgeIncidenceMatrix
//remember mesh.Entities(0) returns a pointer!!! to a cell object! (= lf::mesh::Entity *cell)
auto edges = cell->Subentities(1)
auto edgeOrientations = cell->RelativeOrientations(); 
lf::mesh::Mesh::size_type edge_idx = mesh.Index(*edges(i)); 
lf::mesh::Mesh::size_tyoe cell_idx = mesh.Index(*cell); 
D.coeffRef(cell_idx, edge_idx) += lf::mesh::to_sign(edgeOrientations[i])



