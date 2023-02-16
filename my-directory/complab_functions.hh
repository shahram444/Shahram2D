using namespace plb;
typedef double T;

#define NSDES descriptors::D2Q9Descriptor // Cs2 = 1/3
#define RXNDES descriptors::AdvectionDiffusionD2Q5Descriptor // Cs2 = 1/3
#define thrd DBL_EPSILON
// #define thrd 1e-14

class PressureGradient {
public:
    PressureGradient(T deltaP_, plint nx_) : deltaP(deltaP_), nx(nx_)
    { }
    void operator() (plint iX, plint iY, T& density, Array<T,2>& velocity) const
    {
        velocity.resetToZero();
        density = (T)1 - deltaP*NSDES<T>::invCs2 / (T)(nx-1) * (T)iX;
    }
private:
    T deltaP;
    plint nx;
};

void writeNsVTI(MultiBlockLattice2D<T,NSDES>& lattice, plint iter, std::string nameid)
{
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNx();
    VtkImageOutput2D<T> vtkOut(createFileName(nameid, iter, 7), 1.);

    vtkOut.writeData<float>(*computeVelocityNorm(lattice, Box2D(1,nx-2,0,ny-1)), "velocityNorm", 1.);
    vtkOut.writeData<2,float>(*computeVelocity(lattice, Box2D(1,nx-2,0,ny-1)), "velocity", 1.);
}

void writePorousMediumVTI(MultiScalarField2D<int>& geometry, plint iter, std::string nameid)
{
    const plint nx = geometry.getNx();
    const plint ny = geometry.getNx();
    VtkImageOutput2D<T> vtkOut(createFileName(nameid, iter, 7), 1.);
    vtkOut.writeData<float>(*copyConvert<int,T>(geometry, Box2D(1,nx-2,0,ny-1)), "tag", 1.0);
}

void writeAdvVTI(MultiBlockLattice2D<T,RXNDES>& lattice, plint iter, std::string nameid)
{
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNx();

    VtkImageOutput2D<T> vtkOut(createFileName(nameid, iter, 7), 1.);
    // vtkOut.writeData<T>(*computeDensity(lattice, Box2D(0,nx-1,0,ny-1)), "Density", 1.);
    vtkOut.writeData<T>(*computeDensity(lattice, Box2D(1,nx-2,0,ny-1)), "Density", 1.);
}

void writeScalarVTI(MultiScalarField2D<int>& field)
{
    const plint nx = field.getNx();
    const plint ny = field.getNx();

    VtkImageOutput2D<T> vtkOut("distanceDomain", 1.);
    vtkOut.writeData<float>(*copyConvert<int,T>(field, Box2D (1,nx-2, 0,ny-1)), "tag", 1.0);
}

void readGeometry(std::string fNameIn, MultiScalarField2D<int>& geometry)
{
    pcout << "Reading the geometry file (" << fNameIn << ").\n";
    const plint nx = geometry.getNx();
    const plint ny = geometry.getNy();

    Box2D sliceBox(0,0, 0,ny-1);
    std::unique_ptr< MultiScalarField2D<int> > slice = generateMultiScalarField<int>(geometry, sliceBox);

    plb_ifstream geometryFile(fNameIn.c_str());
    if (!geometryFile.is_open()) {
        pcout << "Error: could not open geometry file " << fNameIn << std::endl;
        exit(EXIT_FAILURE);
    }
    for (plint iX=1; iX<nx-1; ++iX) {
        geometryFile >> *slice;
        if (iX==1) {
            copy(*slice, slice->getBoundingBox(), geometry, Box2D(0,0, 0,ny-1));
            copy(*slice, slice->getBoundingBox(), geometry, Box2D(iX,iX, 0,ny-1));
        }
        else if (iX==nx-2) {
            copy(*slice, slice->getBoundingBox(), geometry, Box2D(iX,iX, 0,ny-1));
            copy(*slice, slice->getBoundingBox(), geometry, Box2D(nx-1,nx-1, 0,ny-1));
        }
        else {
            copy(*slice, slice->getBoundingBox(), geometry, Box2D(iX,iX, 0,ny-1));
        }
    }
}

void saveGeometry(std::string fNameIn, MultiScalarField2D<int>& geometry)
{
    const plint nx = geometry.getNx();
    const plint ny = geometry.getNy();

    pcout << "Save geometry vti file ("<< fNameIn <<").\n";
    VtkImageOutput2D<T> vtkOut(fNameIn, 1.0);
    vtkOut.writeData<float>(*copyConvert<int,T>(geometry, Box2D (1,nx-2, 0,ny-1)), "tag", 1.0);
}

void calculateDistanceFromSolid(MultiScalarField2D<int> distance, plint nodymcs, plint bb, std::vector< std::vector<plint> > &distVec)
{
    const plint nx = distance.getNx();
    const plint ny = distance.getNy();

    for (plint iX=0; iX<nx-1; ++iX) {
        for (plint iY=0; iY<ny-1; ++iY) {
            plint mask = distance.get(iX,iY);
            if (mask == nodymcs) { distVec[iX][iY] = -1; }
            else if (mask == bb) { distVec[iX][iY] = 0; }
            else { distVec[iX][iY] = 1; }
        }
    }
    for (plint iX=0; iX<nx-1; ++iX) {
        for (plint iY=0; iY<ny-1; ++iY) {
            if ( distVec[iX][iY]==1 ) {
                plint lp=1, r=0, dist=0;
                while ( lp==1 ) {
                    ++r;
                    std::vector<plint> vx(r+1), vy(r+1);
                    for (plint tmp=0; tmp<r+1; ++tmp) { vx[tmp]=tmp; vy[tmp]=r-tmp; }
                    for (plint it=0; it<r+1; ++it) {
                        plint xp = iX+vx[it], yp = iY+vy[it], xn = iX-vx[it], yn = iY-vy[it];
                        if (xp>=0 && yp>=0 && xp<nx && yp<ny) {
                            if ( distVec[xp][yp]==0 ) {
                                dist = r; lp = 0; break;
                            }
                        }
                        if (xp>=0 && yn>=0 && xp<nx && yn<ny) {
                            if ( distVec[xp][yn]==0 ) {
                                dist = r; lp = 0; break;
                            }
                        }
                        if (xn>=0 && yp>=0 && xn<nx && yp<ny) {
                            if ( distVec[xn][yp]==0 ) {
                                dist = r; lp = 0; break;
                            }
                        }
                        if (xn>=0 && yn>=0 && xn<nx && yn<ny) {
                            if ( distVec[xn][yn]==0 ) {
                                dist = r; lp=0; break;
                            }
                        }
                    }
                }
                if (lp==0) { distVec[iX][iY] = dist; }
            }
        }
    }
}

void NSdomainSetup(MultiBlockLattice2D<T,NSDES>& lattice, OnLatticeBoundaryCondition2D<T,NSDES>* boundaryCondition, MultiScalarField2D<int>& geometry, T deltaP,
                    T fluidOmega, std::vector<plint> pore, plint bounceback, plint nodymcs, std::vector<std::vector<plint>> bio_dynamics, std::vector<T> permRatio)
{
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();

    Box2D west   (0,0, 0,ny-1);
    Box2D east   (nx-1,nx-1, 0,ny-1);

    // default. initialize the entire domain. may be redundant
    defineDynamics(lattice, lattice.getBoundingBox(), new IncBGKdynamics<T,NSDES>(fluidOmega));

    // pore space
    for (size_t iP=0; iP<pore.size(); ++iP) {
        if (pore[iP] > 0) { defineDynamics(lattice, geometry, new IncBGKdynamics<T,NSDES>(fluidOmega), pore[iP]); }
    }
    // bounce-back boundary
    if (bounceback > 0) {
        defineDynamics(lattice, geometry, new BounceBack<T,NSDES>(), bounceback);
    }
    // no dynamics
    if (nodymcs >= 0) {
        defineDynamics(lattice, geometry, new NoDynamics<T,NSDES>(),nodymcs);
    }
    // microbial material number
    for (size_t iM = 0; iM < bio_dynamics.size(); ++iM) {
        T bioOmega = 1/(permRatio[iM]*(1/fluidOmega-.5)+.5);
        for (size_t iB = 0; iB < bio_dynamics[iM].size(); ++iB) {
            if (bio_dynamics[iM][iB] > 0) {
                if (permRatio[iM] > thrd) { // permeable biofilm
                    defineDynamics(lattice, geometry, new IncBGKdynamics<T,NSDES>(bioOmega), bio_dynamics[iM][iB]);
                }
                else { // impermeable biofilm
                    defineDynamics(lattice, geometry, new BounceBack<T,NSDES>(), bio_dynamics[iM][iB]);
                }
            }
        }
    }

    boundaryCondition->addPressureBoundary0N(west, lattice);
    setBoundaryDensity(lattice, west, (T) 1.);
    boundaryCondition->addPressureBoundary0P(east, lattice);
    setBoundaryDensity(lattice, east, (T) 1. - deltaP*NSDES<T>::invCs2);

    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), PressureGradient(deltaP, nx));

    lattice.initialize();
    delete boundaryCondition;
}

void soluteDomainSetup( MultiBlockLattice2D<T,RXNDES> &lattice, OnLatticeAdvectionDiffusionBoundaryCondition2D<T, RXNDES>* boundaryCondition, MultiScalarField2D<int>& geometry,
                        T substr_bMassOmega, T substrOmega, std::vector<plint> pore, plint bounceback, plint nodymcs, std::vector<std::vector<plint>> bio_dynamics,
                        T rho0, bool left_btype, bool right_btype, T left_BC, T right_BC )
{
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();

    Box2D west   (0,0, 0,ny-1);
    Box2D east   (nx-1,nx-1, 0,ny-1);
    plint processorLevelBC = 1;

    // default. initialize the entire domain. may be redundant
    defineDynamics(lattice, lattice.getBoundingBox(), new AdvectionDiffusionBGKdynamics<T,RXNDES>(substrOmega));

    // pore space
    for (size_t iP=0; iP<pore.size(); ++iP) {
        if (pore[iP] > 0) { defineDynamics(lattice, geometry, new AdvectionDiffusionBGKdynamics<T,RXNDES>(substrOmega), pore[iP]); }
    }
    // bounceback boundary
    if (bounceback > 0) { defineDynamics(lattice, geometry, new BounceBack<T,RXNDES>(), bounceback); }
    // no dynamics
    if (nodymcs >= 0) { defineDynamics(lattice, geometry, new NoDynamics<T,RXNDES>(),nodymcs); }
    // microbial material number
    for (size_t iM = 0; iM < bio_dynamics.size(); ++iM) {
        for (size_t iB = 0; iB < bio_dynamics[iM].size(); ++iB) {
            if (bio_dynamics[iM][iB] > 0) { defineDynamics(lattice, geometry, new AdvectionDiffusionBGKdynamics<T,RXNDES>(substr_bMassOmega), bio_dynamics[iM][iB]); }
        }
    }

    // Set the boundary-conditions
    boundaryCondition->addTemperatureBoundary0N(west, lattice);
    if (left_btype == 0) { setBoundaryDensity(lattice, west, left_BC); }
    else { integrateProcessingFunctional(new FlatAdiabaticBoundaryFunctional2D<T,RXNDES, 0, -1>, west, lattice, processorLevelBC); }
    boundaryCondition->addTemperatureBoundary0P(east, lattice);
    if (right_btype == 0) { setBoundaryDensity(lattice, east, right_BC); }
    else { integrateProcessingFunctional(new FlatAdiabaticBoundaryFunctional2D<T,RXNDES, 0, +1>, east, lattice, processorLevelBC); }

    // Init lattice
    Array<T,2> u0(0., 0.);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), rho0, u0);

    lattice.initialize();
    delete boundaryCondition;
}

void bmassDomainSetup(MultiBlockLattice2D<T,RXNDES> &lattice, OnLatticeAdvectionDiffusionBoundaryCondition2D<T, RXNDES>* boundaryCondition,
                    MultiScalarField2D<int>& geometry, T bioOmegaPore, T bioOmegaFilm, std::vector<plint> pore, plint bounceback, plint nodymcs, std::vector<std::vector<plint>> bio_dynamics,
                    bool left_btype, bool right_btype, T left_BC, T right_BC )
{
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();
    plint processorLevelBC = 1;

    Box2D west   (0,0, 0,ny-1);
    Box2D east   (nx-1,nx-1, 0,ny-1);

    // default. initialize the entire domain. may be redundant
    defineDynamics(lattice, lattice.getBoundingBox(), new AdvectionDiffusionBGKdynamics<T,RXNDES>(bioOmegaPore));

    // pore space
    for (size_t iP=0; iP<pore.size(); ++iP) {
        if (pore[iP] > 0) { defineDynamics(lattice, geometry, new AdvectionDiffusionBGKdynamics<T,RXNDES>(bioOmegaPore), pore[iP]); }
    }
    // bounceback boundary
    if (bounceback > 0) { defineDynamics(lattice, geometry, new BounceBack<T,RXNDES>(), bounceback); }
    // no dynamics
    if (nodymcs >= 0) { defineDynamics(lattice, geometry, new NoDynamics<T,RXNDES>(),nodymcs); }
    // microbial material number
    for (size_t iM = 0; iM < bio_dynamics.size(); ++iM) {
        for (size_t iB = 0; iB < bio_dynamics[iM].size(); ++iB) {
            if (bio_dynamics[iM][iB] > 0) { defineDynamics(lattice, geometry, new AdvectionDiffusionBGKdynamics<T,RXNDES>(bioOmegaFilm), bio_dynamics[iM][iB]); }
        }
    }

    // Set the boundary-conditions
    boundaryCondition->addTemperatureBoundary0N(west, lattice);
    if (left_btype == 0) { setBoundaryDensity(lattice, west, left_BC); }
    else { integrateProcessingFunctional(new FlatAdiabaticBoundaryFunctional2D<T,RXNDES, 0, -1>, west, lattice, processorLevelBC); }
    boundaryCondition->addTemperatureBoundary0P(east, lattice);
    if (right_btype == 0) { setBoundaryDensity(lattice, east, right_BC); }
    else { integrateProcessingFunctional(new FlatAdiabaticBoundaryFunctional2D<T,RXNDES, 0, +1>, east, lattice, processorLevelBC); }

    // Init lattice
    Array<T,2> u0(0., 0.);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), 0., u0);

    lattice.initialize();
    delete boundaryCondition;
}

void scalarDomainDynamicsSetupFromVectors(MultiBlockLattice2D<T,RXNDES> &lattice, MultiScalarField2D<int>& geometry, std::vector<plint> mtrvec, std::vector<T> omegavec)
{
    // default. initialize the entire domain. may be redundant
    defineDynamics(lattice, lattice.getBoundingBox(), new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.));

    if ( mtrvec.size() != omegavec.size() ) {
        pcout << "ERROR: the length of input vectors (mtrvec and omegavec) must be the same.\n";
        exit(EXIT_FAILURE);
    }
    // assign lattice omegas (dynamics) for each mask number
    for (size_t iT=0; iT<mtrvec.size(); ++iT) {
        defineDynamics(lattice, geometry, new AdvectionDiffusionBGKdynamics<T,RXNDES>(omegavec[iT]), mtrvec[iT]);
    }
    // Init lattice
    Array<T,2> jEq(0., 0.);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), 0., jEq);

    lattice.initialize();
}

void scalarDomainDynamicsSetupFromGeometry(MultiBlockLattice2D<T,RXNDES> &lattice, MultiScalarField2D<int>& geometry, plint nx, plint ny)
{
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            plint geom = geometry.get(iX,iY);
            defineDynamics( lattice, iX, iY, new AdvectionDiffusionBGKdynamics<T,RXNDES>((T) geom) );
        }
    }
    // Init lattice
    Array<T,2> jEq(0., 0.);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), 0., jEq);

    lattice.initialize();
}

void gridSearch(MultiScalarField2D<int> geometry, std::vector< std::vector<plint> > &distVec, plint bb, plint solid, std::vector< std::vector<plint> > bio, std::vector<plint> pore)
{
    const plint nx = geometry.getNx();
    const plint ny = geometry.getNy();
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            bool flag0 = 0;
            plint geom = geometry.get(iX,iY);
            for (size_t iB0=0; iB0<bio.size(); ++iB0) {
                for (size_t iB1=0; iB1<bio[iB0].size(); ++iB1) {
                    if (geom == bio[iB0][iB1]) { flag0 = 1; }
                }
            }
            if ( flag0 == 1 ) {
                plint iR = 0; bool flag1 = 0;
                while (flag1 == 0) {
                    ++iR;
                    for (plint rx=0; rx<iR; ++rx) {
                        plint ry=iR-rx;
                        if (iX+rx<nx && iY+ry<ny) {
                            plint mask = geometry.get(iX+rx,iY+ry);
                            for (size_t iP=0; iP<pore.size(); ++iP) {
                                if (mask == pore[iP]) {
                                    flag1 = 1; distVec[iX][iY]=iR; break;
                                }
                            }
                            if (flag1==1) { break; }
                        }
                        if (iX+rx<nx && iY-ry>0) {
                            plint mask = geometry.get(iX+rx,iY-ry);
                            for (size_t iP=0; iP<pore.size(); ++iP) {
                                if (mask == pore[iP]) {
                                    flag1 = 1; distVec[iX][iY]=iR; break;
                                }
                            }
                            if (flag1==1) { break; }
                        }
                        if (iX-rx>0 && iY+ry<ny) {
                            plint mask = geometry.get(iX-rx,iY+ry);
                            for (size_t iP=0; iP<pore.size(); ++iP) {
                                if (mask == pore[iP]) {
                                    flag1 = 1; distVec[iX][iY]=iR; break;
                                }
                            }
                            if (flag1==1) { break; }
                        }
                        if (iX-rx<0 && iY-ry>0) {
                            plint mask = geometry.get(iX-rx,iY-ry);
                            for (size_t iP=0; iP<pore.size(); ++iP) {
                                if (mask == pore[iP]) {
                                    flag1 = 1; distVec[iX][iY]=iR; break;
                                }
                            }
                            if (flag1==1) { break; }
                        }
                    }
                    // if (iR == 2) {
                    //     pcout << "iX = " << iX << ", iY = " << iY << ", iX+1 = " << geometry.get(iX+1,iY) << ", iX-1 = " << geometry.get(iX-1,iY) << ", iY+1 = " << geometry.get(iX,iY+1) << ", iY-1 = " << geometry.get(iX,iY-1) << std::endl;
                    // }
                }
            }
            else if (geom == bb || geom == solid) { distVec[iX][iY]=-1; }
            else { distVec[iX][iY]=0; }
        }
    }
}

void initTotalbFilmLatticeDensity(MultiBlockLattice2D<T,RXNDES> &lattice1, MultiBlockLattice2D<T,RXNDES> &lattice2)
{
    const plint nx = lattice1.getNx();
    const plint ny = lattice1.getNy();
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            T bmass = lattice1.get(iX,iY).computeDensity();
            Array<T,5> g;
            lattice2.get(iX,iY).getPopulations(g);
            g[0]+=(T) (bmass)/3; g[1]+=(T) (bmass)/6; g[2]+=(T) (bmass)/6; g[3]+=(T) (bmass)/6; g[4]+=(T) (bmass)/6;
            lattice2.get(iX,iY).setPopulations(g);
        }
    }
}

void defineMaskLatticeDynamics(MultiBlockLattice2D<T,RXNDES> &lattice1, MultiBlockLattice2D<T,RXNDES> &lattice2, T fbM)
{
    const plint nx = lattice1.getNx();
    const plint ny = lattice1.getNy();
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            T bmass = lattice1.get(iX,iY).computeDensity();
            T omega = 0.;
            if ( bmass > fbM ) { omega = 1.; }
            defineDynamics( lattice2, iX, iY, new AdvectionDiffusionBGKdynamics<T,RXNDES>(omega) );
        }
    }
    // Init lattice
    Array<T,2> jEq(0., 0.);
    initializeAtEquilibrium(lattice2, lattice2.getBoundingBox(), 0., jEq);

    lattice2.initialize();
}

int initialize_complab( char *&main_path, char *&input_path, char *&output_path, char *&ns_filename, std::string &ade_filename, std::string &bio_filename, std::string &geom_filename,
    std::string &mask_filename, bool &read_NS_file, plint &ns_rerun_iT0, T &ns_converge_iT1, T &ns_converge_iT2, plint &ns_maxiTer_1, plint &ns_maxiTer_2, plint &ns_update_interval, plint &ade_update_interval,
    bool &read_ADE_file, plint &ade_rerun_iT0, plint &ade_VTK_iTer, plint &ade_CHK_iTer, T &ade_converge_iT, plint &ade_maxiTer, plint &nx, plint &ny, T &dx, T &dy, T &delta_P, T &tau,
    T &Pe, T &charcs_length, std::vector<T> &solute_poreD, std::vector<T> &solute_bFilmD, std::vector<T> &bMass_poreD, std::vector<T> &bMass_bFilmD, bool &soluteDindex, bool &bmassDindex,
    T &thrd_bFilmFrac, std::vector<T> &vec_permRatio, T &max_bMassRho, std::vector<plint> &pore_dynamics, plint &bounce_back, plint &no_dynamics, std::vector< std::vector<plint> > &bio_dynamics,
    plint &num_of_microbes, plint &num_of_substrates, std::vector<std::string> &vec_subs_names, std::vector<std::string>& vec_microbes_names, std::vector<plint> &solver_type, plint &fd_count,
    plint &lb_count, plint &ca_count, plint &bfilm_count, plint &bfree_count, std::vector< std::string > &vec_mmFileName, plint &mm_count, plint &kns_count, plint &nn_count, std::vector<plint> &reaction_type,
    std::vector< std::vector<T> > &vec_maxUptake, std::vector< std::vector<T> > &vec_maxRelease, std::vector< std::vector<int> > &vec_EX_loc,
    std::vector<T> &vec_c0, std::vector<bool> &left_btype, std::vector<bool> &right_btype, std::vector<T> &vec_leftBC, std::vector<T> &vec_rightBC, std::vector< std::vector<T> > &vec_b0_all,
    std::vector<bool> &bio_left_btype, std::vector<bool> &bio_right_btype, std::vector<T> &bio_leftBC, std::vector<T> &bio_rightBC,
    std::vector< std::vector<T> > &vec_Kc, std::vector< std::vector<T> > &vec_Kc_glpk, std::vector< std::vector<T> > &vec_Kc_kns, std::vector< std::vector<T> > &vec_Kc_nn,
    std::vector<T> &vec_mu, std::vector<T> &vec_mu_glpk, std::vector<T> &vec_mu_kns, std::vector<T> &vec_mu_nn, std::vector<bool> &vec_fixLB, std::vector<bool> &vec_fixC,
    std::vector<bool> &bmass_type, std::vector<T> &vec_b0_free, std::vector< std::vector<T> > &vec_b0_film, std::vector<plint> &vec_sense,
    std::vector<std::vector<T>> &vec_Vmax, std::vector<std::vector<T>> &vec_Vmax_glpk, std::vector<std::vector<T>> &vec_Vmax_kns, std::vector<std::vector<T>> &vec_Vmax_nn,
    std::vector<std::vector<int>> &vec_const_loc, std::vector<std::vector<T>> &vec_const_lb, std::vector<std::vector<T>> &vec_const_ub, bool &track_performance, bool &halfflag)

{
    try {
        std::string fin("CompLaB.xml");
        XMLreader doc(fin);

        // terminate the simulation if inputs are undefined.
        try {
            // LB_numerics
            doc["parameters"]["LB_numerics"]["domain"]["nx"].read(nx); nx+=2;
            doc["parameters"]["LB_numerics"]["domain"]["ny"].read(ny);
            doc["parameters"]["LB_numerics"]["domain"]["dx"].read(dx);
            doc["parameters"]["LB_numerics"]["domain"]["filename"].read(geom_filename);

            // chemistry
            doc["parameters"]["chemistry"]["number_of_substrates"].read(num_of_substrates);
            for (plint iT=0; iT<num_of_substrates; ++iT) {
                T bc0, bc1;
                std::string tmp0, tmp1;
                std::string chemname = "substrate"+std::to_string(iT);
                doc["parameters"]["chemistry"][chemname]["left_boundary_type"].read(tmp0);
                std::transform(tmp0.begin(), tmp0.end(), tmp0.begin(), [](unsigned char c){ return std::tolower(c); });
                if (tmp0.compare("dirichlet")==0) { left_btype.push_back(0); }
                else if (tmp0.compare("neumann")==0) { left_btype.push_back(1); }
                else { pcout << "left_boundary_type (" << tmp0 << ") should be either Dirichlet or Neumann. Terminating the simulation.\n"; return -1; }
                doc["parameters"]["chemistry"][chemname]["right_boundary_type"].read(tmp1);
                std::transform(tmp1.begin(), tmp1.end(), tmp1.begin(), [](unsigned char c){ return std::tolower(c); });
                if (tmp1.compare("dirichlet")==0) { right_btype.push_back(0); }
                else if (tmp1.compare("neumann")==0) { right_btype.push_back(1); }
                else { pcout << "right_boundary_type (" << tmp1 << ") should be either Dirichlet or Neumann. Terminating the simulation.\n"; return -1; }
                doc["parameters"]["chemistry"][chemname]["left_boundary_condition"].read(bc0); vec_leftBC.push_back(bc0);
                doc["parameters"]["chemistry"][chemname]["right_boundary_condition"].read(bc1); vec_rightBC.push_back(bc1);
            }
            if (left_btype.size()!=(unsigned)num_of_substrates) { pcout << "The length of left_boundary_type vector does not match the num_of_substrates. Terminating the simulation.\n"; return -1;}
            if (right_btype.size()!=(unsigned)num_of_substrates) { pcout << "The length of right_boundary_type vector does not match the num_of_substrates. Terminating the simulation.\n"; return -1;}
            if (vec_leftBC.size()!=(unsigned)num_of_substrates) { pcout << "The length of left_boundary_condition vector does not match the num_of_substrates. Terminating the simulation.\n"; return -1;}
            if (vec_rightBC.size()!=(unsigned)num_of_substrates) { pcout << "The length of right_boundary_condition vector does not match the num_of_substrates. Terminating the simulation.\n"; return -1;}

            // microbiology
            doc["parameters"]["microbiology"]["number_of_microbes"].read(num_of_microbes);
            mm_count = 0, kns_count = 0, nn_count = 0, fd_count = 0, ca_count = 0, lb_count = 0;
            for (plint iT=0; iT<num_of_microbes; ++iT) {
                std::string tmp0, tmp1;
                std::vector<T> b0;
                std::string bioname = "microbe"+std::to_string(iT);
                doc["parameters"]["microbiology"][bioname]["reaction_type"].read(tmp0);
                std::transform(tmp0.begin(), tmp0.end(), tmp0.begin(), [](unsigned char c){ return std::tolower(c); });
                if (tmp0.compare("glpk")==0) { reaction_type.push_back(0); ++mm_count; }
                else if (tmp0.compare("kinetics")==0) { reaction_type.push_back(1); ++kns_count; }
                else if (tmp0.compare("nn")==0 || tmp0.compare("neural_network")==0 || tmp0.compare("neural network")==0) { reaction_type.push_back(2); ++nn_count; }
                else if (tmp0.compare("glpk_and_kinetics")==0) { reaction_type.push_back(3); ++mm_count; ++kns_count; }
                else if (tmp0.compare("nn_and_kinetics")==0) { reaction_type.push_back(4); ++nn_count; ++kns_count; }
                else { pcout << "reaction_type " << tmp0 << " is not implemented. The default solver (glpk) is used.\n"; reaction_type.push_back(0); ++mm_count; }
                doc["parameters"]["microbiology"][bioname]["solver_type"].read(tmp1);
                std::transform(tmp1.begin(), tmp1.end(), tmp1.begin(), [](unsigned char c){ return std::tolower(c); });
                if (tmp1.compare("fd")==0 || tmp1.compare("finite difference")==0 || tmp1.compare("finite_difference")==0) { solver_type.push_back(1); ++fd_count; }
                else if (tmp1.compare("ca")==0 || tmp1.compare("cellular automata")==0 || tmp1.compare("cellular_automata")==0) { solver_type.push_back(2); ++ca_count; }
                else if (tmp1.compare("lbm")==0 || tmp1.compare("lattice boltzmann")==0 || tmp1.compare("lattice_boltzmann")==0) { solver_type.push_back(3); ++lb_count; }
                else {pcout << "Palabos IO exception: Element solver_type " << tmp1 << " is not defined. Use either FD, CA, or LBM. Terminating the simulation.\n"; return -1;}
                doc["parameters"]["microbiology"][bioname]["initial_densities"].read(b0); vec_b0_all.push_back(b0);
            }
            if (reaction_type.size() != (unsigned)num_of_microbes) { pcout << "The length of reaction_type vector does not match the num_of_microbes. Terminating the simulation.\n"; return -1; }
            if (solver_type.size() != (unsigned)num_of_microbes) { pcout << "The length of solver_type vector does not match the num_of_microbes. Terminating the simulation.\n"; return -1; }
            if (vec_b0_all.size() != (unsigned)num_of_microbes) { pcout << "The length of initial_densities vector does not match the num_of_microbes. Terminating the simulation.\n"; return -1; }
            if (mm_count > 0) {
                for (plint iT=0; iT<num_of_microbes; ++iT) {
                    if (reaction_type[iT]==0 || reaction_type[iT]==3) {
                        std::string tmp0, tmp4;
                        std::vector<int> tmp1;
                        std::vector<T> tmp2, tmp3;
                        std::string bioname = "microbe"+std::to_string(iT);
                        doc["parameters"]["microbiology"][bioname]["model_filename"].read(tmp0); vec_mmFileName.push_back(tmp0);
                        doc["parameters"]["microbiology"][bioname]["exchange_reaction_indices"].read(tmp1); vec_EX_loc.push_back(tmp1);
                        doc["parameters"]["microbiology"][bioname]["substrate_lower_bounds"].read(tmp2); vec_maxUptake.push_back(tmp2);
                        doc["parameters"]["microbiology"][bioname]["substrate_upper_bounds"].read(tmp3); vec_maxRelease.push_back(tmp3);
                        doc["parameters"]["microbiology"][bioname]["objective_direction"].read(tmp4);
                        std::transform(tmp4.begin(), tmp4.end(), tmp4.begin(), [](unsigned char c){ return std::tolower(c); });
                        if (tmp4.compare("max")==0 || tmp4.compare("maximize")==0 || tmp4.compare("-1")==0 ) { vec_sense.push_back(-1); }
                        else if (tmp4.compare("min")==0 || tmp4.compare("minimize")==0 || tmp4.compare("1")==0) { vec_sense.push_back(1); }
                        else { pcout << "objective_direction (" << tmp4 << ") should be either max or min. Terminating the simulation.\n"; return -1; }
                        if (tmp1.size()!=tmp2.size()) { pcout << "The length of exchange_reaction_indices vector does not match the length of substrate_lower_bounds vector for microbe" << iT
                            << "If one of substrates is not necessary for constraining metabolic model, use any negative number in substrate_lower_bounds vector. Terminating the simulation.\n"; return -1; }
                        if (tmp1.size()!=tmp3.size()) { pcout << "The length of exchange_reaction_indices vector does not match the length of substrate_upper_bounds vector for microbe" << iT
                            << "If one of substrates is not necessary for constraining metabolic model, use any negative number in substrate_upper_bounds vector. Terminating the simulation.\n"; return -1; }
                        if (tmp2.size()!=tmp3.size()) { pcout << "The length of substrate_lower_bounds vector does not match the length of substrate_upper_bounds vector for microbe" << iT <<". Terminating the simulation.\n"; return -1; }
                    }
                }
                if (vec_mmFileName.size()!= (unsigned)mm_count) { pcout << "The length of model_filename vector does not match the number of metabolic models used. Terminating the simulation.\n"; return -1; }
                if (vec_EX_loc.size() != (unsigned)mm_count) { pcout << "The length of exchange_reaction_indices vector does not match the number of metabolic models used. Terminating the simulation.\n"; return -1; }
                if (vec_maxUptake.size() != (unsigned)mm_count) { pcout << "The length of substrate_lower_bounds vector does not match the number of metabolic models used. Terminating the simulation.\n"; return -1; }
                if (vec_maxRelease.size() != (unsigned)mm_count) { pcout << "The length of substrate_upper_bounds vector does not match the number of metabolic models used. Terminating the simulation.\n"; return -1; }
            }
        }
        catch (PlbIOException& exception) {
            pcout << exception.what() << " Terminating the simulation.\n" << std::endl;
            return -1;
        }

        // parameters with default values
        // define paths
        try {
            std::string item;
            doc["parameters"]["path"]["input_path"].read(item);
            input_path = (char *) calloc(item.size()+1,sizeof(char));
            for (size_t i = 0; i < item.size(); ++i) { input_path[i] = item[i]; }
            input_path[item.size()+1] = '\0';
        }
        catch (PlbIOException& exception) {
            std::string item = "input";
            input_path = (char *) calloc(item.size()+1,sizeof(char));
            for (size_t i = 0; i < item.size(); ++i) { input_path[i] = item[i]; }
            input_path[item.size()+1] = '\0';
        }
        try {
            std::string item;
            doc["parameters"]["path"]["output_path"].read(item);
            output_path = (char *) calloc(item.size()+1,sizeof(char));
            for (size_t i = 0; i < item.size(); ++i) { output_path[i] = item[i]; }
            output_path[item.size()+1] = '\0';
        }
        catch (PlbIOException& exception) {
            std::string item = "output";
            output_path = (char *) calloc(item.size()+1,sizeof(char));
            for (size_t i = 0; i < item.size(); ++i) { output_path[i] = item[i]; }
            output_path[item.size()+1] = '\0';
        }

        // LB_numerics
        try { doc["parameters"]["LB_numerics"]["delta_P"].read(delta_P); }
        catch (PlbIOException& exception) { delta_P=0; }
        try {
            std::string tmp;
            doc["parameters"]["LB_numerics"]["track_performance"].read(tmp);
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](unsigned char c){ return std::tolower(c); });
            if (tmp.compare("no")==0 || tmp.compare("false")==0 || tmp.compare("0")==0 ) { track_performance = 0; }
            else if (tmp.compare("yes")==0 || tmp.compare("true")==0 || tmp.compare("1")==0 ) { track_performance = 1; }
            else { pcout << "track_performance (" << tmp << ") should be either true or false. Terminating the simulation.\n"; return -1; }
        }
        catch (PlbIOException& exception) { track_performance=0; }
        try { doc["parameters"]["LB_numerics"]["Peclet"].read(Pe); }
        catch (PlbIOException& exception) { Pe=0; }
        if (delta_P < thrd) {Pe=0;}
        try { doc["parameters"]["LB_numerics"]["tau"].read(tau); }
        catch (PlbIOException& exception) { tau=0.8; }
        try { doc["parameters"]["LB_numerics"]["domain"]["dy"].read(dy); }
        catch (PlbIOException& exception) { dy = dx; }
        try { doc["parameters"]["LB_numerics"]["domain"]["characteristic_length"].read(charcs_length); }
        catch (PlbIOException& exception) {
            charcs_length = 0;
            if (Pe>thrd) {
                pcout << "charcs_length must be defined when for transport simulations (Pe > 0). Terminating the simulation.\n"; return -1;
            }
        }
        try {
            std::string unit;
            doc["parameters"]["LB_numerics"]["domain"]["unit"].read(unit);
            if (unit == "m") {  charcs_length /= dx; /* do nothing */ }
            else if (unit == "mm") { charcs_length /= dx; dx *= 1e-3; }
            else if (unit == "um") { charcs_length /= dx; dx *= 1e-6; }
            else { pcout << "unit ("<< unit << ") must be either m, mm, or um. Terminating the simulation.\n"; return -1; }
        }
        catch (PlbIOException& exception) { charcs_length /= dx; dx *= 1e-6; }
        try { doc["parameters"]["LB_numerics"]["domain"]["material_numbers"]["pore"].read(pore_dynamics); }
        catch (PlbIOException& exception) { pore_dynamics.push_back(2); }
        try { doc["parameters"]["LB_numerics"]["domain"]["material_numbers"]["solid"].read(no_dynamics); }
        catch (PlbIOException& exception) { no_dynamics=0; }
        try { doc["parameters"]["LB_numerics"]["domain"]["material_numbers"]["bounce_back"].read(bounce_back); }
        catch (PlbIOException& exception) { bounce_back=1; }
        bfilm_count=0; bfree_count=0;
        for (plint iT=0; iT<num_of_microbes; ++iT) {
            std::string name = "microbe" + std::to_string(iT);
            try {
                std::vector<plint> tmp;
                doc["parameters"]["LB_numerics"]["domain"]["material_numbers"][name].read(tmp);
                bio_dynamics.push_back(tmp);
                bmass_type.push_back(1);
                vec_b0_film.push_back(vec_b0_all[iT]);
                if (tmp.size()!=vec_b0_all[iT].size()) {pcout << "The length of microbial material_number vector is not consistent with the length of initial_densities vector. Terminating the simulation.\n"; return -1;}
                ++bfilm_count;
            }
            catch (PlbIOException& exception) { bmass_type.push_back(0); vec_b0_free.push_back(vec_b0_all[iT][0]); ++bfree_count; }
        }
        try {
            std::string tmp;
            doc["parameters"]["IO"]["read_NS_file"].read(tmp);
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](unsigned char c){ return std::tolower(c); });
            if (tmp.compare("no")==0 || tmp.compare("false")==0 || tmp.compare("0")==0 ) { read_NS_file = 0; }
            else if (tmp.compare("yes")==0 || tmp.compare("true")==0 || tmp.compare("1")==0 ) { read_NS_file = 1; }
            else { pcout << "read_NS_file (" << tmp << ") should be either true or false. Terminating the simulation.\n"; return -1; }
        }
        catch (PlbIOException& exception) { read_NS_file=0; }
        try {
            doc["parameters"]["LB_numerics"]["iteration"]["ns_rerun_iT0"].read(ns_rerun_iT0);
            if (ns_rerun_iT0 < 0) {
                pcout << "ns_rerun_iT0 ("<< ns_rerun_iT0 << ") must be a positive number. Terminating the simulation.\n";
                return -1;
            }
        }
        catch (PlbIOException& exception) { if (read_NS_file == 1) { pcout << "WARNING: NS checkpoint file is loaded but ns_rerun_iT0 is not provided. Assume no further flow simulation.\n"; ns_rerun_iT0=0; } }
        try { doc["parameters"]["LB_numerics"]["iteration"]["ns_update_interval"].read(ns_update_interval); }
        catch (PlbIOException& exception) { ns_update_interval=1; }
        try { doc["parameters"]["LB_numerics"]["iteration"]["ade_update_interval"].read(ade_update_interval); }
        catch (PlbIOException& exception) { ade_update_interval=1; }
        try { doc["parameters"]["LB_numerics"]["iteration"]["ns_max_iT1"].read(ns_maxiTer_1); }
        catch (PlbIOException& exception) { ns_maxiTer_1=100000; }
        try { doc["parameters"]["LB_numerics"]["iteration"]["ns_max_iT2"].read(ns_maxiTer_2); }
        catch (PlbIOException& exception) { ns_maxiTer_2=100000; }
        try { doc["parameters"]["LB_numerics"]["iteration"]["ns_converge_iT1"].read(ns_converge_iT1); }
        catch (PlbIOException& exception) { ns_converge_iT1=1e-8; }
        try { doc["parameters"]["LB_numerics"]["iteration"]["ns_converge_iT2"].read(ns_converge_iT2); }
        catch (PlbIOException& exception) { ns_converge_iT2=1e-6; }
        try {
            doc["parameters"]["LB_numerics"]["iteration"]["ade_rerun_iT0"].read(ade_rerun_iT0);
            if (ade_rerun_iT0 < 0) {
                pcout << "ade_rerun_iT0 ("<< ade_rerun_iT0 << ") must be a positive number. Terminating the simulation.\n";
                return -1;
            }
        }
        catch (PlbIOException& exception) { if (read_ADE_file == 1) { pcout << "WARNING: ADE checkpoint file is loaded but ade_rerun_iT0 is not provided. Assume no further flow simulation.\n"; ade_rerun_iT0=0; } }
        try { doc["parameters"]["LB_numerics"]["iteration"]["ade_max_iT"].read(ade_maxiTer); }
        catch (PlbIOException& exception) { ade_maxiTer=10000000; }
        try { doc["parameters"]["LB_numerics"]["iteration"]["ade_converge_iT"].read(ade_converge_iT); }
        catch (PlbIOException& exception) { ade_converge_iT=1e-8; }

        // chemistry
        soluteDindex = 0;
        for (plint iT=0; iT<num_of_substrates; ++iT) {
            T D0, D1, c0;
            std::string chemname = "substrate"+std::to_string(iT);
            try { doc["parameters"]["chemistry"][chemname]["name_of_substrates"].read(vec_subs_names); }
            catch (PlbIOException& exception) { vec_subs_names.push_back("substrate_"+std::to_string(iT)); }
            try { doc["parameters"]["chemistry"][chemname]["substrate_diffusion_coefficients"]["in_pore"].read(D0); solute_poreD.push_back(D0);}
            catch (PlbIOException& exception) { solute_poreD.push_back(1e-9); }
            try { doc["parameters"]["chemistry"][chemname]["substrate_diffusion_coefficients"]["in_biofilm"].read(D1); solute_bFilmD.push_back(D1);}
            catch (PlbIOException& exception) { solute_bFilmD.push_back(2e-10); }
            if ( std::abs(D1-D0)> thrd ) { soluteDindex = 1; }
            try { doc["parameters"]["chemistry"][chemname]["initial_concentration"].read(c0); vec_c0.push_back(c0); }
            catch (PlbIOException& exception) { vec_c0.push_back(0.0); }
        }
        if (vec_c0.size()!=(unsigned)num_of_substrates) { pcout << "The length of initial_concentration vector does not match the num_of_substrates. Terminating the simulation.\n"; return -1;}
        if (solute_poreD.size()!=(unsigned)num_of_substrates) { pcout << "The length of substrate_diffusion_coefficients in_pore vector does not match the num_of_substrates. Terminating the simulation.\n"; return -1; }
        if (solute_bFilmD.size()!=(unsigned)num_of_substrates) { pcout << "The length of substrate_diffusion_coefficients in_biofilm vector does not match the num_of_substrates. Terminating the simulation.\n"; return -1; }
        if (vec_subs_names.size()!=(unsigned)num_of_substrates) { pcout << "The length of name_of_substrates vector does not match the num_of_substrates. Terminating the simulation.\n"; return -1; }
        if (mm_count > 0) {
            try {
                std::vector<std::string> tmp;
                doc["parameters"]["chemistry"]["fix_concentration"].read(tmp);
                if (tmp.size()!=(unsigned)num_of_substrates) { pcout << "The length of fix_concentration vector does not match the num_of_substrates. Terminating the simulation.\n"; return -1;}
                for (plint iT=0; iT<num_of_substrates; ++iT) {
                    std::transform(tmp[iT].begin(), tmp[iT].end(), tmp[iT].begin(), [](unsigned char c){ return std::tolower(c); });
                    if (tmp[iT].compare("no")==0 || tmp[iT].compare("false")==0 || tmp[iT].compare("0")==0 ) { vec_fixC.push_back(0); }
                    else if (tmp[iT].compare("yes")==0 || tmp[iT].compare("true")==0 || tmp[iT].compare("1")==0 ) { vec_fixC.push_back(1); }
                    else { pcout << "Element " << iT << " in fix_concentration (" << tmp[iT] << ") should be either true or false. Terminating the simulation.\n"; return -1; }
                }
            }
            catch (PlbIOException& exception) { for (plint iT=0; iT<num_of_substrates; ++iT) {vec_fixC.push_back(0);}  }
            try {
                std::vector<std::string> tmp;
                doc["parameters"]["chemistry"]["fix_lower_bounds"].read(tmp);
                if (tmp.size()!=(unsigned)num_of_substrates) { pcout << "The length of fix_lower_bounds vector does not match the num_of_substrates. Terminating the simulation.\n"; return -1;}
                for (plint iT=0; iT<num_of_substrates; ++iT) {
                    std::transform(tmp[iT].begin(), tmp[iT].end(), tmp[iT].begin(), [](unsigned char c){ return std::tolower(c); });
                    if (tmp[iT].compare("no")==0 || tmp[iT].compare("false")==0 || tmp[iT].compare("0")==0 ) { vec_fixLB.push_back(0); }
                    else if (tmp[iT].compare("yes")==0 || tmp[iT].compare("true")==0 || tmp[iT].compare("1")==0 ) { vec_fixLB.push_back(1); }
                    else { pcout << "Element " << iT << "in fix_lower_bounds (" << tmp[iT] << ") should be either true or false. Terminating the simulation.\n"; return -1; }
                }
            }
            catch (PlbIOException& exception) { for (plint iT=0; iT<num_of_substrates; ++iT) {vec_fixLB.push_back(0);}  }
        }

        // microbiology
        bmassDindex = 0;
        for (plint iT=0; iT<num_of_microbes; ++iT) {
            std::string bioname = "microbe"+std::to_string(iT);
            try {
                std::string tmp;
                doc["parameters"]["microbiology"][bioname]["name_of_microbes"].read(tmp);
                vec_microbes_names.push_back(tmp);
            }
            catch (PlbIOException& exception) { vec_microbes_names.push_back(bioname); }
            try {
                T mu;
                doc["parameters"]["microbiology"][bioname]["decay_coefficient"].read(mu);
                vec_mu.push_back(mu);
            }
            catch (PlbIOException& exception) { vec_mu.push_back(0.); }
            try {
                std::string tmp0, tmp1;
                doc["parameters"]["microbiology"][bioname]["left_boundary_type"].read(tmp0); tmp1=tmp0;
                std::transform(tmp0.begin(), tmp0.end(), tmp0.begin(), [](unsigned char c){ return std::tolower(c); });
                if (tmp0.compare("dirichlet")==0) { bio_left_btype.push_back(0); }
                else if (tmp0.compare("neumann")==0) { bio_left_btype.push_back(1); }
                else { pcout << "left_boundary_type (" << tmp1 << ") should be either Dirichlet or Neumann. Terminating the simulation.\n"; return -1; }
            }
            catch (PlbIOException& exception) { bio_left_btype.push_back(1); }
            try {
                std::string tmp0, tmp1;
                doc["parameters"]["microbiology"][bioname]["right_boundary_type"].read(tmp0); tmp1=tmp0;
                std::transform(tmp0.begin(), tmp0.end(), tmp0.begin(), [](unsigned char c){ return std::tolower(c); });
                if (tmp0.compare("dirichlet")==0) { bio_right_btype.push_back(0); }
                else if (tmp0.compare("neumann")==0) { bio_right_btype.push_back(1); }
                else { pcout << "right_boundary_type (" << tmp1 << ") should be either Dirichlet or Neumann. Terminating the simulation.\n"; return -1; }
            }
            catch (PlbIOException& exception) { bio_right_btype.push_back(1); }
            try {
                T bc;
                doc["parameters"]["microbiology"][bioname]["left_boundary_condition"].read(bc);
                bio_leftBC.push_back(bc);
            }
            catch (PlbIOException& exception) { bio_leftBC.push_back(0.); }
            try {
                T bc;
                doc["parameters"]["microbiology"][bioname]["right_boundary_condition"].read(bc);
                bio_rightBC.push_back(bc);
            }
            catch (PlbIOException& exception) { bio_rightBC.push_back(0.); }
            try {
                T D0;
                doc["parameters"]["microbiology"][bioname]["biomass_diffusion_coefficients"]["in_pore"].read(D0);
                bMass_poreD.push_back(D0);
            }
            catch (PlbIOException& exception) {
                if ( solver_type[iT]==1 ) { pcout << exception.what() << " for microbe" << iT << ". It must be defined when solver_type is Finite Difference. Terminating the simulation.\n"; return -1; }
                else if ( solver_type[iT]==3 ) { pcout << exception.what() << " It must be defined when solver_type is Lattice Boltzmann. Terminating the simulation.\n"; return -1; }
                else { bMass_poreD.push_back(-99); }
            }
            try {
                T D0;
                doc["parameters"]["microbiology"][bioname]["biomass_diffusion_coefficients"]["in_biofilm"].read(D0);
                bMass_bFilmD.push_back(D0);
            }
            catch (PlbIOException& exception) {
                if ( solver_type[iT]==1 ) { pcout << exception.what() << " for microbe" << iT << ". It must be defined when solver_type is Finite Difference. Terminating the simulation.\n"; return -1; }
                else if ( solver_type[iT]==3 ) { pcout << exception.what() << " It must be defined when solver_type is Lattice Boltzmann. Terminating the simulation.\n"; return -1; }
                else { bMass_bFilmD.push_back(-99); }
            }
            if ( (bMass_poreD[iT]>0&&bMass_poreD[iT]>0) && (std::abs(bMass_poreD[iT]-bMass_bFilmD[iT])>thrd) ) { bmassDindex = 1; }
            try {
                T nu0;
                doc["parameters"]["microbiology"][bioname]["viscosity_ratio_in_biofilm"].read(nu0);
                if (nu0 > thrd) { vec_permRatio.push_back(1/nu0); }
                else { vec_permRatio.push_back(-99); }
            }
            catch (PlbIOException& exception) {
                if (solver_type[iT]==2 ) { pcout << exception.what() << " It must be defined when solver_type is Cellular Automata. Terminating the simulation.\n"; return -1; }
            }
            try {
                std::vector<T> tmp;
                doc["parameters"]["microbiology"][bioname]["half_saturation_constants"].read(tmp);
                if (tmp.size()!=(unsigned)num_of_substrates) { pcout << "The length of half_saturation_constants should be equal to the number_of_substates. Terminating the simulation.\n"; }
                else { vec_Kc.push_back(tmp); }
            }
            catch (PlbIOException& exception) {
                if (reaction_type[iT]==0 ) { pcout << exception.what() << " half_saturation_constants must be defined when reaction_type is GLPK. Terminating the simulation.\n"; return -1; }
                else if (reaction_type[iT]==2 ) { pcout << exception.what() << " half_saturation_constants must be defined when reaction_type is NN. Terminating the simulation.\n"; return -1; }
                else { vec_Kc.push_back( std::vector<T> (1,-99) ); }
            }
            try {
                std::vector<T> vmax;
                doc["parameters"]["microbiology"][bioname]["maximum_uptake_flux"].read(vmax);
                vec_Vmax.push_back(vmax);
            }
            catch (PlbIOException& exception) { vec_Vmax.push_back( std::vector<T> (num_of_substrates,0) ); }
            try {
                std::vector<int> tmp;
                doc["parameters"]["microbiology"][bioname]["constraint_indices"].read(tmp);
                vec_const_loc.push_back(tmp);
            }
            catch (PlbIOException& exception) { vec_const_loc.push_back( std::vector<int> (1,-99) ); }
            try {
                std::vector<T> tmp;
                doc["parameters"]["microbiology"][bioname]["constraint_lower_bounds"].read(tmp);
                vec_const_lb.push_back(tmp);
            }
            catch (PlbIOException& exception) { vec_const_lb.push_back( std::vector<T> (1,-99) ); }
            try {
                std::vector<T> tmp;
                doc["parameters"]["microbiology"][bioname]["constraint_upper_bounds"].read(tmp);
                vec_const_ub.push_back(tmp);
            }
            catch (PlbIOException& exception) { vec_const_ub.push_back( std::vector<T> (1,-99) ); }
        }
        try { doc["parameters"]["microbiology"]["thrd_biofilm_fraction"].read(thrd_bFilmFrac); }
        catch (PlbIOException& exception) { if (ca_count>0 ) { pcout << exception.what() << " It must be defined when solver_type is Cellular Automata. Terminating the simulation.\n"; return -1; } }
        if (vec_microbes_names.size()!=(unsigned)num_of_microbes) { pcout << "The length of name_of_microbes vector does not match the num_of_microbes. Terminating the simulation.\n"; return -1; }
        if (vec_mu.size()!=(unsigned)num_of_microbes) { pcout << "The length of decay_coefficient vector does not match the num_of_microbes. Terminating the simulation.\n"; return -1; }
        if (bMass_poreD.size()!=(unsigned)num_of_microbes) { pcout << "The length of biomass_diffusion_coefficients in_pore vector does not match the num_of_microbes. Terminating the simulation.\n"; return -1; }
        if (bMass_bFilmD.size()!=(unsigned)num_of_microbes) { pcout << "The length of biomass_diffusion_coefficients in_biofilm vector does not match the num_of_microbes. Terminating the simulation.\n"; return -1; }
        if (vec_permRatio.size()!=(unsigned)bfilm_count) { pcout << "The length of biofilm_permeability_ratio vector does not match the length of the material number vector for microbes. Terminating the simulation.\n"; return -1; }

        try { doc["parameters"]["microbiology"]["maximum_biomass_density"].read(max_bMassRho); }
        catch (PlbIOException& exception) {
            if (ca_count==1 ) { pcout << exception.what() << " It must be defined when solver_type is Cellular Automata. Terminating the simulation.\n"; return -1; }
            else { max_bMassRho = 999999999.; }
        }
        try {
            std::string tmp;
            doc["parameters"]["microbiology"]["CA_method"].read(tmp);
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](unsigned char c){ return std::tolower(c); });
            if (tmp.compare("fraction")==0 || tmp.compare("0")==0 || tmp.compare("no")==0 ) { halfflag=0; }
            else if (tmp.compare("half")==0 || tmp.compare("1")==0 || tmp.compare("yes")==0 ) { halfflag=1; }
            else { pcout << "halfflag (" << tmp << ") should be either half or fraction. Terminating the simulation.\n"; return -1; }
        }
        catch (PlbIOException& exception) { halfflag=0; }

        // IO
        try {
            std::string tmp;
            doc["parameters"]["IO"]["read_ADE_file"].read(tmp);
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](unsigned char c){ return std::tolower(c); });
            if (tmp.compare("no")==0 || tmp.compare("false")==0 || tmp.compare("0")==0 ) { read_ADE_file = 0; }
            else if (tmp.compare("yes")==0 || tmp.compare("true")==0 || tmp.compare("1")==0 ) { read_ADE_file = 1; }
            else { pcout << "read_ADE_file (" << tmp << ") should be either true or false. Terminating the simulation.\n"; return -1; }
        }
        catch (PlbIOException& exception) { read_ADE_file=0; }
        try {
            std::string item;
            doc["parameters"]["IO"]["ns_filename"].read(item);
            ns_filename = (char *) calloc(item.size()+1,sizeof(char));
            for (size_t i = 0; i < item.size(); ++i) { ns_filename[i] = item[i]; }
            ns_filename[item.size()+1] = '\0';
        }
        catch (PlbIOException& exception) {
            std::string item = "nsLattice";
            ns_filename = (char *) calloc(item.size()+1,sizeof(char));
            for (size_t i = 0; i < item.size(); ++i) { ns_filename[i] = item[i]; }
            ns_filename[item.size()+1] = '\0';
        }
        try { doc["parameters"]["IO"]["mask_filename"].read(mask_filename); }
        catch (PlbIOException& exception) { mask_filename="maskLattice"; }
        try { doc["parameters"]["IO"]["subs_filename"].read(ade_filename); }
        catch (PlbIOException& exception) { ade_filename="subsLattice"; }
        try { doc["parameters"]["IO"]["bio_filename"].read(bio_filename); }
        catch (PlbIOException& exception) { bio_filename="bioLattice"; }
        try { doc["parameters"]["IO"]["save_VTK_interval"].read(ade_VTK_iTer); }
        catch (PlbIOException& exception) { ade_VTK_iTer=1000; }
        try { doc["parameters"]["IO"]["save_CHK_interval"].read(ade_CHK_iTer); }
        catch (PlbIOException& exception) { ade_CHK_iTer=1000000; }

        if (vec_Kc.size()>0) {
            for (plint iT=0; iT<num_of_microbes; ++iT) {
                std::vector<T> Vmax0=vec_Vmax[iT]; std::vector<T> Kc0=vec_Kc[iT]; T mu0=vec_mu[iT];
                if (reaction_type[iT]==0) { vec_Kc_glpk.push_back(Kc0); vec_mu_glpk.push_back(mu0); vec_Vmax_glpk.push_back(Vmax0); }
                else if (reaction_type[iT]==1) { vec_Kc_kns.push_back(Kc0); vec_mu_kns.push_back(mu0); vec_Vmax_kns.push_back(Vmax0); }
                else if (reaction_type[iT]==2) { vec_Kc_nn.push_back(Kc0); vec_mu_nn.push_back(mu0); vec_Vmax_nn.push_back(Vmax0); }
                else if (reaction_type[iT]==3) { vec_Kc_glpk.push_back(Kc0); vec_mu_glpk.push_back(mu0); vec_Vmax_glpk.push_back(Vmax0); vec_Kc_kns.push_back(Kc0); vec_mu_kns.push_back(mu0); vec_Vmax_kns.push_back(Vmax0); }
                else if (reaction_type[iT]==4) { vec_Kc_nn.push_back(Kc0); vec_mu_nn.push_back(mu0); vec_Vmax_nn.push_back(Vmax0); vec_Kc_kns.push_back(Kc0); vec_mu_kns.push_back(mu0); vec_Vmax_kns.push_back(Vmax0); }
                else { pcout << "ERROR: undefined reaction_type. Terminating the simulation.\n"; return -1;}
            }
        }
    }
    catch (PlbIOException& exception) {
        pcout << exception.what() << " Terminating the simulation.\n" << std::endl;
        return -1;
    }
    return 0;
}

int load_metabolic_models(std::string input_path, plint num_of_microbes, std::vector<std::string> vec_mmFileName,
                          std::vector<std::vector<std::vector<T>>> &vec3_S, std::vector<std::vector<T>> &vec2_b,
                          std::vector<std::vector<T>> &vec2_c, std::vector<std::vector<T>> &vec2_lb, std::vector<std::vector<T>> &vec2_ub, std::vector<plint> &vec1_objLoc,
                          std::vector<std::vector<int>> vec2_const_loc, std::vector<std::vector<T>> vec2_const_lb, std::vector<std::vector<T>> vec2_const_ub)
{
    pcout << "loading metabolic models ..." << std::endl;
    for (plint iM=0; iM<num_of_microbes; ++iM) {
        try {
            std::size_t found = vec_mmFileName[iM].find_last_of(".");
            std::string ifile;
            if (found == std::string::npos) { ifile = input_path+vec_mmFileName[iM]+".xml"; }
            else {
                if ( vec_mmFileName[iM].substr(found+1).compare("xml")==0 ) {
                  ifile = input_path+vec_mmFileName[iM];
                }
                else {
                    pcout << "ERROR: metabolic input file extension must be xml. Terminating the simulation.\n";
                    return -1;
                }
            }
            std::string fin(ifile);
            XMLreader doc(fin);

            plint nrow, ncol, obj;
            std::vector<T> s1,b,c,lb,ub;
            doc["Metabolic_Model"]["nmet"].read(nrow);
            doc["Metabolic_Model"]["nrxn"].read(ncol);
            try { doc["Metabolic_Model"]["objLoc"].read(obj); }
            catch (PlbIOException& exception) { obj = 0; }
            doc["Metabolic_Model"]["S"].read(s1);
            doc["Metabolic_Model"]["b"].read(b);
            doc["Metabolic_Model"]["c"].read(c);
            doc["Metabolic_Model"]["lb"].read(lb);
            doc["Metabolic_Model"]["ub"].read(ub);
            if (vec2_const_loc[iM][0] > 0) {
                for (size_t iV=0; iV<vec2_const_loc[iM].size(); ++iV) {
                    int loc = vec2_const_loc[iM][iV];
                    if (vec2_const_lb[iM][iV]>0) { lb[loc] = vec2_const_lb[iM][iV]; }
                    if (vec2_const_ub[iM][iV]>0) { ub[loc] = vec2_const_ub[iM][iV]; }
                }
            }

            std::vector< std::vector<T> > s2(nrow);
            for (plint irow=0; irow<nrow; ++irow) {
                s2[irow] = std::vector<T> (ncol);
                for (plint icol=0; icol<ncol; ++icol) {
                    s2[irow][icol]=s1[irow*ncol+icol];
                }
            }
            vec3_S[iM]=s2;
            vec2_b[iM]=b;
            vec2_c[iM]=c;
            vec2_lb[iM]=lb;
            vec2_ub[iM]=ub;
            vec1_objLoc[iM]=obj;
        }
        catch (PlbIOException& exception) {
            pcout << exception.what() << std::endl;
            return -1;
        }
    }
    return 0;
}

// void prep_glpk(std::string str_inputDir, plint num_of_microbes, std::vector<std::string> vec_mmFileName, std::vector<int> &vec_nrow, std::vector<int> &vec_ncol, std::vector<int> &vec_method,
//               std::vector<int> &vec_isMIP, std::vector<glp_prob *> &vec_lp, std::vector<glp_smcp> &vec_sParam, std::vector<glp_iocp> &vec_iParam, std::vector<plint> &vec_objLoc)
// {
//     std::vector < std::vector < std::vector<T> > > S(num_of_microbes);
//     std::vector < std::vector<T> > vec_b(num_of_microbes), vec_c(num_of_microbes), vec_lb(num_of_microbes), vec_ub(num_of_microbes);

//     load_metabolic_models(str_inputDir, num_of_microbes, vec_mmFileName, S, vec_b, vec_c, vec_lb, vec_ub, vec_objLoc);

//     for (plint iM = 0; iM < num_of_microbes; ++iM) {
//         int nrow = S[iM].size();
//         int ncol = S[iM][0].size();
//         vec_nrow.push_back(nrow);
//         vec_ncol.push_back(ncol);
//     }

//     pcout << "Setting up the glpk environments ... " << std::endl;
//     for (plint iM = 0; iM < num_of_microbes; ++iM) {
//         // setting up the glpk environments

//         // ctype = A row array containing the sense of each constraint in the constraint matrix.
//         // 'F' Free (unbounded) variable (the constraint is ignored).
//         // 'U' Variable with upper bound ( A(i,:)*x <= b(i)).
//         // 'S' Fixed Variable (A(i,:)*x = b(i)).
//         // 'L' Variable with lower bound (A(i,:)*x >= b(i)).
//         // 'D' Double-bounded variable (A(i,:)*x >= -b(i) and A(i,:)*x <= b(i)).
//         char ctype[vec_nrow[iM]+1];
//         memset( ctype, 'S', (vec_nrow[iM]+1)*sizeof(char) );

//         // vartype = A column array containing the types of the variables.
//         // 'C' Continuous variable.
//         // 'I' Integer variable
//         // 'B' Binary variable
//         char vtype[vec_ncol[iM]+1];
//         memset( vtype, 'C', (vec_ncol[iM]+1)*sizeof(char) );

//         //-- Sense of optimization (maximize = -1, minimize = 1).
//         plint sense = -1;

//         // lpsolver: Select which solver to use.
//         // This flag will be ignored if the problem is a MIP problem .
//         // 1 - Revised simplex method.
//         // 2 - Interior point method.
//         // 3 - Simplex method with exact arithmatic.
//         plint lpsolver = 1;

//         // save_pb : Save a copy of the original problem to file name specified below
//         // 0 - no save
//         // 1 - cplex (glpk_output.lp)
//         // 2 - fixedmps (glpk_output.mps)
//         // 3 - freemps  (glpk_output.mps)
//         // 4 - plain (glpk_output.txt)
//         plint save_pb = 0;

//         // Create an empty LP/MILP object
//         vec_lp.push_back(glp_create_prob ());
//         // glp_prob *lp = glp_create_prob ();

//         int method, isMIP;
//         glp_smcp sParam;
//         glp_iocp iParam;

//         initialize_glpk(iM, vec_lp[iM], S[iM], vec_b[iM], vec_c[iM], vec_lb[iM], vec_ub[iM], ctype, vtype, sense, lpsolver, save_pb, method, isMIP, sParam, iParam);

//         vec_method.push_back(method);
//         vec_isMIP.push_back(isMIP);
//         vec_sParam.push_back(sParam);
//         vec_iParam.push_back(iParam);
//     }
// }
