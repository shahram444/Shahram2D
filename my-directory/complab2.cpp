#include <glpk.h>
#include "glpkcpp.hh"
#include "complab_functions.hh"
#include "complab_processors.hh"
#include "axiv_processor.hh"

#include <chrono>
#include <string>
#include <iostream>
#include <cstring>
#include <vector>
#include <sys/stat.h>

int main(int argc, char **argv) {

    plbInit(&argc, &argv);
    ImageWriter<T> image("leeloo");
    auto t_start_0 = std::chrono::high_resolution_clock::now();

    // asserted variables
    plint mm_count=0, nn_count=0, kns_count=0, fd_count=0, lb_count=0, ca_count=0, bfilm_count=0, bfree_count=0;
    char *main_path = (char*)malloc(100 * sizeof(char));
    getcwd(main_path, 100 * sizeof(char));
    // char *main_path = (char*)malloc(100 * sizeof(char));
    char *input_path = (char*)malloc(100 * sizeof(char));
    char *output_path = (char*)malloc(100 * sizeof(char));
    char *ns_filename = (char*)malloc(100 * sizeof(char));
    plint nx, ny, num_of_microbes, num_of_substrates;
    T dx, dy, deltaP, Pe, charcs_length;
    std::string geom_filename, mask_filename;
    std::vector<bool> vec_left_btype, vec_right_btype, bio_left_btype, bio_right_btype;
    std::vector<T> vec_c0, vec_b0_free, vec_left_bcondition, vec_right_bcondition, bio_left_bcondition, bio_right_bcondition, vec_permRatio;
    std::vector< std::vector<int> > vec_EX_loc, vec_const_loc;
    std::vector< std::vector<T> > vec_b0_all, vec_b0_film, vec_Kc, vec_Kc_glpk, vec_Kc_kns, vec_Kc_nn, vec_maxUptake, vec_maxRelease, vec_Vmax, vec_Vmax_glpk, vec_Vmax_kns, vec_Vmax_nn, vec_const_ub, vec_const_lb;
    std::vector< std::string > vec_mmFileName;

    // variables with default values
    std::string ade_filename, bio_filename;
    bool read_NS_file=0, read_ADE_file=0, soluteDindex=0, bmassDindex=0, track_performance=0., halfflag=0;
    plint no_dynamics=0, bounce_back=1, ns_rerun_iT0=0, ns_update_interval=1, ade_update_interval=1,
        ns_maxiTer_1, ns_maxiTer_2, ade_rerun_iT0=0, ade_maxiTer=10000000, ade_VTI_iTer=1000, ade_CHK_iTer=1000000;
    T tau=0.8, max_bMassRho=1., ns_converge_iT1=1e-8, ns_converge_iT2=1e-4, ade_converge_iT=1e-8, thrd_bFilmFrac;
    std::vector<bool> vec_fixC, vec_fixLB, bmass_type;
    std::vector<plint> pore_dynamics, solver_type, reaction_type, vec_sense;
    std::vector<T> vec_solute_poreD, vec_solute_bFilmD, vec_bMass_poreD, vec_bMass_bFilmD, vec_mu, vec_mu_glpk, vec_mu_kns, vec_mu_nn;
    std::vector<std::string> vec_subs_names, vec_microbes_names;
    std::vector< std::vector<plint> > bio_dynamics;

    std::string str_mainDir=main_path;
    if (std::to_string(str_mainDir.back()).compare("/")!=0) { str_mainDir+="/"; }
    std::srand(std::time(nullptr));

    //-------------------- load the input file and initialize all input parameters --------------------//
    int erck = 0;
    try {
        erck=initialize_complab( main_path, input_path, output_path, ns_filename, ade_filename, bio_filename, geom_filename, mask_filename,
        read_NS_file, ns_rerun_iT0, ns_converge_iT1, ns_converge_iT2, ns_maxiTer_1, ns_maxiTer_2, ns_update_interval, ade_update_interval,
        read_ADE_file, ade_rerun_iT0, ade_VTI_iTer, ade_CHK_iTer, ade_converge_iT, ade_maxiTer, nx, ny, dx, dy, deltaP, tau,
        Pe, charcs_length, vec_solute_poreD, vec_solute_bFilmD, vec_bMass_poreD, vec_bMass_bFilmD, soluteDindex, bmassDindex, thrd_bFilmFrac, vec_permRatio, max_bMassRho,
        pore_dynamics, bounce_back, no_dynamics, bio_dynamics, num_of_microbes, num_of_substrates, vec_subs_names, vec_microbes_names,
        solver_type, fd_count, lb_count, ca_count, bfilm_count, bfree_count, vec_mmFileName, mm_count, kns_count, nn_count, reaction_type, vec_maxUptake, vec_maxRelease, vec_EX_loc,
        vec_c0, vec_left_btype, vec_right_btype, vec_left_bcondition, vec_right_bcondition, vec_b0_all, bio_left_btype, bio_right_btype, bio_left_bcondition, bio_right_bcondition,
        vec_Kc, vec_Kc_glpk, vec_Kc_kns, vec_Kc_nn, vec_mu, vec_mu_glpk, vec_mu_kns, vec_mu_nn, vec_fixLB, vec_fixC, bmass_type, vec_b0_free, vec_b0_film, vec_sense,
        vec_Vmax, vec_Vmax_glpk, vec_Vmax_kns, vec_Vmax_nn, vec_const_loc, vec_const_lb, vec_const_ub, track_performance, halfflag);
    }
    catch (PlbIOException& exception) {
        pcout << exception.what() << " Terminating the simulation.\n" << std::endl;
        return -1;
    }
    if (erck!=0) { return -1; }
    // plint maxiTwhile=util::roundToInt(charcs_length);

    std::string  str_inputDir=input_path, str_outputDir=output_path;
    if (std::to_string(str_inputDir.back()).compare("/")!=0) { str_inputDir+="/"; }
    if (std::to_string(str_outputDir.back()).compare("/")!=0) { str_outputDir+="/"; }

    // ------------------------------------ load metabolic models ------------------------------------ //

    std::vector<int> vec_method(mm_count), vec_isMIP(mm_count), vec_ncol(mm_count), vec_nrow(mm_count);
    std::vector<glp_prob *> vec_lp(mm_count);
    std::vector<glp_smcp> vec_sParam(mm_count);
    std::vector<glp_iocp> vec_iParam(mm_count);
    std::vector<plint> vec_objLoc(mm_count);

    if (mm_count > 0) {
        // prep_glpk does not work for some reasons. For now, make it explicit
        // prep_glpk(str_inputDir, mm_count, vec_mmFileName, vec_nrow, vec_ncol, vec_method, vec_isMIP, vec_lp, vec_sParam, vec_iParam, vec_objLoc);
        std::vector < std::vector < std::vector<T> > > S(mm_count);
        std::vector < std::vector<T> > vec_b(mm_count), vec_c(mm_count), vec_lb(mm_count), vec_ub(mm_count);

        load_metabolic_models(str_inputDir, mm_count, vec_mmFileName, S, vec_b, vec_c, vec_lb, vec_ub, vec_objLoc, vec_const_loc, vec_const_lb, vec_const_ub);

        for (plint iM = 0; iM < mm_count; ++iM) {
            int nrow = S[iM].size();
            int ncol = S[iM][0].size();
            vec_nrow[iM]=nrow;
            vec_ncol[iM]=ncol;
        }

        pcout << "Setting up the glpk environments ... " << std::endl;
        for (plint iM = 0; iM < mm_count; ++iM) {
            // ctype = A row array containing the sense of each constraint in the constraint matrix.
            // 'F' Free (unbounded) variable (the constraint is ignored).
            // 'U' Variable with upper bound ( A(i,:)*x <= b(i)). L
            // 'S' Fixed Variable (A(i,:)*x = b(i)). E
            // 'L' Variable with lower bound (A(i,:)*x >= b(i)). G
            // 'D' Double-bounded variable (A(i,:)*x >= -b(i) and A(i,:)*x <= b(i)).
            char ctype[vec_nrow[iM]+1];
            memset( ctype, 'S', (vec_nrow[iM]+1)*sizeof(char) );

            // vartype = A column array containing the types of the variables.
            // 'C' Continuous variable.
            // 'I' Integer variable
            // 'B' Binary variable
            char vtype[vec_ncol[iM]+1];
            memset( vtype, 'C', (vec_ncol[iM]+1)*sizeof(char) );

            // lpsolver: Select which solver to use.
            // This flag will be ignored if the problem is a MIP problem .
            // 1 - Revised simplex method.
            // 2 - Interior point method.
            // 3 - Simplex method with exact arithmetic.
            plint lpsolver = 1;

            // save_pb : Save a copy of the original problem to file name specified below
            // 0 - no save
            // 1 - cplex (glpk_output.lp)
            // 2 - fixedmps (glpk_output.mps)
            // 3 - freemps  (glpk_output.mps)
            // 4 - plain (glpk_output.txt)
            plint save_pb = 0;

            // Create an empty LP/MILP object
            glp_prob *lp = glp_create_prob ();

            int method, isMIP;
            glp_smcp sParam;
            glp_iocp iParam;

            initialize_glpk(iM, lp, S[iM], vec_b[iM], vec_c[iM], vec_lb[iM], vec_ub[iM], ctype, vtype, vec_sense[iM], lpsolver, save_pb, method, isMIP, sParam, iParam);

            vec_lp[iM]=lp;
            vec_method[iM]=method;
            vec_isMIP[iM]=isMIP;
            vec_sParam[iM]=sParam;
            vec_iParam[iM]=iParam;
        }
    }

    /*  =================================== NS Lattice setup ===================================  */

    struct stat statStruct;
    stat(output_path, &statStruct);

    pcout << "CompLaB main directory = " << str_mainDir << std::endl;
    pcout << "CompLaB input directory = " << main_path << "/" << input_path << std::endl;
    pcout << "CompLaB output directory = " << main_path << "/" << output_path << std::endl << std::endl;
    if (S_ISDIR(statStruct.st_mode)) {} else { mkdir(output_path, 0777); }
    global::directories().setOutputDir(str_outputDir);

    // NS parameters
    // T DarcyOutletUx, DarcyMiddleUx, DarcyInletUx;
    T PoreMeanU=0, PoreMaxUx=0;
    plint iT = 0;
    T nsLatticeTau = tau;
    T nsLatticeOmega = 1 / nsLatticeTau;
    T nsLatticeNu = NSDES<T>::cs2*(nsLatticeTau-0.5);
    char *ns_read_filename = strcat(strdup(str_inputDir.c_str()),ns_filename);

    MultiScalarField2D<int> geometry(nx,ny);
    readGeometry(str_inputDir+geom_filename, geometry);
    saveGeometry("inputGeom.vti", geometry);

    MultiScalarField2D<int> distanceDomain(nx,ny);
    distanceDomain = geometry;
    std::vector< std::vector<plint> > distVec(nx);
    for (plint iX=0; iX<nx; ++iX) { distVec[iX]=std::vector<plint> (ny); }
    calculateDistanceFromSolid(distanceDomain, no_dynamics, bounce_back, distVec);
    applyProcessingFunctional(new createDistanceDomain<int> (distVec), distanceDomain.getBoundingBox(), distanceDomain);
    if (track_performance == 0) { writeScalarVTI(distanceDomain); }

    MultiScalarField2D<int> ageDomain(nx,ny);
    ageDomain = geometry;
    applyProcessingFunctional(new createAgeDomain<int> (pore_dynamics, bounce_back, no_dynamics), ageDomain.getBoundingBox(), ageDomain);
    if (track_performance == 1) { pcout << "Performance tracker has been activated. Skipping all the non-essential IO.\n"; }

    pcout << "\nInitializing the fluid lattice (deltaP = " << deltaP << ").\n";
    MultiBlockLattice2D<T,NSDES> nsLattice(nx, ny, new IncBGKdynamics<T,NSDES>(nsLatticeOmega));
    // MultiBlockLattice2D<T,NSDES> nsLattice(nx, ny, new IncMRTdynamics<T,NSDES>(nsLatticeOmega));
    util::ValueTracer<T> ns_convg1(1.0,1000.0,ns_converge_iT1);
    NSdomainSetup( nsLattice, createLocalBoundaryCondition2D<T,NSDES>(), geometry, deltaP, nsLatticeOmega, pore_dynamics, bounce_back, no_dynamics, bio_dynamics, vec_permRatio);

    //-------------------- NS Lattice main loop --------------------//
    auto t_NS_start = std::chrono::high_resolution_clock::now();
    if ( Pe == 0 ) { pcout << "Peclet number is set to 0. Skipping a lattice Boltzmann flow solver.\n"; }
    else {
        pcout << "nsLatticeTau = " << nsLatticeTau <<", nsLatticeOmega = " << nsLatticeOmega << ", nsLatticeNu = " << nsLatticeNu << std::endl;
        pcout << "\n========== LBM NS simulation begins ==========\n";
        if (read_NS_file == 1 && track_performance == 0) {
            pcout << "run continuous simulation for nsLattice." << std::endl;
            pcout << "Load binary block for nsLattice." << std::endl;
            try { loadBinaryBlock(nsLattice, strcat(ns_read_filename,".chk")); }
            catch (PlbIOException& exception) { pcout << exception.what() << ". Terminating the simulation.\n" << std::endl; return -1; }

            if (ns_rerun_iT0==0) { pcout << "Use the existing checkpoint file for nsLattice." << std::endl; }
            else {
                pcout << "run the main nsLattice loop to a new steady state" << std::endl;
                iT = ns_rerun_iT0;
                for (; iT < ns_maxiTer_1; ++iT) {
                    nsLattice.collideAndStream();
                    ns_convg1.takeValue(getStoredAverageEnergy(nsLattice),true);
                    if (ns_convg1.hasConverged()) { break; }
                }
            }
        }
        else {
            pcout << "Run a new simulation for nsLattice" << std::endl;
            for (; iT < ns_maxiTer_1; ++iT) {
                nsLattice.collideAndStream();
                ns_convg1.takeValue(getStoredAverageEnergy(nsLattice),true);
                if (ns_convg1.hasConverged()) { break; }
            }
        }
        if ( read_NS_file==0 || (read_NS_file==1 && ns_rerun_iT0>0) ) {
            pcout << std::endl << "flow calc finished at iT = " << iT << std::endl;
            if (track_performance == 0) {
                pcout << "Writing velocity VTI... \n";
                writeNsVTI(nsLattice,ns_maxiTer_1,"nsLatticeFinal1_");
                pcout << "Writing checkpoint... \n";
                saveBinaryBlock(nsLattice, str_outputDir+ns_filename+".chk");
            }
        }

        // calculate the mean pore velocity. this is necessary for ADE simulation
        if (bfilm_count > 0) {
            plint totalCount = 0;
            T totalVel = 0;
            for (size_t iT = 0; iT < pore_dynamics.size(); ++iT) {
                plint poreCount = MaskedScalarCounts(Box2D (1,nx-2,0,ny-1), geometry, pore_dynamics[iT]);
                totalCount += poreCount;
                totalVel += ( computeAverage( *computeVelocityNorm(nsLattice, Box2D (1,nx-2,0,ny-1)), geometry, pore_dynamics[iT] ) * poreCount );
            }
            for (plint iT0 = 0; iT0 < bfilm_count; ++iT0) {
                plint bFilmCount = 0;
                for (size_t iT1 = 0; iT1 < bio_dynamics[iT0].size(); ++iT1) {
                    bFilmCount += MaskedScalarCounts(Box2D (1,nx-2,0,ny-1), geometry, bio_dynamics[iT0][iT1]);
                }
                totalCount += bFilmCount;
                totalVel += ( computeAverage( *computeVelocityNorm(nsLattice, Box2D (1,nx-2,0,ny-1)), geometry, bio_dynamics[iT0][0] ) * bFilmCount );
            }
            PoreMeanU = totalVel / totalCount;
        }
        else {
            PoreMeanU = computeAverage( *computeVelocityNorm(nsLattice, Box2D (1,nx-2,0,ny-1)) );
        }

        PoreMaxUx = computeMax( *computeVelocityComponent(nsLattice, Box2D (1,nx-2,0,ny-1), 0) );
        T Ma = PoreMaxUx/sqrt(RXNDES<T>::cs2);
        pcout << "CFL number (= maximum local lattice velocity)= " << PoreMaxUx << ".\n";
        pcout << "Mach number = " << Ma << ".\n";
        if (Ma > 1) { pcout << "Ma must be << 1. Terminating the simulation." << std::endl; return -1; }
        // DarcyOutletUx = computeAverage( *computeVelocityComponent(nsLattice, Box2D (nx-2,nx-2, 0,ny-1)) );
        // DarcyMiddleUx = computeAverage( *computeVelocityComponent(nsLattice, Box2D ((nx-1)/2,(nx-1)/2, 0,ny-1)) );
        // DarcyInletUx = computeAverage( *computeVelocityComponent(nsLattice, Box2D (1,1, 0,ny-1)) );
        // pcout << "Outlet Darcy Ux = " << DarcyOutletUx << " lu/ts" << std::endl;
        // pcout << "Middle Darcy Ux = " << DarcyMiddleUx << " lu/ts" << std::endl;
        // pcout << "Inlet Darcy Ux = " << DarcyInletUx << " lu/ts" << std::endl;
    }
    T nstime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now()-t_NS_start).count();

    if ( ade_maxiTer == 0) {
        pcout << "ade_max_iTer is set to 0. Terminating the simulation.\n";
        return 0;
    }

    /*  =================================== rxn Lattice setup  ===================================  */

    T refNu, refTau;
    if (Pe > thrd) {
        refNu = PoreMeanU * charcs_length / Pe;
        refTau = refNu * RXNDES<T>::invCs2 + 0.5;
        if (refTau > 2) {
            pcout << "Reference relaxation time is > 2 (refTau = " << refTau << "). Consider reducing it for numerical accuracy by reducing average flow velocity (e.g. reduce delta_P).\n";
            return -1;
        }
        else if (refTau <= 0.5) {
            pcout << "Reference relaxation time does not satisfy a necessary stability condition for the BGK operator. (tau must be > 0.5, but refTau = " << refTau << ").\n";
            pcout << "Consider increasing average flow velocity (e.g. increase delta_P).\n";
            return -1;
        }
    }
    else {
        refTau = tau;
        refNu = RXNDES<T>::cs2 * (refTau - 0.5);
    }
    T refOmega = 1/refTau;
    T ade_dt = refNu * dx * dx / vec_solute_poreD[0];

    std::vector<T> substrNUinPore(num_of_substrates), substrTAUinPore(num_of_substrates), substrOMEGAinPore(num_of_substrates), substrOMEGAinbFilm(num_of_substrates);
    for (plint iS = 0; iS < num_of_substrates; ++iS) {
        if (iS == 0) {
            substrNUinPore[iS]=refNu; substrTAUinPore[iS]=refTau; substrOMEGAinPore[iS]=refOmega;
        }
        else {
            substrNUinPore[iS]   = substrNUinPore[0]*vec_solute_poreD[iS]/vec_solute_poreD[0];
            substrTAUinPore[iS]  = substrNUinPore[iS]*RXNDES<T>::invCs2+0.5;
            substrOMEGAinPore[iS]= 1/substrTAUinPore[iS];
        }
        substrOMEGAinbFilm[iS] = 1/(refNu*vec_solute_bFilmD[iS]/vec_solute_poreD[0]*RXNDES<T>::invCs2+0.5);
    }

    std::vector<T> bioNUinPore(num_of_microbes), bioTAUinPore(num_of_microbes), bioOMEGAinPore(num_of_microbes), bioOMEGAinbFilm(num_of_microbes), bioTAUinbFilm(num_of_microbes);
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        if (vec_bMass_poreD[iM] > 0) {
            bioNUinPore[iM] = refNu * vec_bMass_poreD[iM] / vec_solute_poreD[0];
            bioTAUinPore[iM] = bioNUinPore[iM] * RXNDES<T>::invCs2 + 0.5;
            bioOMEGAinPore[iM] = 1/bioTAUinPore[iM];
        }
        else { bioNUinPore[iM] = 0.; bioTAUinPore[iM] = 0.; bioOMEGAinPore[iM] = 0.; }
        if (vec_bMass_bFilmD[iM] > 0) {
            bioOMEGAinbFilm[iM] = 1/(refNu*vec_bMass_bFilmD[iM]/vec_bMass_poreD[iM]*RXNDES<T>::invCs2+0.5);
            bioTAUinbFilm[iM] = 1/bioOMEGAinbFilm[iM];
        }
        else { bioOMEGAinbFilm[iM] = 0.; bioTAUinbFilm[iM] = 0.; }
    }

    pcout << "\nInitializing the reaction lattices... \n";
    pcout << "substrTAUinPore = ";
    for (plint iS = 0; iS < num_of_substrates; ++iS) {
        pcout << substrTAUinPore[iS] << " ";
    }
    pcout << std::endl;
    pcout << "substrTAUinbFilm = ";
    for (plint iS = 0; iS < num_of_substrates; ++iS) {
        pcout << 1/substrOMEGAinbFilm[iS] << " ";
    }
    pcout << std::endl;
    pcout << "bioTAUinPore = ";
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        pcout << bioTAUinPore[iM] << " ";
        if ( vec_bMass_poreD[iM] > 0 && bioTAUinPore[iM] <= 0.5 ) {
            pcout << "\nRelaxation time for biomass lattices does not satisfy a necessary stability condition for the BGK operator.\n";
            pcout << "tau must be > 0.5, but bioTAUinPore = " << bioTAUinPore[iM] << ". Consider increasing biomass diffusivity.\n";
            pcout << "Terminating the simulation.\n";
            return -1;
        }
    }
    pcout << std::endl;

    if (Pe > thrd) {
        pcout << "Peclet Number (meanU) = " << PoreMeanU * charcs_length / refNu << ", Grid Peclet Number (maxU) = " << PoreMaxUx / refNu << std::endl;
    }
    pcout << "ade_dt = " << ade_dt <<  " s/ts" << std::endl;

    // vector of substrate lattices
    MultiBlockLattice2D<T,RXNDES> substrLattice(nx, ny, new AdvectionDiffusionBGKdynamics<T,RXNDES>(refOmega));
    std::vector< MultiBlockLattice2D<T,RXNDES> > vec_substr_lattices(num_of_substrates, substrLattice);
    for (plint iS = 0; iS < num_of_substrates; ++iS) {
        soluteDomainSetup( vec_substr_lattices[iS], createLocalAdvectionDiffusionBoundaryCondition2D<T,RXNDES>(), geometry,
                           substrOMEGAinbFilm[iS], substrOMEGAinPore[iS], pore_dynamics, bounce_back, no_dynamics, bio_dynamics,
                           vec_c0[iS], vec_left_btype[iS], vec_right_btype[iS], vec_left_bcondition[iS], vec_right_bcondition[iS] );
    }

    // vector of biomass lattices
    MultiBlockLattice2D<T,RXNDES> initbFilmLattice(nx, ny, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.));
    MultiBlockLattice2D<T,RXNDES> copybFilmLattice(nx, ny, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.));
    MultiBlockLattice2D<T,RXNDES> initbFreeLattice(nx, ny, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.));
    MultiBlockLattice2D<T,RXNDES> copybFreeLattice(nx, ny, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.));
    std::vector< MultiBlockLattice2D<T,RXNDES> > vec_bFilm_lattices(bfilm_count, initbFilmLattice);
    std::vector< MultiBlockLattice2D<T,RXNDES> > vec_bFcopy_lattices(bfilm_count, copybFilmLattice);
    std::vector< MultiBlockLattice2D<T,RXNDES> > vec_bFree_lattices(bfree_count, initbFreeLattice);
    std::vector< MultiBlockLattice2D<T,RXNDES> > vec_bPcopy_lattices(bfree_count, copybFreeLattice);

    plint tmpIT0=0, tmpIT1=0;
    std::vector<plint> loctrack;
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        if (bmass_type[iM]==1) {
            bmassDomainSetup( vec_bFilm_lattices[tmpIT0], createLocalAdvectionDiffusionBoundaryCondition2D<T,RXNDES>(), geometry, bioOMEGAinPore[iM], bioOMEGAinbFilm[iM],
                              pore_dynamics, bounce_back, no_dynamics, bio_dynamics, bio_left_btype[iM], bio_right_btype[iM], bio_left_bcondition[iM], bio_right_bcondition[iM] );
            bmassDomainSetup( vec_bFcopy_lattices[tmpIT0], createLocalAdvectionDiffusionBoundaryCondition2D<T,RXNDES>(), geometry, 0., 0.,
                              pore_dynamics, bounce_back, no_dynamics, bio_dynamics, bio_left_btype[iM], bio_right_btype[iM], bio_left_bcondition[iM], bio_right_bcondition[iM] );
            loctrack.push_back(tmpIT0); ++tmpIT0;
        }
        else {
            if (solver_type[iM]==3) { // lb diffusion
                soluteDomainSetup( vec_bFree_lattices[tmpIT1], createLocalAdvectionDiffusionBoundaryCondition2D<T,RXNDES>(), geometry, bioOMEGAinbFilm[iM], bioOMEGAinPore[iM],
                       pore_dynamics, bounce_back, no_dynamics, bio_dynamics, vec_b0_free[tmpIT1], bio_left_btype[iM], bio_right_btype[iM], bio_left_bcondition[iM], bio_right_bcondition[iM] );
                bmassDomainSetup( vec_bPcopy_lattices[tmpIT1], createLocalAdvectionDiffusionBoundaryCondition2D<T,RXNDES>(), geometry, 0., 0.,
                                  pore_dynamics, bounce_back, no_dynamics, bio_dynamics, bio_left_btype[iM], bio_right_btype[iM], bio_left_bcondition[iM], bio_right_bcondition[iM] );
            }
            else if (solver_type[iM]==1) { // finite difference
                pcout <<  "ERROR: Finite difference for solute and biomass diffusion is not implemented yet." << std::endl;
                return -1;
                bmassDomainSetup( vec_bFree_lattices[tmpIT1], createLocalAdvectionDiffusionBoundaryCondition2D<T,RXNDES>(), geometry, bioOMEGAinPore[iM], bioOMEGAinbFilm[iM],
                                  pore_dynamics, bounce_back, no_dynamics, bio_dynamics, bio_left_btype[iM], bio_right_btype[iM], bio_left_bcondition[iM], bio_right_bcondition[iM] );
                bmassDomainSetup( vec_bPcopy_lattices[tmpIT1], createLocalAdvectionDiffusionBoundaryCondition2D<T,RXNDES>(), geometry, 0., 0.,
                                  pore_dynamics, bounce_back, no_dynamics, bio_dynamics, bio_left_btype[iM], bio_right_btype[iM], bio_left_bcondition[iM], bio_right_bcondition[iM] );
            }
            loctrack.push_back(tmpIT1); ++tmpIT1;
        }
    }
    MultiBlockLattice2D<T,RXNDES> totalbFilmLattice(nx, ny, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.));
    bmassDomainSetup( totalbFilmLattice, createLocalAdvectionDiffusionBoundaryCondition2D<T,RXNDES>(), geometry, bioOMEGAinPore[0], bioOMEGAinbFilm[0],
                      pore_dynamics, bounce_back, no_dynamics, bio_dynamics, bio_left_btype[0], bio_right_btype[0], bio_left_bcondition[0], bio_right_bcondition[0] );
    bmassDomainSetup( copybFilmLattice, createLocalAdvectionDiffusionBoundaryCondition2D<T,RXNDES>(), geometry, 0., 0.,
                      pore_dynamics, bounce_back, no_dynamics, bio_dynamics, bio_left_btype[0], bio_right_btype[0], bio_left_bcondition[0], bio_right_bcondition[0] );
    // MultiBlockLattice2D<T,RXNDES> totalbFreeLattice(nx, ny, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.));
    // bmassDomainSetup( totalbFreeLattice, createLocalAdvectionDiffusionBoundaryCondition2D<T,RXNDES>(), geometry, bioOMEGAinPore[0], bioOMEGAinbFilm[0],
    //                   pore_dynamics, bounce_back, no_dynamics, bio_dynamics, bio_left_btype[0], bio_right_btype[0], bio_left_bcondition[0], bio_right_bcondition[0] );

    // define initial biomass for biofilm lattices
    for (plint iM = 0; iM < bfilm_count; ++iM) {
        applyProcessingFunctional(new initializeScalarLattice<T,RXNDES,int> (vec_b0_film[iM], bio_dynamics[iM]), vec_bFilm_lattices[iM].getBoundingBox(), vec_bFilm_lattices[iM], geometry);
        std::vector<T> vec_b1(vec_b0_film[iM].size(),0.);
        applyProcessingFunctional(new initializeScalarLattice<T,RXNDES,int> (vec_b1, bio_dynamics[iM]), vec_bFcopy_lattices[iM].getBoundingBox(), vec_bFcopy_lattices[iM], geometry);
        initTotalbFilmLatticeDensity( vec_bFilm_lattices[iM], totalbFilmLattice);
    }

    // mask and distance lattices storing material numbers and distance from solid surface
    MultiBlockLattice2D<T,RXNDES> maskLattice(nx, ny, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.));
    MultiBlockLattice2D<T,RXNDES> ageLattice(nx, ny, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.));
    MultiBlockLattice2D<T,RXNDES> distLattice(nx, ny, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.));
    // copy geometry material numbers to maskLattice
    defineMaskLatticeDynamics( totalbFilmLattice, maskLattice, thrd_bFilmFrac );
    applyProcessingFunctional(new CopyGeometryScalar2maskLattice2D<T,RXNDES,int>(bio_dynamics), maskLattice.getBoundingBox(), maskLattice, geometry);
    applyProcessingFunctional(new CopyGeometryScalar2ageLattice2D<T,RXNDES,int>(), ageLattice.getBoundingBox(), ageLattice, ageDomain);
    applyProcessingFunctional(new CopyGeometryScalar2distLattice2D<T,RXNDES,int>(), distLattice.getBoundingBox(), distLattice, distanceDomain);

    MultiScalarField2D<int> latticeX(nx,ny);
    MultiScalarField2D<int> latticeY(nx,ny);
    applyProcessingFunctional(new latticeXY2D<T,RXNDES,int> (0), maskLattice.getBoundingBox(), maskLattice, latticeX);
    applyProcessingFunctional(new latticeXY2D<T,RXNDES,int> (1), maskLattice.getBoundingBox(), maskLattice, latticeY);
    if (track_performance == 0) {
        saveGeometry("latticeX", latticeX);
        saveGeometry("latticeY", latticeY);
        pcout << "saving subdomain iX and iY coordinates.\n";
    }

    // pointer vector of substrate lattices
    std::vector< MultiBlockLattice2D<T, RXNDES>* > substrate_lattices;
    for (plint iS = 0; iS < num_of_substrates; ++iS) {
        substrate_lattices.push_back( &vec_substr_lattices[iS] );
    }
    substrate_lattices.push_back( &maskLattice );

    std::vector< MultiBlockLattice2D<T, RXNDES>* > planktonic_lattices;
    for (size_t iP = 0; iP < vec_bFree_lattices.size(); ++iP) {
        planktonic_lattices.push_back( &vec_bFree_lattices[iP] );
    }
    planktonic_lattices.push_back( &maskLattice );

    std::vector< MultiBlockLattice2D<T, RXNDES>* > ptr_glpk_lattices, ptr_kns_lattices, ptr_nn_lattices;
    for (plint iS = 0; iS < num_of_substrates; ++iS) {
        ptr_glpk_lattices.push_back( &vec_substr_lattices[iS] );
        ptr_kns_lattices.push_back( &vec_substr_lattices[iS] );
        ptr_nn_lattices.push_back( &vec_substr_lattices[iS] );
    }
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        if (bmass_type[iM]==1) {
            if (reaction_type[iM]==0) { ptr_glpk_lattices.push_back( &vec_bFilm_lattices[loctrack[iM]] ); }
            else if (reaction_type[iM]==1) { ptr_kns_lattices.push_back( &vec_bFilm_lattices[loctrack[iM]] ); }
            else if (reaction_type[iM]==2) { ptr_nn_lattices.push_back( &vec_bFilm_lattices[loctrack[iM]] ); }
            else if (reaction_type[iM]==3) { ptr_glpk_lattices.push_back( &vec_bFilm_lattices[loctrack[iM]] ); ptr_kns_lattices.push_back( &vec_bFilm_lattices[loctrack[iM]] ); }
            else if (reaction_type[iM]==4) { ptr_nn_lattices.push_back( &vec_bFilm_lattices[loctrack[iM]] ); ptr_kns_lattices.push_back( &vec_bFilm_lattices[loctrack[iM]] ); }
        }
        else {
            if (reaction_type[iM]==0) { ptr_glpk_lattices.push_back( &vec_bFree_lattices[loctrack[iM]] ); }
            else if (reaction_type[iM]==1) { ptr_kns_lattices.push_back( &vec_bFree_lattices[loctrack[iM]] ); }
            else if (reaction_type[iM]==2) { ptr_nn_lattices.push_back( &vec_bFree_lattices[loctrack[iM]] ); }
            else if (reaction_type[iM]==3) { ptr_glpk_lattices.push_back( &vec_bFree_lattices[loctrack[iM]] ); ptr_kns_lattices.push_back( &vec_bFree_lattices[loctrack[iM]] ); }
            else if (reaction_type[iM]==4) { ptr_nn_lattices.push_back( &vec_bFree_lattices[loctrack[iM]] ); ptr_kns_lattices.push_back( &vec_bFree_lattices[loctrack[iM]] ); }
        }
    }
    ptr_glpk_lattices.push_back( &maskLattice );
    ptr_kns_lattices.push_back( &maskLattice );
    ptr_nn_lattices.push_back( &maskLattice );

    std::vector< MultiBlockLattice2D<T, RXNDES>* > ptr_ca_lattices;
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        if (solver_type[iM]==2) {
            if (bmass_type[iM]==1) { ptr_ca_lattices.push_back( &vec_bFilm_lattices[loctrack[iM]] ); }
            else { pcout << "ERROR: Cellular Automata can be applied only to sessile biomass. Terminating the simulation.\n"; return -1; }
        }
    }
    for (plint iM = 0; iM < num_of_microbes; ++iM) { if (solver_type[iM]==2) { ptr_ca_lattices.push_back( &vec_bFcopy_lattices[loctrack[iM]] ); } }
    ptr_ca_lattices.push_back( &totalbFilmLattice );
    ptr_ca_lattices.push_back( &maskLattice );
    ptr_ca_lattices.push_back( &ageLattice );
    plint caLlen = ptr_ca_lattices.size();
    if (2*ca_count+3 != caLlen) { pcout << "The length of ptr_ca_lattices = " << caLlen << " is not what it is supposed to be ("<< 2*ca_count+2 << "). Terminating the simulation.\n"; return -1; }

    std::vector< MultiBlockLattice2D<T, RXNDES>* > ptr_fd_lattices;
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        if (solver_type[iM]==1) {
            if (bmass_type[iM]==1) { ptr_fd_lattices.push_back( &vec_bFilm_lattices[loctrack[iM]] ); }
            else { ptr_fd_lattices.push_back( &vec_bFree_lattices[loctrack[iM]] ); }
        }
    }
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        if (solver_type[iM]==1) {
            if (bmass_type[iM]==1) { ptr_fd_lattices.push_back( &vec_bFcopy_lattices[loctrack[iM]] ); }
            else { ptr_fd_lattices.push_back( &vec_bPcopy_lattices[loctrack[iM]] ); }
        }
    }
    ptr_fd_lattices.push_back( &maskLattice );
    plint fdLlen = ptr_fd_lattices.size();
    if (2*fd_count+1 != fdLlen) { pcout << "The length of ptr_fd_lattices = " << fdLlen << "  is not what it is supposed to be (" << 2*fd_count+1 << "). Terminating the simulation.\n"; return -1; }

    std::vector< MultiBlockLattice2D<T, RXNDES>* > ageNdistance_lattices;
    ageNdistance_lattices.push_back( &ageLattice );
    ageNdistance_lattices.push_back( &distLattice );
    ageNdistance_lattices.push_back( &totalbFilmLattice );

    pcout << "\nbmass_type = ";
    for (size_t iT=0; iT<bmass_type.size(); ++iT) {
        if (bmass_type[iT]==1) { pcout << "biofilm "; }
        else { pcout << "planktonic "; }
    }
    pcout << std::endl;
    pcout << "solver_type = ";
    for (size_t iT=0; iT<solver_type.size(); ++iT) {
        if (solver_type[iT]==1) { pcout << "fd "; }
        else if (solver_type[iT]==2) { pcout << "ca "; }
        else { pcout << "lbm "; }
    }
    pcout << std::endl;

    // Couple the two physics (NS and RXN)
    if (Pe > thrd) {
        for (plint iS = 0; iS < num_of_substrates; ++iS) {
            latticeToPassiveAdvDiff(nsLattice, vec_substr_lattices[iS], vec_substr_lattices[iS].getBoundingBox());
        }
        // if biomass transport is simulated through the LB method
        tmpIT0=0;
        for (plint iM = 0; iM < num_of_substrates; ++iM) {
            if (solver_type[iM] == 3) {
                latticeToPassiveAdvDiff(nsLattice, vec_bFree_lattices[tmpIT0], vec_bFree_lattices[tmpIT0].getBoundingBox());
                ++tmpIT0;
            }
        }
        if (tmpIT0!=bfree_count) { pcout << "SOMETHING IS WRONG.\n" << std::endl; }
        if (track_performance == 0) {
            pcout << "Stabilizing the ADE lattices after coupling the NS lattice...\n";
            for (plint iT=0; iT<10000; ++iT) {
                for (plint iS = 0; iS < num_of_substrates; ++iS) {
                    vec_substr_lattices[iS].collideAndStream();
                }
                for (size_t iM = 0; iM < vec_bFree_lattices.size(); ++iM) {
                    vec_bFree_lattices[iM].collideAndStream();
                }
            }
            for (plint iS = 0; iS < num_of_substrates; ++iS) { applyProcessingFunctional(new stabilizeADElattice<T,RXNDES,int> (vec_c0[iS], pore_dynamics, bio_dynamics), vec_substr_lattices[iS].getBoundingBox(), vec_substr_lattices[iS], geometry); }
            for (size_t iM = 0; iM < vec_bFree_lattices.size(); ++iM) { applyProcessingFunctional(new stabilizeADElattice<T,RXNDES,int> (vec_b0_free[iM], pore_dynamics, bio_dynamics), vec_bFree_lattices[iM].getBoundingBox(), vec_bFree_lattices[iM], geometry); }
        }
    }

/*  ================================= rxn Lattice main loop  =================================  */
    iT = 0;
    if ( read_ADE_file==1) {
        if (ade_rerun_iT0>0 ) {
            pcout << "read binary ADE files" << std::endl;
            for (plint iS = 0; iS < num_of_substrates; ++iS) { loadBinaryBlock(vec_substr_lattices[iS], str_outputDir+ade_filename+"_"+std::to_string(iS)); }
            tmpIT0=0; tmpIT1=0;
            for (plint iM = 0; iM < num_of_microbes; ++iM) {
                if (bmass_type[iM]==1) { loadBinaryBlock(vec_bFilm_lattices[tmpIT0], str_outputDir+bio_filename+"_"+std::to_string(iM)); ++tmpIT0; }
                else { loadBinaryBlock(vec_bFree_lattices[tmpIT1], str_outputDir+bio_filename+"_"+std::to_string(iM)); ++tmpIT1; }
            }
            iT = ade_rerun_iT0;
            pcout << "ADE binary files successfully loaded" << std::endl;
        }
        else {
            pcout << "WARNING: number of input files should be equal to the number of substrates and microbes" << std::endl;
            return -1;
        }
    }
    T fbatime = 0, catime = 0, anntime = 0, adetime = 0, knstime = 0;

    auto t_start_ade = std::chrono::high_resolution_clock::now();
    auto t_start_fba = t_start_ade;
    auto t_start_ann = t_start_ade;
    auto t_start_ca = t_start_ade;
    auto t_start_kns = t_start_ade;
    pcout << "\n===================== LBM ADE simulation begins =====================\n\n";
    util::ValueTracer<T> ns_convg2(1.0,1000.0,ns_converge_iT2);
    plint old_totMask = util::roundToInt( computeAverage(*computeDensity(maskLattice))*nx*ny );
    bool ns_saturate=0, percolationFlag=0;

    for (; iT < ade_maxiTer; ++iT) {
        // ========================= save VTI files ========================= //
        if ( ade_VTI_iTer > 0 && iT % ade_VTI_iTer == 0 ) {
            pcout << "Iteration = " << iT << "; current_simulation_time = " << iT*ade_dt << " seconds" << std::endl;
            if ( track_performance == 0 ) {
                for (plint iS = 0; iS < num_of_substrates; ++iS) {
                    if (mm_count>0) { if (vec_fixC[iS] == 0) { writeAdvVTI( vec_substr_lattices[iS], iT, ade_filename+std::to_string(iS)+"_" ); } }
                    else { writeAdvVTI( vec_substr_lattices[iS], iT, ade_filename+std::to_string(iS)+"_" ); }
                }
                tmpIT0=0;tmpIT1=0;
                for (plint iM = 0; iM < num_of_microbes; ++iM) {
                    if (bmass_type[iM]==1) { writeAdvVTI( vec_bFilm_lattices[tmpIT0], iT, bio_filename+std::to_string(iM)+"_" ); ++tmpIT0; }
                    else { writeAdvVTI( vec_bFree_lattices[tmpIT1], iT, bio_filename+std::to_string(iM)+"_" ); ++tmpIT1; }
                }
                if (Pe > thrd) { writeNsVTI(nsLattice,iT,"nsLattice_"); }
                writeAdvVTI( maskLattice, iT, mask_filename+"_" );
                writeAdvVTI( ageLattice, iT, "ageLattice_" );
                pcout << "Writing ADE VTI files... \n";
            }
            T itvtime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now()-t_start_ade).count();
            adetime += itvtime;
            pcout << "(Wall-clock) Time elapsed: " << itvtime/1000 << " seconds." << std::endl;
            t_start_ade = std::chrono::high_resolution_clock::now();
        }
        // ===================== save checkpoint files ====================== //
        if (ade_CHK_iTer > 0 && iT % ade_CHK_iTer == 0 && iT > 0 && track_performance == 0) {
            tmpIT0=0;tmpIT1=0;
            for (plint iS = 0; iS < num_of_substrates; ++iS) {
                if (vec_fixC.size() > 0) { if (vec_fixC[iS] == 0) { saveBinaryBlock( vec_substr_lattices[iS], str_outputDir + ade_filename + std::to_string(iS) + "_" + std::to_string(iT) + ".chk"); } }
                else { saveBinaryBlock( vec_substr_lattices[iS], str_outputDir + ade_filename + std::to_string(iS) + "_" + std::to_string(iT) + ".chk"); }
            }
            for (plint iM = 0; iM < num_of_microbes; ++iM) {
                if (bmass_type[iM]==1) { saveBinaryBlock( vec_bFilm_lattices[tmpIT0], str_outputDir + bio_filename + std::to_string(iM) + "_" + std::to_string(iT) + ".chk"); ++tmpIT0; }
                else { saveBinaryBlock( vec_bFree_lattices[tmpIT1], str_outputDir + bio_filename + std::to_string(iM) + "_" + std::to_string(iT) + ".chk"); ++tmpIT1; }
            }
            saveBinaryBlock( maskLattice, str_outputDir + mask_filename + "_" + std::to_string(iT) + ".chk");
            pcout << "Writing checkpoint files... \n";
        }

        // =================== solute and bFree lattice collision ==================== //
        for (plint iS = 0; iS < num_of_substrates; ++iS) {
            if (mm_count>0) { if (vec_fixC[iS] == 0 && vec_fixLB[iS] == 0) { vec_substr_lattices[iS].collide(); } }
            else { vec_substr_lattices[iS].collide(); }
        }
        if (lb_count > 0) { // biomass LB diffusion (collision)
            for (plint iM = 0; iM < num_of_microbes; ++iM) {
                if (solver_type[iM]==3) {
                    if (bmass_type[iM]==1) {vec_bFilm_lattices[loctrack[iM]].collide(); }
                    else {vec_bFree_lattices[loctrack[iM]].collide(); }
                }
            }
        }

        // ============================ run reaction term ============================ //
        if (mm_count > 0) {
            if (track_performance == 1) { t_start_fba = std::chrono::high_resolution_clock::now(); }
            applyProcessingFunctional(new run_FBA<T,RXNDES> (nx, num_of_substrates, mm_count, ade_dt, max_bMassRho, vec_mu_glpk, vec_maxUptake, vec_maxRelease, vec_fixLB, vec_fixC,
                                                            vec_Kc_glpk, vec_EX_loc, vec_nrow, vec_ncol, vec_lp, vec_method, vec_isMIP, vec_sParam, vec_iParam, vec_objLoc, vec_Vmax_glpk),
                                                            vec_substr_lattices[0].getBoundingBox(), ptr_glpk_lattices);
            if (track_performance == 1) { fbatime+=std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now()-t_start_fba).count(); }
        }
        if (nn_count > 0) {
            if (track_performance == 1) { t_start_ann = std::chrono::high_resolution_clock::now(); }

            applyProcessingFunctional(new run_NN<T,RXNDES> (nx, num_of_substrates, nn_count, no_dynamics, bounce_back, ade_dt, reaction_type, vec_Vmax_nn, vec_Kc_nn, vec_mu_nn, 1),
                                                                 vec_substr_lattices[0].getBoundingBox(), ptr_nn_lattices);

            if (track_performance == 1) { anntime+=std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now()-t_start_ann).count(); }
        }
        if (kns_count > 0) {
            if (track_performance == 1) { t_start_kns = std::chrono::high_resolution_clock::now(); }
            applyProcessingFunctional(new run_kinetics<T,RXNDES> (nx, num_of_substrates, kns_count, ade_dt, vec_Kc_kns, vec_mu_kns, no_dynamics, bounce_back),
                                                                 vec_substr_lattices[0].getBoundingBox(), ptr_kns_lattices);
            if (track_performance == 1) { knstime+=std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now()-t_start_kns).count(); }
        }
        // if (py_count > 0) {
        //     applyProcessingFunctional(new runFBA_cobrapy<T,RXNDES> (iT,nx, num_of_substrates, num_of_microbes, ade_dt, 10., vec_mu, vec_maxUptake, vec_maxRelease, vec_fixLB, vec_fixC, vec_Kc, vec_EX_loc, vec_model, pyFileName),
        //                                                          vec_substr_lattices[0].getBoundingBox(), ptr_py_lattices);
        // }

        // ======================= biomass expansion ======================= //
        if (ca_count > 0) { // cellular automata
            applyProcessingFunctional(new updateLocalMaskNtotalLattices2D<T,RXNDES> (caLlen, bounce_back, no_dynamics, bio_dynamics, pore_dynamics, thrd_bFilmFrac, max_bMassRho), vec_bFilm_lattices[0].getBoundingBox(), ptr_ca_lattices);
            // T sumT0 = computeSum (*computeDensity(totalbFilmLattice));
            T globalBmax = computeMax( *computeDensity(totalbFilmLattice) );
            if (std::isnan(globalBmax) == 1) {
                pcout << "ERROR: biomass density goes wrong.\n";
                return -1;
            }
            plint whilecount=0;
            if (track_performance == 1) { t_start_ca = std::chrono::high_resolution_clock::now(); }
            while ( globalBmax-max_bMassRho > 1e-12 ) {
                for (plint iM=0; iM<bfilm_count; ++iM) { vec_bFcopy_lattices[iM]=copybFilmLattice; }
                if (halfflag == 0) { applyProcessingFunctional(new pushExcessBiomass2D<T,RXNDES> (max_bMassRho, nx, ny, 1, caLlen, no_dynamics, bounce_back, pore_dynamics), vec_bFilm_lattices[0].getBoundingBox(), ptr_ca_lattices); }
                else { applyProcessingFunctional(new halfPushExcessBiomass2D<T,RXNDES> (max_bMassRho, nx, ny, 1, caLlen, no_dynamics, bounce_back, pore_dynamics), vec_bFilm_lattices[0].getBoundingBox(), ptr_ca_lattices); }
                applyProcessingFunctional(new pullExcessBiomass2D<T,RXNDES> (nx, ny, 1, caLlen), vec_bFilm_lattices[0].getBoundingBox(), ptr_ca_lattices);
                applyProcessingFunctional(new updateLocalMaskNtotalLattices2D<T,RXNDES> (caLlen, bounce_back, no_dynamics, bio_dynamics, pore_dynamics, thrd_bFilmFrac, max_bMassRho), vec_bFilm_lattices[0].getBoundingBox(), ptr_ca_lattices);
                globalBmax = computeMax( *computeDensity(totalbFilmLattice) );
                if (whilecount%50 == 0) {
                    plint diff = 1;
                    plint whilecount1 = 0;
                    while (diff != 0) {
                        plint old_totAge = util::roundToInt( computeAverage(*computeDensity(ageLattice))*nx*ny );
                        applyProcessingFunctional(new updateAgeDistance2D<T,RXNDES>(max_bMassRho, nx, ny), ageLattice.getBoundingBox(), ageNdistance_lattices);
                        plint new_totAge = util::roundToInt( computeAverage(*computeDensity(ageLattice))*nx*ny );
                        diff = new_totAge-old_totAge;
                        ++whilecount1;
                        if (whilecount1 > 1000) {
                            pcout << "Iteration = " << iT << "; current_simulation_time = " << iT*ade_dt << " seconds" << std::endl;
                            if ( track_performance == 0 ) {
                                for (plint iS = 0; iS < num_of_substrates; ++iS) {
                                    if (mm_count>0) { if (vec_fixC[iS] == 0) { writeAdvVTI( vec_substr_lattices[iS], iT, ade_filename+std::to_string(iS)+"_" ); } }
                                    else { writeAdvVTI( vec_substr_lattices[iS], iT, ade_filename+std::to_string(iS)+"_" ); }
                                }
                                tmpIT0=0;tmpIT1=0;
                                for (plint iM = 0; iM < num_of_microbes; ++iM) {
                                    if (bmass_type[iM]==1) { writeAdvVTI( vec_bFilm_lattices[tmpIT0], iT, bio_filename+std::to_string(iM)+"_" ); ++tmpIT0; }
                                    else { writeAdvVTI( vec_bFree_lattices[tmpIT1], iT, bio_filename+std::to_string(iM)+"_" ); ++tmpIT1; }
                                }
                                if (Pe > thrd) { writeNsVTI(nsLattice,iT,"nsLattice_"); }
                                writeAdvVTI( maskLattice, iT, mask_filename+"_" );
                                writeAdvVTI( ageLattice, iT, "ageLattice_" );
                                pcout << "Writing ADE VTI files... \n";
                            }
                            pcout << "Stuck in the while loop 2. Terminating the simulation.\n";
                            exit(EXIT_FAILURE);
                        }
                    }
                }
                if (whilecount > 1000) {
                    pcout << "Iteration = " << iT << "; current_simulation_time = " << iT*ade_dt << " seconds" << std::endl;
                    if ( track_performance == 0 ) {
                        for (plint iS = 0; iS < num_of_substrates; ++iS) {
                            if (mm_count>0) { if (vec_fixC[iS] == 0) { writeAdvVTI( vec_substr_lattices[iS], iT, ade_filename+std::to_string(iS)+"_" ); } }
                            else { writeAdvVTI( vec_substr_lattices[iS], iT, ade_filename+std::to_string(iS)+"_" ); }
                        }
                        tmpIT0=0;tmpIT1=0;
                        for (plint iM = 0; iM < num_of_microbes; ++iM) {
                            if (bmass_type[iM]==1) { writeAdvVTI( vec_bFilm_lattices[tmpIT0], iT, bio_filename+std::to_string(iM)+"_" ); ++tmpIT0; }
                            else { writeAdvVTI( vec_bFree_lattices[tmpIT1], iT, bio_filename+std::to_string(iM)+"_" ); ++tmpIT1; }
                        }
                        if (Pe > thrd) { writeNsVTI(nsLattice,iT,"nsLattice_"); }
                        writeAdvVTI( maskLattice, iT, mask_filename+"_" );
                        writeAdvVTI( ageLattice, iT, "ageLattice_" );
                        pcout << "Writing ADE VTI files... \n";
                    }
                    pcout << "Stuck in the while loop. Terminating the simulation.\n";
                    exit(EXIT_FAILURE);
                }
                ++whilecount;
            }
            if (track_performance == 1) { catime+=std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now()-t_start_ca).count(); }
        }
        if (fd_count > 0) { // biomass diffusion: finite difference
            applyProcessingFunctional(new updateLocalMaskNtotalLattices2D<T,RXNDES> (fdLlen, bounce_back, no_dynamics, bio_dynamics, pore_dynamics, thrd_bFilmFrac, max_bMassRho), vec_bFilm_lattices[0].getBoundingBox(), ptr_fd_lattices);
            for (plint iM=0; iM<bfilm_count; ++iM) { vec_bFcopy_lattices[iM]=vec_bFilm_lattices[iM]; }
            for (plint iP=0; iP<bfree_count; ++iP) { vec_bPcopy_lattices[iP]=vec_bFree_lattices[iP]; }
            applyProcessingFunctional(new fdDiffusion2D<T,RXNDES> (nx, ny, fdLlen, 1, bioNUinPore[0]), vec_bFilm_lattices[0].getBoundingBox(), ptr_fd_lattices);
            applyProcessingFunctional(new updateLocalMaskNtotalLattices2D<T,RXNDES> (fdLlen, bounce_back, no_dynamics, bio_dynamics, pore_dynamics, thrd_bFilmFrac, max_bMassRho), vec_bFilm_lattices[0].getBoundingBox(), ptr_fd_lattices);
        }

        // ======================= update flow field and lattice dynamics =======================

        plint new_totMask = util::roundToInt( computeAverage(*computeDensity(maskLattice))*nx*ny );
        bool maskflag=0;
        if (std::abs(old_totMask-new_totMask)>0) { maskflag = 1; old_totMask = new_totMask; }
        if ( maskflag==1 ) {
            if (iT %ade_update_interval == 0) {
                if ( soluteDindex == 1 ) { applyProcessingFunctional(new updateSoluteDynamics2D<T,RXNDES> (num_of_substrates, bounce_back, no_dynamics, pore_dynamics, substrOMEGAinbFilm, substrOMEGAinPore), vec_substr_lattices[0].getBoundingBox(), substrate_lattices); }
                if ( bmassDindex == 1 ) { applyProcessingFunctional(new updateBiomassDynamics2D<T,RXNDES> ((plint) vec_bFree_lattices.size(), bounce_back, no_dynamics, pore_dynamics, bioOMEGAinbFilm, bioOMEGAinPore), vec_bFree_lattices[0].getBoundingBox(), planktonic_lattices); }
            }
            if (iT % ns_update_interval == 0) {
                applyProcessingFunctional(new updateAgeDistance2D<T,RXNDES>(max_bMassRho, nx, ny), ageLattice.getBoundingBox(), ageNdistance_lattices);
                //update ns lattice
                if (Pe > thrd && ns_saturate == 0) {
                    applyProcessingFunctional(new updateNsLatticesDynamics2D<T,NSDES,T,RXNDES> (nsLatticeOmega, vec_permRatio[0], pore_dynamics, no_dynamics, bounce_back), nsLattice.getBoundingBox(), nsLattice, maskLattice);
                    for (plint iT2 = 0; iT2 <ns_maxiTer_2; ++iT2) {
                        nsLattice.collideAndStream();
                        ns_convg2.takeValue(getStoredAverageEnergy(nsLattice),true);
                        if (ns_convg2.hasConverged()) { break; }
                        if ( iT2 == (ns_maxiTer_2-1) ) { ns_saturate = 1;}
                    }
                    if (ns_saturate == 1) {
                        T outletvel = computeAverage( *computeVelocityComponent(nsLattice, Box2D(nx-2,nx-2, 0,ny-1), 0) );
                        if (outletvel > thrd) { ns_saturate = 0; }
                        else { pcout << "\nThe simulation has reached a percolation limit. Terminating the simulation at iT = " << iT << ".\n"; percolationFlag = 1;}
                    }
                    for (plint iS = 0; iS < num_of_substrates; ++iS) { latticeToPassiveAdvDiff(nsLattice, vec_substr_lattices[iS], vec_substr_lattices[iS].getBoundingBox()); }
                    // if biomass transport is simulated through the LB method
                    if ( lb_count > 0 ) {
                        for (plint iM = 0; iM < num_of_microbes; ++iM) {
                            if (solver_type[iM]==3) {
                                if (bmass_type[iM]==1) { latticeToPassiveAdvDiff(nsLattice, vec_bFilm_lattices[loctrack[iM]], vec_bFilm_lattices[loctrack[iM]].getBoundingBox()); }
                                else { latticeToPassiveAdvDiff(nsLattice, vec_bFree_lattices[loctrack[iM]], vec_bFree_lattices[loctrack[iM]].getBoundingBox()); }
                            }
                        }
                    }
                }
            }
        }

        // solute lattice stream
        for (plint iS = 0; iS < num_of_substrates; ++iS) {
            if (mm_count>0) { if (vec_fixC[iS] == 0 && vec_fixLB[iS] == 0) { vec_substr_lattices[iS].stream(); } }
            else { vec_substr_lattices[iS].stream(); }
        }
        if (lb_count > 0) { // biomass LB diffusion (collision)
            for (plint iM = 0; iM < num_of_microbes; ++iM) {
                if (solver_type[iM]==3) {
                    if (bmass_type[iM]==1) {vec_bFilm_lattices[loctrack[iM]].stream(); }
                    else {vec_bFree_lattices[loctrack[iM]].stream(); }
                }
            }
        }
        if (percolationFlag == 1) {break;}
    }
    pcout << "End of simulation at iteration " << iT << std::endl << std::endl;

/*  ================================= Finalize the simulation  =================================  */

    T TET = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now()-t_start_0).count()/1000;
    pcout << "Total elapsed time: " << TET << " seconds, " << TET/60 << " minutes, and " << TET/3600 << " hours."<< std::endl;
    if (track_performance == 1) {
        adetime += std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now()-t_start_ade).count();
        nstime /= 1000; adetime /= 1000; fbatime /= 1000; catime /= 1000; anntime /= 1000; knstime /= 1000;
        pcout << "Total time consumed by NS: " << nstime << " seconds, " << nstime/60 << " minutes, and " << nstime/3600 << " hours."<< std::endl;
        pcout << "Total time consumed by ADE: " << adetime << " seconds, " << adetime/60 << " minutes, and " << adetime/3600 << " hours."<< std::endl;
        if (mm_count > 0) { pcout << "Total time consumed by FBA: " << fbatime << " seconds, " << fbatime/60 << " minutes, and " << fbatime/3600 << " hours."<< std::endl; }
        if (ca_count > 0) { pcout << "Total time consumed by CA: " << catime << " seconds, " << catime/60 << " minutes, and " << catime/3600 << " hours."<< std::endl; }
        if (nn_count > 0) { pcout << "Total time consumed by ANN: " << anntime << " seconds, " << anntime/60 << " minutes, and " << anntime/3600 << " hours."<< std::endl; }
        if (kns_count > 0) { pcout << "Total time consumed by KNS: " << knstime << " seconds, " << knstime/60 << " minutes, and " << knstime/3600 << " hours."<< std::endl; }
    }
    else {
        pcout << "Writing VTI and CHK files ..." << std::endl << std::endl;
        for (plint iS = 0; iS < num_of_substrates; ++iS) {
            writeAdvVTI( vec_substr_lattices[iS], iT, ade_filename+std::to_string(iS)+"_" );
            saveBinaryBlock( vec_substr_lattices[iS], str_outputDir+ade_filename+std::to_string(iS)+"_"+std::to_string(iT)+".chk");
        }
        tmpIT0=0;tmpIT1=0;
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (bmass_type[iM]==1) { writeAdvVTI( vec_bFilm_lattices[tmpIT0], iT, bio_filename+std::to_string(iM)+"_" ); ++tmpIT0; }
            else { writeAdvVTI( vec_bFree_lattices[tmpIT1], iT, bio_filename+std::to_string(iM)+"_" ); ++tmpIT1; }
        }
        writeAdvVTI( maskLattice, iT, mask_filename+"_" );
        saveBinaryBlock( maskLattice, str_outputDir+mask_filename+"_"+std::to_string(iT)+".chk");
        if (Pe > thrd) {
            writeNsVTI(nsLattice,iT,"nsLattice_");
            saveBinaryBlock(nsLattice, str_outputDir+ns_filename+".chk");
        }
    }

    // finalize glpk
    if (mm_count > 0) {
        pcout << "\nFinalize GLPK" << std::endl;
        // for (plint iT = 0; iT < num_of_microbes; ++iT) { glp_delete_prob(vec_lp[iT]); }
        /* this shouldn't be nessiary with glp_deleted_prob, but try it if we have weird behavior again... */
        glp_free_env();
    }

    pcout << "\nSimulation Finished!" << std::endl << std::endl;
    return 0;
}