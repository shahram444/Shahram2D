#include <climits> //INT_MAX
#include <cfloat> //DBL_MAX
#include "palabos2D.h"
#include "palabos2D.hh"

using namespace plb;

#define NIntP  21
#define NRealP 11
// control parameters not defined in glpk.h
#define LPX_K_PREPROCESS 401   /* preprocessing  */
#define LPX_K_RATIO_TEST 402   /* ratio test */

int set_glpk (glp_prob *lp, int sense, int nrow, int ncol, int nz, int &method, 
            double *c, int *rn, int *cn, double *a, double *b, char *ctype,
            bool *freeLB, double *lb, bool *freeUB, double *ub, int *vartype, int &isMIP, int lpsolver, glp_smcp &sParam, glp_iocp &iParam,
            int save_pb, char *save_filename, int *glpIntParam, int *IParam, double *glpRealParam, int *RParam)
{
  //Redirect standard output
  glp_term_hook (NULL, NULL); // default

  //-- Set the sense of optimization (the optimization direction flag)
  if (sense == 1) { // sense is the first argument of the glpk function
    glp_set_obj_dir (lp, GLP_MIN); // minimization
  }
  else {
    glp_set_obj_dir (lp, GLP_MAX); // maximization
  }

  //-- Define the number of unknowns and their domains.
  
  // column (structural) variables
  glp_add_cols (lp, ncol);
  for (int iT = 0; iT < ncol; ++iT) {
    //-- Define type of the structural variables
    if (!freeLB[iT] && !freeUB[iT]) {
      if ( lb[iT] == ub[iT] ) {
        glp_set_col_bnds (lp, iT+1, GLP_FX, lb[iT], ub[iT]); // GLP_FX: Fixed variable
      }
      else {
        glp_set_col_bnds (lp, iT+1, GLP_DB, lb[iT], ub[iT]); // GLP_DB: Double-bounded variable
      }
    }
    else {
      if (!freeLB[iT] && freeUB[iT]) {
        glp_set_col_bnds (lp, iT+1, GLP_LO, lb[iT], ub[iT]); // Variable with lower bound
      }
      else {
        if (freeLB[iT] && !freeUB[iT]) {
          glp_set_col_bnds (lp, iT+1, GLP_UP, lb[iT], ub[iT]); // Variable with upper bound
        }
        else {
          glp_set_col_bnds (lp, iT+1, GLP_FR, lb[iT], ub[iT]); // Free (unbounded) variable
        }
      }
    }
  
    // -- Set the objective coefficient of the corresponding
    // -- structural variable. No constant term is assumed.
    glp_set_obj_coef(lp, iT+1, c[iT]);

    if (isMIP) {
      glp_set_col_kind (lp, iT+1, vartype[iT]);
    }
  }

  // row (auxiliary) variables
  int typx = 0;
  glp_add_rows (lp, nrow);
  for (int iT = 0; iT < nrow; ++iT)
  {
  /*  If the i-th row has no lower bound (types F,U),
      the corrispondent parameter will be ignored.
      If the i-th row has no upper bound (types F,L),
      the corrispondent parameter will be ignored.
      If the i-th row is of S type,
      the i-th LB is used, but the i-th UB is ignored. */
    switch (ctype[iT])
    {
      // Free bound
      case 'F': typx = GLP_FR; break;
      // upper bound
      case 'U': typx = GLP_UP; break;
      // lower bound
      case 'L': typx = GLP_LO; break;
      // fixed constraint
      case 'S': typx = GLP_FX; break;
      // double-bounded variable
      case 'D': typx = GLP_DB; break;
    }

    glp_set_row_bnds (lp, iT+1, typx, b[iT], b[iT]);
    // glp_set_row_bnds (lp, iT+1, typx, 0, 0);

    // if (b[iT] != 0) {
    //   std::cout << b[iT] << " " << iT << std::endl;
    // }

  }

  // Load the vectorized constraint matrix a
  glp_load_matrix (lp, nz, rn, cn, a);

  // Save problem
  if ( save_pb == 1 ){
    if (glp_write_lp (lp, NULL, save_filename) != 0) {
      pcout << "glpkcpp: unable to write the problem."<< std::endl;
      exit(EXIT_FAILURE);
    }
  }
  else if ( save_pb == 2) {
    if (glp_write_mps (lp, GLP_MPS_DECK, NULL, save_filename) != 0) {
      pcout << "glpkcpp: unable to write the problem."<< std::endl;
      exit(EXIT_FAILURE);
    }
  }
  else if ( save_pb == 3) {
    if (glp_write_mps (lp, GLP_MPS_FILE, NULL, save_filename) != 0) {
      pcout << "glpkcpp: unable to write the problem."<< std::endl;
      exit(EXIT_FAILURE);
    }
  }
  else if ( save_pb == 4 ) {
    if (lpx_print_prob (lp, save_filename) != 0) {
      pcout << "glpkcpp: unable to write the problem."<< std::endl;
      exit(EXIT_FAILURE);
    }
  }

  //-- scale the problem data (if required) (see manual page 33)
  // In GLPK the scaling means a linear transformation applied to the constraint matrix to improve its numerical properties
  if (glpIntParam[1] && (! glpIntParam[16] || lpsolver != 1)) {
    switch ( glpIntParam[1] ) {
      case ( 1 ): glp_scale_prob( lp, GLP_SF_GM ); break; // geometric mean scaling
      case ( 2 ): glp_scale_prob( lp, GLP_SF_EQ ); break; // equilibration scaling
      case ( 3 ): glp_scale_prob( lp, GLP_SF_GM | GLP_SF_EQ ); break; 
      case ( 4 ): glp_scale_prob( lp, GLP_SF_2N ); break; // round scale factors to nearest power of two
      default :
        pcout << "glpkcpp: unrecognized scaling option" << std::endl;
        exit(EXIT_FAILURE);
    }
  }
  else {
    /* do nothing? or unscale?
        glp_unscale_prob( lp );
    */
  }

  //-- build advanced initial LP basis (if required)
  if (lpsolver == 1 && ! glpIntParam[16])
    glp_adv_basis (lp, 0);
  
  //-- set control parameters for simplex (lpsolver == 1) / exact (lpsolver == 3) method (see manual page 41)
  if (lpsolver == 1 || lpsolver == 3) {
    //remap of control parameters for simplex method
    sParam.msg_lev=glpIntParam[0];  // message level

    // simplex method: primal/dual
    switch ( glpIntParam[2] ) {
      case 0: sParam.meth=GLP_PRIMAL; break;
      case 1: sParam.meth=GLP_DUAL;   break;
      case 2: sParam.meth=GLP_DUALP;  break;
      default: 
        pcout << "glpkcpp: unrecognized primal/dual method."<< std::endl;
        exit(EXIT_FAILURE);
    }

    // pricing technique
    if (glpIntParam[3]==0) sParam.pricing=GLP_PT_STD;
    else sParam.pricing=GLP_PT_PSE;

    // ratio test    
    if (glpIntParam[20]==0) sParam.r_test = GLP_RT_STD;
    else sParam.r_test=GLP_RT_HAR;

    //tollerances
    sParam.tol_bnd=glpRealParam[1]; // primal feasible tollerance
    sParam.tol_dj=glpRealParam[2];  // dual feasible tollerance
    sParam.tol_piv=glpRealParam[3]; // pivot tollerance
    sParam.obj_ll=glpRealParam[4];  // lower limit
    sParam.obj_ul=glpRealParam[5];  // upper limit

    // iteration limit
    if (glpIntParam[5]==-1) sParam.it_lim=INT_MAX;
    else sParam.it_lim=glpIntParam[5];   

    // time limit
    if (glpRealParam[6]==-1) sParam.tm_lim=INT_MAX;
    else sParam.tm_lim=(int) glpRealParam[6]; 
    sParam.out_frq=glpIntParam[7];  // output frequency
    sParam.out_dly=(int) glpRealParam[7]; // output delay
    // presolver
    if (glpIntParam[16]) sParam.presolve=GLP_ON;
    else sParam.presolve=GLP_OFF;
  }
  else {
    for (int iT = 0; iT < NIntP; ++iT) {
      // skip assinging ratio test 
      if ( iT == 18 || iT == 20) {
        continue; 
      }
      // or
      lpx_set_int_parm (lp, IParam[iT], glpIntParam[iT]);
    }   

    for (int iT = 0; iT < NRealP; ++iT) {
      lpx_set_real_parm (lp, RParam[iT], glpRealParam[iT]);
    }
  }
  
  // set MIP params if MIP....
  if (isMIP) {
    method = 'I';
   
    switch (glpIntParam[0]) { //message level
      case 0:  iParam.msg_lev = GLP_MSG_OFF;   break;
      case 1:  iParam.msg_lev = GLP_MSG_ERR;   break;
      case 2:  iParam.msg_lev = GLP_MSG_ON;    break;
      case 3:  iParam.msg_lev = GLP_MSG_ALL;   break;
      default: pcout << "__glpk__: msg_lev bad param" << std::endl;
               exit(EXIT_FAILURE);
    }
    switch (glpIntParam[14]) { //branching param
      case 0:  iParam.br_tech = GLP_BR_FFV;    break;
      case 1:  iParam.br_tech = GLP_BR_LFV;    break;
      case 2:  iParam.br_tech = GLP_BR_MFV;    break;
      case 3:  iParam.br_tech = GLP_BR_DTH;    break;
      default: pcout << "__glpk__: branch bad param" << std::endl;
               exit(EXIT_FAILURE);
    }
    switch (glpIntParam[15]) { //backtracking heuristic
      case 0:  iParam.bt_tech = GLP_BT_DFS;    break;
      case 1:  iParam.bt_tech = GLP_BT_BFS;    break;
      case 2:  iParam.bt_tech = GLP_BT_BLB;    break;
      case 3:  iParam.bt_tech = GLP_BT_BPH;    break;
      default: pcout << "__glpk__: backtrack bad param" << std::endl;
               exit(EXIT_FAILURE);
    }

    if ( glpRealParam[8] > 0.0 && glpRealParam[8] < 1.0 ) {
      iParam.tol_int = glpRealParam[8];  // absolute tolorence
    }
    else {
      pcout << "__glpk__: tolint must be between 0 and 1" << std::endl;
      exit(EXIT_FAILURE);
    }

    iParam.tol_obj = glpRealParam[9];  // relative tolarence
    iParam.mip_gap = glpRealParam[10]; // realative gap tolerance

    // set time limit for mip
    if ( glpRealParam[6] < 0.0 || glpRealParam[6] > 1e6 ) {
      iParam.tm_lim = INT_MAX;
    }
    else {
      iParam.tm_lim = (int)(1000.0 * glpRealParam[6] );
    }

    // Choose Cutsets for mip
    // shut all cuts off, then start over....
    iParam.gmi_cuts = GLP_OFF; 
    iParam.mir_cuts = GLP_OFF; 
    iParam.cov_cuts = GLP_OFF; 
    iParam.clq_cuts = GLP_OFF;

    switch( glpIntParam[17] ) {
      case 0: break; 
      case 1: iParam.gmi_cuts = GLP_ON; break;
      case 2: iParam.mir_cuts = GLP_ON; break;
      case 3: iParam.cov_cuts = GLP_ON; break;
      case 4: iParam.clq_cuts = GLP_ON; break;
      case 5: iParam.clq_cuts = GLP_ON; 
              iParam.gmi_cuts = GLP_ON; 
              iParam.mir_cuts = GLP_ON;  
              iParam.cov_cuts = GLP_ON; 
              iParam.clq_cuts = GLP_ON; break;
      default: pcout << "__glpk__: cutset bad param" << std::endl;
               exit(EXIT_FAILURE);
    }

    switch( glpIntParam[18] ) { // pre-processing for mip
        case 0: iParam.pp_tech = GLP_PP_NONE; break;
        case 1: iParam.pp_tech = GLP_PP_ROOT; break;
        case 2: iParam.pp_tech = GLP_PP_ALL;  break;
        default: pcout << "__glpk__: pprocess bad param" << std::endl;
                 exit(EXIT_FAILURE);
    }

    if (glpIntParam[16])  iParam.presolve=GLP_ON; 
    else                  iParam.presolve=GLP_OFF; 

    if (glpIntParam[19])  iParam.binarize = GLP_ON; 
    else                  iParam.binarize = GLP_OFF;

  }
  else {
     /* Choose simplex method ('S') 
     or interior point method ('T') 
     or Exact method          ('E') 
     to solve the problem  */
    switch (lpsolver) {
      case 1: method = 'S'; break;
      case 2: method = 'T'; break;
      case 3: method = 'E'; break;
      default: 
            pcout << "__glpk__:  lpsolver != lpsovler" << std::endl;
            exit(EXIT_FAILURE);
    }
  }

  return 0;
}


int run_glpk (glp_prob *lp, int method, int isMIP, int nrow, int ncol, int lpsolver, glp_smcp sParam, glp_iocp iParam,
             std::vector<double> &xmin, double &fmin, int &status, std::vector<int> location, std::vector<double> lb, 
             std::vector<double> ub, std::vector<double> &lambda, std::vector<double> &redcosts, double &time, double &mem)

{
  clock_t t_start = clock();

  for (size_t iL = 0; iL < location.size(); ++iL) {
    if (location[iL] >= 0 ) {
      glp_set_col_bnds (lp, location[iL]+1, GLP_DB, lb[iL], ub[iL]); // GLP_DB: Double-bounded variable
    }
  }

  // now run the problem...
  int errnum;
  switch (method) {
    case 'I': errnum = glp_intopt( lp, &iParam ); // solve MIP problem with the branch-and-cut method
              errnum += 100; //this is to avoid ambiguity in the return codes.
              break;

    case 'S': errnum = glp_simplex(lp, &sParam); // solve LP problem with the primal or dual simplex method
              errnum += 100; //this is to avoid ambiguity in the return codes.
              break;

    case 'T': errnum = glp_interior(lp, NULL ); break; // solve LP problem with the interior-point method

    case 'E': errnum = glp_exact(lp, &sParam); break; // solve LP problem in exact arithmetic

    default:  /*xassert (method != method); */
        pcout << "__glpk__: method != method" << std::endl;
        exit(EXIT_FAILURE);
  }

  /* errnum assumes the following results:
     errnum = 0   <=> No errors
     errnum = 109 <=> Iteration limit exceeded.
     errnum = 2   <=> Numerical problems with basis matrix. */
  if (errnum == LPX_E_OK || errnum==100 || errnum ==109 || errnum == 108) {
    // Get status and object value
    if (isMIP) {
      status = glp_mip_status (lp);
      fmin = glp_mip_obj_val (lp);
    }
    else {
      if (lpsolver == 1 || lpsolver == 3) {
        status = glp_get_status (lp);
        fmin = glp_get_obj_val (lp);
      }
      else {
        status = glp_ipt_status (lp);
        fmin = glp_ipt_obj_val (lp);
      }
    }

    // Get optimal solution (if exists)
    if (isMIP) {
      for (int iT = 0; iT < ncol; ++iT) {
        xmin[iT] = glp_mip_col_val (lp, iT+1);
      }
    }
    else {
      /* Primal values */
      for (int iT = 0; iT < ncol; ++iT) {
        if (lpsolver == 1 || lpsolver == 3) {
          xmin[iT] = glp_get_col_prim (lp, iT+1);
        }
        else {
          xmin[iT] = glp_ipt_col_prim (lp, iT+1);
        }
      }
      /* Dual values */
      for (int iT = 0; iT < nrow; ++iT) {
        if (lpsolver == 1 || lpsolver == 3) {
          lambda[iT] = glp_get_row_dual (lp, iT+1);
        }
        else {
          lambda[iT] = glp_ipt_row_dual (lp, iT+1);
        }
      }
      /* Reduced costs */
      for (int iT = 0; iT < ncol; ++iT) {
        if (lpsolver == 1 || lpsolver == 3) {
          redcosts[iT] = glp_get_col_dual (lp, iT+1);
        }
        else {
          redcosts[iT] = glp_ipt_col_dual (lp, iT+1);
        }
      }
    }
    
    time = (clock () - t_start) / CLOCKS_PER_SEC;

    glp_long tpeak;
    glp_mem_usage(NULL, NULL, NULL, &tpeak);
    mem=(double)(4294967296.0 * tpeak.hi + tpeak.lo) / (1024);      

    return 0;
  }
  else {
    // printf("errnum is %d\n", errnum);
  }
  status = errnum;

  return errnum;
}

void initialize_glpk(int iter,glp_prob *lp, std::vector< std::vector<double> > S, std::vector<double> vec_b, std::vector<double> vec_c,
                    std::vector<double> vec_lb, std::vector<double> vec_ub, char *ctype, char *vtype,
                    int sense, int lpsolver, int save_pb, int &method, int &isMIP, glp_smcp &sParam, glp_iocp &iParam) {

  std::string siter = std::to_string(iter+1);
  if (iter == 0) siter += "st ";
  else if (iter == 1) siter += "nd ";
  else if (iter == 2) siter += "rd ";
  else siter += "th ";
  // pcout << siter << "calling glpk_initialize function ";

  // check if the glpk version is right
  std::string glpk_ver = "";
  for (size_t iT = 0; iT < sizeof(glp_version()); ++iT) {
    glpk_ver += glp_version()[iT];
  }
  // pcout << "(GLPK Version " << glpk_ver << ")." << std::endl;
  if (!glpk_ver.compare("4.39")) {
    pcout << "The code is written for the glpk version 4.39. Compatibility is not guaranteed" << std::endl;
  }

  // Integer/Real Param Defaluts
  // int glpIntParam[NIntP] = { 0,1,0,1,0,INT_MAX,INT_MAX,200,1,2,0,1,0,0,3,2, presolver, 0,2,0,1 };
  int glpIntParam[NIntP] = { 0,1,0,1,0,INT_MAX,INT_MAX,200,1,2,0,1,0,0,3,2,1,0,2,0,1 };
  double glpRealParam[NRealP] = {0.07,1e-7,1e-7,1e-10,-DBL_MAX,DBL_MAX,INT_MAX,0.0,1e-5,1e-7,0.0};

  // Integer/Real Param Names
  int IParam[NIntP] = {LPX_K_MSGLEV,LPX_K_SCALE,LPX_K_DUAL,LPX_K_PRICE,LPX_K_ROUND,LPX_K_ITLIM,LPX_K_ITCNT,
    LPX_K_OUTFRQ,LPX_K_MPSINFO,LPX_K_MPSOBJ,LPX_K_MPSORIG,LPX_K_MPSWIDE,LPX_K_MPSFREE,LPX_K_MPSSKIP,
    LPX_K_BRANCH,LPX_K_BTRACK,LPX_K_PRESOL,LPX_K_USECUTS,LPX_K_PREPROCESS,LPX_K_BINARIZE,LPX_K_RATIO_TEST};
  int RParam[NRealP] = {LPX_K_RELAX,LPX_K_TOLBND,LPX_K_TOLDJ,LPX_K_TOLPIV,LPX_K_OBJLL,LPX_K_OBJUL,LPX_K_TMLIM,
    LPX_K_OUTDLY,LPX_K_TOLINT,LPX_K_TOLOBJ,LPX_K_MIPGAP};

  const int nrow = S.size();    // number of substrates
  const int ncol = S[0].size(); // number of reactions

  double c[ncol], lb[ncol], ub[ncol];
  memset(  c, 0, ncol*sizeof(int) );
  memset( lb, 0, ncol*sizeof(int) );
  memset( ub, 0, ncol*sizeof(int) );

  for (int iT = 0; iT < ncol; ++iT) {
     c[iT] =  vec_c[iT];
    lb[iT] = vec_lb[iT];
    ub[iT] = vec_ub[iT];
  }

  int *rn, *cn;
  rn = (int *) calloc(nrow*ncol+1,sizeof(int));
  cn = (int *) calloc(nrow*ncol+1,sizeof(int));

  double  *a;
  a = (double *) calloc(nrow*ncol+1,sizeof(double));

  int nz = 0;
  for (int iT0 = 0; iT0 < nrow; ++iT0) {
      for (int iT1 = 0; iT1 < ncol; ++iT1) {
        if (S[iT0][iT1] != 0) {
          ++nz;
          rn[nz] = iT0 + 1;
          cn[nz] = iT1 + 1;
           a[nz] = S[iT0][iT1]; // a matrix containing the constraints coefficients.
      }
    }
  }

  double  *b;
  b = (double *) calloc(nrow,sizeof(double)); // calloc initializes the memory to 0
  for (int iT0 = 0; iT0 < nrow; ++iT0) {
    b[iT0] = vec_b[iT0];
  }

  //-- freeLB/UB arguments, default: Free
  bool *freeLB, *freeUB;
  freeLB = (bool *) calloc(ncol,sizeof(bool));
  freeUB = (bool *) calloc(ncol,sizeof(bool));
  for (int iT = 0; iT < ncol; iT++) {
    if ( isinf(lb[iT]) ) {
      freeLB[iT] = 1;
    }
    else {
      freeLB[iT] = 0;
    }
    if ( isinf(ub[iT]) ) {
      freeUB[iT] = 1;
    }
    else {
      freeUB[iT] = 0;
    }
  }

  //-- vartype arguments, default: Free
  int *vartype;
  isMIP = 0;
  vartype = (int *) calloc(ncol+1,sizeof(int));
  for (int iT = 0; iT < (ncol+1) ; ++iT) {
    switch ( vtype[iT] ) {
      case 'I': vartype[iT] = GLP_IV; isMIP = 1; break;
      case 'B': vartype[iT] = GLP_BV; isMIP = 1; break;
      default : vartype[iT] = GLP_CV;
    }
  }

  //-- Save option
  char save_filename[15];
  strcpy(save_filename,"glpk_output");
  if (save_pb > 0) {
    char save_filetype[5]; // .mps; .txt; .lp
    if (save_pb == 1) {
      strcpy(save_filetype,".lp");
      if (isMIP == 1) {
        pcout << "WARNING: save file type incompatible with the input vtype" << std::endl;
      }
    }
    else if (save_pb == 2 || save_pb == 3) {
      strcpy(save_filetype,".mps");
      if (isMIP == 0) {
        pcout << "WARNING: save file type incompatible with the input vtype" << std::endl;
      }
    }
    else if (save_pb == 4) {
      strcpy(save_filetype,".mps");
    }
    else {
      pcout << "ERROR: unspecified save_pb" << std::endl;
      exit(EXIT_FAILURE);
    }
    strcat(save_filename,save_filetype);
  }

  // glp_iocp iParam;
  // glp_smcp sParam;
  glp_init_iocp(&iParam);
  glp_init_smcp(&sParam);

  int errchk = set_glpk (lp, sense, nrow, ncol, nz, method, c, rn, cn, a, b,
          ctype, freeLB, lb, freeUB, ub, vartype, isMIP, lpsolver, sParam, iParam,
          save_pb, save_filename, glpIntParam, IParam, glpRealParam, RParam);

  if (!errchk) {
    pcout << siter << "metabolic model has been initialized successfully." << std::endl;
  }
  else {
    pcout << "ERROR in set_glpk" << std::endl;
    exit(EXIT_FAILURE);
  }
}
