// nonlocal biomass reaction
template<typename T, template<typename U> class Descriptor>
class localReactionLattices2D : public LatticeBoxProcessingFunctional2D<T,Descriptor>
{
public:
    localReactionLattices2D(T dt_) : dt(dt_)
    {}
    // lattices[0] = solute 1 concentration
    // lattices[1] = surface-associated biomass 1
    // lattices[2] = mask field lattice
    virtual void process(Box2D domain, std::vector<BlockLattice2D<T,Descriptor>*> lattices) {
        Dot2D offset_01 = computeRelativeDisplacement(*lattices[0],*lattices[1]);
        Dot2D offset_02 = computeRelativeDisplacement(*lattices[0],*lattices[2]);
        for (plint iX0=domain.x0; iX0<=domain.x1; ++iX0) {
            plint iX2 = iX0 + offset_02.x;
            for (plint iY0=domain.y0; iY0<=domain.y1; ++iY0) {
                plint iY2 = iY0 + offset_02.y;
                plint mask = util::roundToInt( lattices[2]->get(iX2,iY2).computeDensity() );
                if ( (mask > 1 && mask < 5) || mask > 6 ) { // biomass
                    Array<T,5> g;
                    plint iX1 = iX0 + offset_01.x; plint iY1 = iY0 + offset_01.y;

                    T c0 = lattices[0]->get(iX0,iY0).computeDensity(); // (mM)
                    T B2 = lattices[1]->get(iX1,iY1).computeDensity(); // sessile (gdw/L)

                    call_biomassGrowth(c0, B2, dt);

                    lattices[0]->get(iX0,iY0).getPopulations(g);
                    g[0]+=(T) c0/3; g[1]+=(T) c0/6; g[2]+=(T) c0/6; g[3]+=(T) c0/6; g[4]+=(T) c0/6;
                    lattices[0]->get(iX0,iY0).setPopulations(g);

                    lattices[1]->get(iX1,iY1).getPopulations(g);
                    g[0]+=(T) B2/3; g[1]+=(T) B2/6; g[2]+=(T) B2/6; g[3]+=(T) B2/6; g[4]+=(T) B2/6;
                    lattices[1]->get(iX1,iY1).setPopulations(g);
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    virtual localReactionLattices2D<T,Descriptor>* clone() const {
        return new localReactionLattices2D<T,Descriptor>(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
        modified[1] = modif::staticVariables;
        modified[2] = modif::nothing;
    }
private:
    T dt;
};

// nonlocal biomass reaction
template<typename T, template<typename U> class Descriptor>
class updateTOTbMassLattices2D : public LatticeBoxProcessingFunctional2D<T,Descriptor>
{
public:
    updateTOTbMassLattices2D(plint length_) : length(length_)
    {}
    // lattices[0~(#ofbM-1)] = original biomass lattices
    // lattices[#ofbM~(len-3)] = copy biomass lattices
    // lattices[len-2] = total biomass lattice
    // lattices[len-1] = mask lattice
    virtual void process(Box2D domain, std::vector<BlockLattice2D<T,Descriptor>*> lattices) {
        std::vector<Dot2D> vec_offset;
        for (plint iL=0; iL<length; ++iL) {
            vec_offset.push_back(computeRelativeDisplacement(*lattices[0],*lattices[iL]));
        }
        plint numbM = (length-2)/2;
        for (plint iX0=domain.x0; iX0<=domain.x1; ++iX0) {
            for (plint iY0=domain.y0; iY0<=domain.y1; ++iY0) {
                plint iXm = iX0 + vec_offset[length-1].x; plint iYm = iY0 + vec_offset[length-1].y;
                plint mask = util::roundToInt( lattices[length-1]->get(iXm,iYm).computeDensity() );
                if ( mask > 2 ) {
                    T bmass = 0;
                    for (plint iB = 0; iB < numbM; ++iB) {
                        plint iXb = iX0 + vec_offset[iB].x; plint iYb = iY0 + vec_offset[iB].y;
                        bmass += lattices[iB]->get(iXb,iYb).computeDensity();
                    }
                    Array<T,5> g;
                    g[0]=(T) (bmass-1)/3; g[1]=g[2]=g[3]=g[4]=(T) (bmass-1)/6;
                    plint iXt = iX0 + vec_offset[length-2].x; plint iYt = iY0 + vec_offset[length-2].y;
                    lattices[length-2]->get(iXt,iYt).setPopulations(g);
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    virtual updateTOTbMassLattices2D<T,Descriptor>* clone() const {
        return new updateTOTbMassLattices2D<T,Descriptor>(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (plint iL=0; iL<length; ++iL) {
            if (iL == length-2) {
                modified[iL] = modif::staticVariables;
            }
            else {
                modified[iL] = modif::nothing;
            }
        }
    }
private:
    plint length;
};

// Update local dynamics of the mask lattice
template<typename T1, template<typename U> class Descriptor1, typename T2, template<typename U> class Descriptor2>
class updateMaskLatticeDynamics2D : public BoxProcessingFunctional2D_LL<T1,Descriptor1,T2,Descriptor2>
{
public:
    updateMaskLatticeDynamics2D(T bMf_, T Bmax_, plint pore_, plint solid_, plint bb_)
    : bMthrd(bMf_*Bmax_), pore(pore_), solid(solid_), bb(bb_)
    {}
    virtual void process(Box2D domain, BlockLattice2D<T1,Descriptor1>& lattice1, BlockLattice2D<T2,Descriptor2>& lattice2) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint mask = util::roundToInt( lattice1.get(iX,iY).computeDensity() );
                if (mask!=pore && mask!=solid && mask!=bb ) {
                    T bmass = lattice2.get(iX,iY).computeDensity();
                    plint omega = util::roundToInt( lattice2.get(iX,iY).getDynamics().getOmega() );
                    if (bmass >= bMthrd && omega == 0) {
                        lattice2.get(iX,iY).getDynamics().setOmega(1.);
                    }
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulk;
    }
    virtual updateMaskLatticeDynamics2D<T1,Descriptor1,T2,Descriptor2>* clone() const {
        return new updateMaskLatticeDynamics2D<T1,Descriptor1,T2,Descriptor2>(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::dynamicVariables;
    }
private:
    T bMthrd;
    plint pore, solid, bb;
};

// lattice divisions for MPI
template<typename T1, template<typename U> class Descriptor, typename T2>
class latticeXY2D : public BoxProcessingFunctional2D_LS<T1,Descriptor,T2>
{
public:
    latticeXY2D(bool biId_): biId(biId_)
    {}
    virtual void process(Box2D domain, BlockLattice2D<T1,Descriptor>& lattice, ScalarField2D<T2> &field) {
        Dot2D offset = computeRelativeDisplacement(lattice, field);
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            plint iX1 = iX + offset.x;
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint iY1 = iY + offset.y;
                if (biId == 0) {
                    field.get(iX1,iY1) = iX;
                }
                else {
                    field.get(iX1,iY1) = iY;
                }
            }
        }
    }
    virtual latticeXY2D<T1,Descriptor,T2>* clone() const {
        return new latticeXY2D<T1,Descriptor,T2>(*this);
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }
private:
    bool biId;
};

// absolute lattice indices
template<typename T1, template<typename U> class Descriptor, typename T2>
class latticeAbsoluteXY2D : public BoxProcessingFunctional2D_LS<T1,Descriptor,T2>
{
public:
    latticeAbsoluteXY2D(bool biId_): biId(biId_)
    {}
    virtual void process(Box2D domain, BlockLattice2D<T1,Descriptor>& lattice, ScalarField2D<T2> &field) {
        Dot2D offset = computeRelativeDisplacement(lattice, field);
        Dot2D absoluteOffset = lattice.getLocation();
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            plint absX = iX + absoluteOffset.x;
            plint iX1 = iX + offset.x;
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint absY = iY + absoluteOffset.y;
                plint iY1 = iY + offset.y;
                if (biId == 0) {
                    field.get(iX1,iY1) = absX;
                }
                else {
                    field.get(iX1,iY1) = absY;
                }
            }
        }
    }
    virtual latticeAbsoluteXY2D<T1,Descriptor,T2>* clone() const {
        return new latticeAbsoluteXY2D<T1,Descriptor,T2>(*this);
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }
private:
    bool biId;
};

// initialize scalar biomass lattice
template<typename T1, template<typename U> class Descriptor, typename T2>
class initializeMaskedScalarLattice : public BoxProcessingFunctional2D_LS<T1,Descriptor,T2>
{
public:
    initializeMaskedScalarLattice(T rho0_, plint id_): rho0(rho0_), id(id_)
    {}
    virtual void process(Box2D domain, BlockLattice2D<T1,Descriptor>& lattice, ScalarField2D<T2> &field) {
        Dot2D offset = computeRelativeDisplacement(lattice, field);
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            plint iX1 = iX + offset.x;
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint iY1 = iY + offset.y;
                T2 mask = field.get(iX1,iY1);
                if (mask == id) {
                    Array<T,5> g;
                    g[0]=(T) (rho0-1)/3; g[1]=g[2]=g[3]=g[4]=(T) (rho0-1)/6;
                    lattice.get(iX,iY).setPopulations(g);
                }
            }
        }
    }
    virtual initializeMaskedScalarLattice<T1,Descriptor,T2>* clone() const {
        return new initializeMaskedScalarLattice<T1,Descriptor,T2>(*this);
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
        modified[1] = modif::nothing;
    }
private:
    T rho0;
    plint id;
};


/* ===============================================================================================================
   ========================================== REDUCTIVE DATA PROCESSORS ==========================================
   =============================================================================================================== */

// count the number of a certain mask
template<typename T1, template<typename U1> class Descriptor1>
class BoxScalarSelectedSumFunctional2D : public ReductiveBoxProcessingFunctional2D_L<T1,Descriptor1>
{
public:
    BoxScalarSelectedSumFunctional2D(plint mask_) : countId(this->getStatistics().subscribeSum()), mask(mask_)
    {}
    virtual void process(Box2D domain, BlockLattice2D<T1,Descriptor1>& lattice) {
        BlockStatistics& statistics = this->getStatistics();
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint tmpMask = util::roundToInt(lattice.get(iX,iY).computeDensity());
                if ( tmpMask == mask ) {
                    statistics.gatherSum(countId, (int) 1);
                }
            }
        }
    }
    virtual BoxScalarSelectedSumFunctional2D<T1,Descriptor1>* clone() const {
        return new BoxScalarSelectedSumFunctional2D<T1,Descriptor1>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
    }
    T getCount() const {
        // The sum is internally computed on floating-point values.
        // If T is integer, the value must be rounded at the end.
        double doubleSum = this->getStatistics().getSum(countId);
        if (std::numeric_limits<T>::is_integer) {
            return (T) util::roundToInt(doubleSum);
        }
        return (T) doubleSum;
    }
private:
    plint countId;
    plint mask;
};

template<typename T1, template<typename U1> class Descriptor1>
T countLatticeMaskNumbers(Box2D domain, MultiBlockLattice2D<T1,Descriptor1>& lattice, plint mask) {
    BoxScalarSelectedSumFunctional2D<T1,Descriptor1> functional = BoxScalarSelectedSumFunctional2D<T1,Descriptor1> ( mask );
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getCount();
}

// sum up lattice omegas
template<typename T1, template<typename U1> class Descriptor1>
class SumLatticeCellOmegas2D : public ReductiveBoxProcessingFunctional2D_L<T1,Descriptor1>
{
public:
    SumLatticeCellOmegas2D() : omegasum(this->getStatistics().subscribeSum())
    {}
    virtual void process(Box2D domain, BlockLattice2D<T1,Descriptor1>& lattice) {
        BlockStatistics& statistics = this->getStatistics();
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint omega = util::roundToInt( lattice.get(iX,iY).getDynamics().getOmega() );
                statistics.gatherSum(omegasum, omega);
            }
        }
    }
    virtual SumLatticeCellOmegas2D<T1,Descriptor1>* clone() const {
        return new SumLatticeCellOmegas2D<T1,Descriptor1>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
    }
    T getSum() const {
        // The sum is internally computed on floating-point values.
        // If T is integer, the value must be rounded at the end.
        double doubleSum = this->getStatistics().getSum(omegasum);
        if (std::numeric_limits<T>::is_integer) {
            return (T) util::roundToInt(doubleSum);
        }
        return (T) doubleSum;
    }
private:
    plint omegasum;
};

template<typename T1, template<typename U1> class Descriptor1>
T sumUpLatticeOmegas(Box2D domain, MultiBlockLattice2D<T1,Descriptor1>& lattice) {
    SumLatticeCellOmegas2D<T1,Descriptor1> functional = SumLatticeCellOmegas2D<T1,Descriptor1> ();
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getSum();
}

// calculate the maximum density with a mask lattice
template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
class MaskedBoxLatticeMaxFunctional2D : public ReductiveBoxProcessingFunctional2D_LL<T1,Descriptor1,T2,Descriptor2>
{
public:
    MaskedBoxLatticeMaxFunctional2D(plint pore_, plint solid_, plint bb_) : maxLatticeId(this->getStatistics().subscribeMax()), pore(pore_), solid(solid_), bb(bb_)
    {}
    virtual void process(Box2D domain, BlockLattice2D<T1,Descriptor1> &lattice0, BlockLattice2D<T2,Descriptor2> &lattice1) {
        BlockStatistics& statistics = this->getStatistics();
        Dot2D offset_01 = computeRelativeDisplacement(lattice0,lattice1);
        for (plint iX0=domain.x0; iX0<=domain.x1; ++iX0) {
            for (plint iY0=domain.y0; iY0<=domain.y1; ++iY0) {
                plint mask = util::roundToInt( lattice1.get(iX0+offset_01.x,iY0+offset_01.y).computeDensity() );
                // if (mask > 2) {
                if (mask!=pore && mask!=solid && mask!=bb) {
                    T max = lattice0.get(iX0,iY0).computeDensity();
                    statistics.gatherMax(maxLatticeId, max);
                }
            }
        }
    }
    virtual MaskedBoxLatticeMaxFunctional2D<T1,Descriptor1,T2,Descriptor2>* clone() const {
        return new MaskedBoxLatticeMaxFunctional2D<T1,Descriptor1,T2,Descriptor2>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::nothing;
    }
    T getMaxLattice() const {
        // The sum is internally computed on floating-point values. If T is
        //   integer, the value must be rounded at the end.
        double doubleMax = this->getStatistics().getMax(maxLatticeId);
        if (std::numeric_limits<T>::is_integer) {
            return (T) util::roundToInt(doubleMax);
        }
        return (T) doubleMax;
    }
private:
    plint maxLatticeId, pore, solid, bb;
};

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
T computeMaskedLatticeMax(Box2D domain, MultiBlockLattice2D<T1,Descriptor1>& lattice0, MultiBlockLattice2D<T2,Descriptor2>& lattice1, plint pore, plint solid, plint bb) {
    MaskedBoxLatticeMaxFunctional2D<T1,Descriptor1,T2,Descriptor2> functional = MaskedBoxLatticeMaxFunctional2D<T1,Descriptor1,T2,Descriptor2>(pore, solid, bb);
    applyProcessingFunctional(functional, domain, lattice0, lattice1);
    return functional.getMaxLattice();
}

// calculate the minimum density with a mask lattice
template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
class MaskedBoxLatticeMinFunctional2D : public ReductiveBoxProcessingFunctional2D_LL<T1,Descriptor1,T2,Descriptor2>
{
public:
    MaskedBoxLatticeMinFunctional2D(plint mask1_, plint mask2_) : maxLatticeId(this->getStatistics().subscribeMax()), mask1(mask1_), mask2(mask2_)
    {}
    virtual void process(Box2D domain, BlockLattice2D<T1,Descriptor1> &lattice0, BlockLattice2D<T2,Descriptor2> &lattice1) {
        BlockStatistics& statistics = this->getStatistics();
        Dot2D offset_01 = computeRelativeDisplacement(lattice0,lattice1);
        for (plint iX0=domain.x0; iX0<=domain.x1; ++iX0) {
            for (plint iY0=domain.y0; iY0<=domain.y1; ++iY0) {
                plint tmpMask = lattice1.get(iX0+offset_01.x,iY0+offset_01.y).computeDensity();
                if (tmpMask == mask1 || tmpMask == mask2) {
                    T min = -lattice0.get(iX0,iY0).computeDensity();
                    statistics.gatherMax(maxLatticeId, min);
                }
            }
        }
    }
    virtual MaskedBoxLatticeMinFunctional2D<T1,Descriptor1,T2,Descriptor2>* clone() const {
        return new MaskedBoxLatticeMinFunctional2D<T1,Descriptor1,T2,Descriptor2>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::nothing;
    }
    T getMinLattice() const {
        // The sum is internally computed on floating-point values. If T is
        //   integer, the value must be rounded at the end.
        double doubleMin = - this->getStatistics().getMax(maxLatticeId);
        if (std::numeric_limits<T>::is_integer) {
            return (T) util::roundToInt(doubleMin);
        }
        return (T) doubleMin;
    }
private:
    plint maxLatticeId;
    plint mask1;
    plint mask2;
};

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
T computeMaskedLatticeMin(Box2D domain, MultiBlockLattice2D<T1,Descriptor1>& lattice0, MultiBlockLattice2D<T2,Descriptor2>& lattice1, plint mask1, plint mask2) {
    MaskedBoxLatticeMinFunctional2D<T1,Descriptor1,T2,Descriptor2> functional = MaskedBoxLatticeMinFunctional2D<T1,Descriptor1,T2,Descriptor2>(mask1, mask2);
    applyProcessingFunctional(functional, domain, lattice0, lattice1);
    return functional.getMinLattice();
}
