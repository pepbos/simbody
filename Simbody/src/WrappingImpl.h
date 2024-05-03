#ifndef SimTK_SIMBODY_WRAPPING_PATH_IMPL_H_
#define SimTK_SIMBODY_WRAPPING_PATH_IMPL_H_

#include "SimTKcommon/internal/State.h"
#include "SimTKmath.h"
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/Wrapping.h"
#include "simbody/internal/common.h"
#include "simmath/internal/ContactGeometry.h"
#include <functional>
#include <stdexcept>

namespace SimTK
{

//==============================================================================
//                                ??? IMPL
//==============================================================================
class CurveSegment::Impl
{

//------------------------------------------------------------------------------
//                        Some public types and aliases
//------------------------------------------------------------------------------
public:
    using FrenetFrame = ContactGeometry::FrenetFrame;
    using Variation   = ContactGeometry::GeodesicVariation;
    using Correction  = ContactGeometry::GeodesicCorrection;

    // Instance level cache entry.
    struct LocalGeodesicSample
    {
        LocalGeodesicSample(Real l, FrenetFrame K) : length(l), frame(K)
        {}

        Real length;
        FrenetFrame frame;
    };

    // Instance level cache entry.
    struct InstanceEntry
    {
        bool isActive() const
        {
            return status == Status::Ok;
        }

        FrenetFrame K_P{};
        FrenetFrame K_Q{};

        Real length = NaN;

        Variation dK_P{};
        Variation dK_Q{};

        std::vector<LocalGeodesicSample> samples;
        double sHint = NaN;

        Vec3 trackingPointOnLine{NaN, NaN, NaN};
        Status status = Status::Lifted;
    };

    // Position level cache: Curve in ground frame.
    struct PosInfo
    {
        Transform X_GS{};

        FrenetFrame KP{};
        FrenetFrame KQ{};

        Variation dKP{};
        Variation dKQ{};

        SpatialVec unitForce;
    };
//------------------------------------------------------------------------------

public:
    Impl()                           = delete;
    Impl(const Impl&)                = delete;
    Impl& operator=(const Impl&)     = delete;

    ~Impl()                          = default;
    Impl(Impl&&) noexcept            = default;
    Impl& operator=(Impl&&) noexcept = default;

    // TODO you would expect the constructor to take the index as well here?
    Impl(
        CableSpan path,
        const MobilizedBody& mobod,
        const Transform& X_BS,
        ContactGeometry geometry,
        Vec3 initPointGuess);

//------------------------------------------------------------------------------
//                          Stage realizations
//------------------------------------------------------------------------------
    // Allocate state variables and cache entries.
    void realizeTopology(State& s)
    {
        // Allocate an auto-update discrete variable for the last computed
        // geodesic.
        Value<InstanceEntry>* cache = new Value<InstanceEntry>();

        cache->upd().length              = 0.;
        cache->upd().status              = Status::Lifted;
        cache->upd().trackingPointOnLine = getContactPointHint();

        m_InstanceIx = updSubsystem().allocateAutoUpdateDiscreteVariable(
            s,
            Stage::Report,
            cache,
            Stage::Position);

        PosInfo posInfo{};
        m_PosIx = updSubsystem().allocateCacheEntry(
            s,
            Stage::Position,
            Stage::Infinity,
            new Value<PosInfo>(posInfo));
    }

    void realizePosition(const State& s, Vec3 prevPointG, Vec3 nextPointG) const
    {
        if (getSubsystem().isCacheValueRealized(s, m_PosIx)) {
            throw std::runtime_error(
                "expected not realized when calling realizePosition");
        }

        if (getInstanceEntry(s).status == Status::Disabled) {
            return;
        }

        // Compute tramsform from local surface frame to ground.
        const Transform& X_GS = calcSurfaceFrameInGround(s);

        {
            // Transform the prev and next path points to the surface frame.
            const Vec3 prevPoint_S = X_GS.shiftBaseStationToFrame(prevPointG);
            const Vec3 nextPoint_S = X_GS.shiftBaseStationToFrame(nextPointG);

            // Detect liftoff, touchdown and potential invalid configurations.
            // TODO this doesnt follow the regular invalidation scheme...
            // Grab the last geodesic that was computed.
            assertSurfaceBounds(prevPoint_S, nextPoint_S);

            calcTouchdownIfNeeded(prevPoint_S, nextPoint_S, updInstanceEntry(s));
            calcLiftoffIfNeeded(prevPoint_S, nextPoint_S, updInstanceEntry(s));

            getSubsystem().markDiscreteVarUpdateValueRealized(s, m_InstanceIx);
        }

        // At this point we have a valid geodesic in surface frame.
        const InstanceEntry& ie = getInstanceEntry(s);

        // Start updating the position level cache.
        PosInfo& pos            = updPosInfo(s);

        // Transform geodesic in local surface coordinates to ground.
        {
            // Store the local geodesic in ground frame.
            pos.X_GS = X_GS;

            // Store the the local geodesic in ground frame.
            pos.KP = X_GS.compose(ie.K_P);
            pos.KQ = X_GS.compose(ie.K_Q);

            pos.dKP[0] = X_GS.R() * ie.dK_P[0];
            pos.dKP[1] = X_GS.R() * ie.dK_P[1];

            pos.dKQ[0] = X_GS.R() * ie.dK_Q[0];
            pos.dKQ[1] = X_GS.R() * ie.dK_Q[1];
        }

        getSubsystem().markCacheValueRealized(s, m_PosIx);
    }

    void realizeCablePosition(const State& s) const;

    void invalidatePosEntry(const State& state) const
    {
        getSubsystem().markCacheValueNotRealized(state, m_PosIx);
    }

//------------------------------------------------------------------------------
//                  Parameter & model components access
//------------------------------------------------------------------------------
    const CableSpan& getCable() const
    {
        return m_Path;
    }

    CurveSegmentIndex getIndex() const
    {
        return m_Index;
    }

    void setIndex(CurveSegmentIndex ix)
    {
        m_Index = ix;
    }

    // Set the user defined point that controls the initial wrapping path.
    // Point is in surface coordinates.
    void setContactPointHint(Vec3 contactPointHint_S)
    {
        m_ContactPointHint_S = contactPointHint_S;
    }

    // Get the user defined point that controls the initial wrapping path.
    Vec3 getContactPointHint() const
    {
        return m_ContactPointHint_S;
    }

    const MobilizedBody& getMobilizedBody() const {return m_Mobod;}

    CableSubsystem& updSubsystem() { return *m_Subsystem; }
    const CableSubsystem& getSubsystem() const { return *m_Subsystem; }

    const DecorativeGeometry& getDecoration() const { return m_Decoration; }

//------------------------------------------------------------------------------
//                         Cache entry access
//------------------------------------------------------------------------------
    const InstanceEntry& getInstanceEntry(const State& s) const
    {
        const CableSubsystem& subsystem = getSubsystem();
        if (!subsystem.isDiscreteVarUpdateValueRealized(s, m_InstanceIx)) {
            updInstanceEntry(s) = getPrevInstanceEntry(s);
            subsystem.markDiscreteVarUpdateValueRealized(s, m_InstanceIx);
        }
        return Value<InstanceEntry>::downcast(
            subsystem.getDiscreteVarUpdateValue(s, m_InstanceIx));
    }

    const PosInfo& getPosInfo(const State& s) const
    {
        return Value<PosInfo>::downcast(
            getSubsystem().getCacheEntry(s, m_PosIx));
    }

    Transform calcSurfaceFrameInGround(const State& s) const
    {
        return m_Mobod.getBodyTransform(s).compose(m_X_BS);
    }

    int calcPathPoints(const State& s, std::vector<Vec3>& points, int nSamples=0) const;

    void calcUnitForce(const State& s, SpatialVec& unitForce_G) const
    {
        const PosInfo& posInfo = getPosInfo(s);
        const Vec3& x_BG = m_Mobod.getBodyOriginLocation(s);

        // Contact point moment arms in ground.
        const Vec3 r_P = posInfo.KP.p() - x_BG;
        const Vec3 r_Q = posInfo.KQ.p() - x_BG;

        // Tangent directions at contact points in ground.
        const UnitVec3& t_P = posInfo.KP.R().getAxisUnitVec(TangentAxis);
        const UnitVec3& t_Q = posInfo.KQ.R().getAxisUnitVec(TangentAxis);

        unitForce_G[0] = r_Q % t_Q - r_P % t_P;
        unitForce_G[1] = t_Q - t_P;
    }

    // Apply the correction to the initial condition of the geodesic, and
    // shoot a new geodesic, updating the cache variable.
    void applyGeodesicCorrection(const State& s, const Correction& c) const
    {
        // Get the previous geodesic.
        const InstanceEntry& cache = getInstanceEntry(s);
        const Variation& dK_P      = cache.dK_P;
        const FrenetFrame& K_P     = cache.K_P;

        // TODO This is part of the unit test block below...
        const FrenetFrame K0_P = cache.K_P;
        const FrenetFrame K0_Q = cache.K_Q;
        const Variation dK0_P  = cache.dK_P;
        const Variation dK0_Q  = cache.dK_Q;
        // TODO end of part belongning to unit test below..

        // Get corrected initial conditions.
        const Vec3 v     = dK_P[1] * c;
        const Vec3 point = K_P.p() + v;

        const Vec3 w       = dK_P[0] * c;
        const UnitVec3 t   = K_P.R().getAxisUnitVec(TangentAxis);
        const Vec3 tangent = t + cross(w, t);

        // Take the length correction, and add to the current length.
        const Real dl = c[3]; // Length increment is the last element.
        const Real length =
            std::max(cache.length + dl, 0.); // Clamp length to be nonnegative.

        // Shoot the new geodesic.
        shootNewGeodesic(
            point,
            tangent,
            length,
            cache.sHint,
            updInstanceEntry(s));

        getSubsystem().markDiscreteVarUpdateValueRealized(s, m_InstanceIx);

        // TODO Below code should be moved to a unit test...
        const bool DO_UNIT_TEST_HERE = true;
        if (DO_UNIT_TEST_HERE)
        {
            const Real delta = c.norm();
            const Real eps   = delta / 10.;
            auto AssertAxis  = [&](const Rotation& R0,
                    const Rotation& R1,
                    const Vec3& w,
                    CoordinateAxis axis) -> bool {
                const UnitVec3 a0        = R0.getAxisUnitVec(axis);
                const UnitVec3 a1        = R1.getAxisUnitVec(axis);
                const Vec3 expected_diff = cross(w, a0);
                const Vec3 got_diff      = a1 - a0;
                const bool isOk          = (expected_diff - got_diff).norm() < eps;
                if (!isOk) {
                    std::cout << "    a0 = " << R0.transpose() * a0 << "\n";
                    std::cout << "    a1 = " << R0.transpose() * a1 << "\n";
                    std::cout << "    expected diff = "
                        << R0.transpose() * expected_diff / delta << "\n";
                    std::cout << "    got      dt_Q = "
                        << R0.transpose() * got_diff / delta << "\n";
                    std::cout << "    err           = "
                        << (expected_diff - got_diff).norm() / delta << "\n";
                }
                return isOk;
            };

            auto AssertFrame = [&](const Transform& K0,
                    const Transform& K1,
                    const Variation& dK0) -> bool {
                const Vec3 dx_got      = K1.p() - K0.p();
                const Vec3 dx_expected = dK0[1] * c;
                bool isOk              = (dx_expected - dx_got).norm() < eps;
                if (!isOk) {
                    std::cout << "Apply variation c = " << c << "\n";
                    std::cout << "    x0 = " << K0.p() << "\n";
                    std::cout << "    x1 = " << K1.p() << "\n";
                    std::cout << "    expected dx_Q = "
                        << K0.R().transpose() * dx_expected / delta << "\n";
                    std::cout << "    got      dx_Q = "
                        << K0.R().transpose() * dx_got / delta << "\n";
                    std::cout << "    err           = "
                        << (dx_expected - dx_got).norm() / delta << "\n";
                    std::cout << "WARNING: Large deviation in final position\n";
                }
                const Vec3 w                       = dK0[0] * c;
                std::array<CoordinateAxis, 3> axes = {
                    TangentAxis,
                    NormalAxis,
                    BinormalAxis};
                for (size_t i = 0; i < 3; ++i) {
                    isOk &= AssertAxis(K0.R(), K1.R(), w, axes.at(i));
                    if (!isOk) {
                        std::cout << "FAILED axis " << i << "\n";
                    }
                }

                return isOk;
            };

            if (delta > 1e-10) {
                if (!AssertFrame(K0_P, getInstanceEntry(s).K_P, dK0_P)) {
                    throw std::runtime_error("Start frame variation check failed");
                }
                if (!AssertFrame(K0_Q, getInstanceEntry(s).K_Q, dK0_Q)) {
                    throw std::runtime_error("End frame variation check failed");
                }
            }
        }
        // End of unit test code block...

        invalidatePosEntry(s);
    }

    Vec3 calcInitialContactPoint(const State& s) const
    {
        const InstanceEntry& ic = getInstanceEntry(s);
        if (!ic.isActive()) {
            throw std::runtime_error(
                "Invalid contact point: Curve is not active");
        }
        const Vec3& x_PS = ic.K_P.p();
        return calcSurfaceFrameInGround(s).shiftFrameStationToBase(x_PS);
    }

    Vec3 calcFinalContactPoint(const State& s) const
    {
        const InstanceEntry& ic = getInstanceEntry(s);
        if (!ic.isActive()) {
            throw std::runtime_error(
                "Invalid contact point: Curve is not active");
        }
        const Vec3& x_QS = ic.K_Q.p();
        return calcSurfaceFrameInGround(s).shiftFrameStationToBase(x_QS);
    }

    void calcDecorativeGeometryAndAppend(
        const State& s,
        Array_<DecorativeGeometry>& decorations) const
    {
        const InstanceEntry& cache = getInstanceEntry(s);
        if (!cache.isActive()) {return;}

        const Transform& X_GS      = calcSurfaceFrameInGround(s);
        Vec3 a = X_GS.shiftFrameStationToBase(cache.K_P.p());
        for (size_t i = 1; i < cache.samples.size(); ++i) {
            const Vec3 b =
                X_GS.shiftFrameStationToBase(cache.samples.at(i).frame.p());
            decorations.push_back(
                DecorativeLine(a, b).setColor(Purple).setLineThickness(3));
            a = b;
        }
    }

private:

    InstanceEntry& updInstanceEntry(const State& state) const
    {
        return Value<InstanceEntry>::updDowncast(
            getSubsystem().updDiscreteVarUpdateValue(state, m_InstanceIx));
    }

    const InstanceEntry& getPrevInstanceEntry(const State& state) const
    {
        return Value<InstanceEntry>::downcast(
            getSubsystem().getDiscreteVariable(state, m_InstanceIx));
    }

    InstanceEntry& updPrevInstanceEntry(State& state) const
    {
        return Value<InstanceEntry>::updDowncast(
            getSubsystem().updDiscreteVariable(state, m_InstanceIx));
    }

    PosInfo& updPosInfo(const State& state) const
    {
        return Value<PosInfo>::updDowncast(
            getSubsystem().updCacheEntry(state, m_PosIx));
    }

//------------------------------------------------------------------------------
//                          Local geodesic helpers
//------------------------------------------------------------------------------
    // Assert that previous and next point lie above the surface. Points are in
    // surface coordinates.
    void assertSurfaceBounds(const Vec3& prevPoint_S, const Vec3& nextPoint_S)
        const;

    // Attempt to compute the point of touchdown on the surface.
    void calcTouchdownIfNeeded(
        const Vec3& prevPoint_S,
        const Vec3& nextPoint_S,
        InstanceEntry& cache) const;

    void calcLiftoffIfNeeded(
        const Vec3& prevPoint_S,
        const Vec3& nextPoint_S,
        InstanceEntry& cache) const;

    // Compute a new geodesic in surface frame coordinates.
    void shootNewGeodesic(
        Vec3 point_S,
        Vec3 tangent_S,
        Real length,
        Real stepSizeHint,
        InstanceEntry& cache) const;
//------------------------------------------------------------------------------

    // TODO Required for accessing the cache variable?
    CableSubsystem* m_Subsystem; // The subsystem this segment belongs to.
    CableSpan m_Path;            // The path this segment belongs to.
    CurveSegmentIndex m_Index;   // The index in its path.

    // MobilizedBody that surface is attached to.
    MobilizedBody m_Mobod;
    // Surface to body transform.
    Transform m_X_BS;

    // Decoration TODO should this be here?
    ContactGeometry m_Geometry;
    DecorativeGeometry m_Decoration;

    // TOPOLOGY CACHE
    CacheEntryIndex m_PosIx;
    DiscreteVariableIndex m_InstanceIx;

    Vec3 m_ContactPointHint_S{NaN, NaN, NaN};

    size_t m_ProjectionMaxIter = 10;
    Real m_ProjectionAccuracy  = 1e-10;

    Real m_IntegratorAccuracy = 1e-8;

    Real m_TouchdownAccuracy = 1e-4;
    size_t m_TouchdownIter   = 10;
};

//==============================================================================
//                         PATH :: IMPL
//==============================================================================
// TODO
// - handout CurveSegment index
// - add other cache variables: velocity, acceleration, force
// - names: isValid, ContactStationOnBody, Entry (not info/cache)
class CableSpan::Impl
{
public:
    Impl(
        CableSubsystem& subsystem,
        MobilizedBody originBody,
        Vec3 originPoint,
        MobilizedBody terminationBody,
        Vec3 terminationPoint) :
        m_Subsystem(&subsystem),
        m_OriginBody(originBody), m_OriginPoint(originPoint),
        m_TerminationBody(terminationBody), m_TerminationPoint(terminationPoint)
    {}

    int getNumCurveSegments() const
    {
        return m_CurveSegments.size();
    }

    const CurveSegment& getCurveSegment(CurveSegmentIndex ix) const
    {
        return m_CurveSegments[ix];
    }

    CurveSegmentIndex adoptSegment(const CurveSegment& segment)
    {
        m_CurveSegments.push_back(segment);
        return CurveSegmentIndex(m_CurveSegments.size() - 1);
    }

    // Position level cache entry.
    struct PosInfo
    {
        Vec3 xO{NaN, NaN, NaN};
        Vec3 xI{NaN, NaN, NaN};

        Real l = NaN;

        size_t loopIter = 0;
    };

    // Velocity level cache entry.
    struct VelInfo
    {
        Real lengthDot = NaN;
    };

    // Allocate state variables and cache entries.
    void realizeTopology(State& state);
    void realizePosition(const State& state) const;
    void realizeVelocity(const State& state) const;
    void invalidateTopology()
    {
        getSubsystem().invalidateSubsystemTopologyCache();
    }
    void invalidatePositionLevelCache(const State& state) const;

    const PosInfo& getPosInfo(const State& state) const;
    const VelInfo& getVelInfo(const State& state) const;

    void applyBodyForces(
        const State& state,
        Real tension,
        Vector_<SpatialVec>& bodyForcesInG) const;

    int calcDecorativeGeometryAndAppend(
        const State& state,
        Stage stage,
        Array_<DecorativeGeometry>& decorations) const;

    void calcPathPoints(const State& state, std::vector<Vec3>& points) const
    {
        points.push_back(getPosInfo(state).xO);
        for (const CurveSegment& curve : m_CurveSegments) {
            curve.getImpl().calcPathPoints(state, points);
        }
        points.push_back(getPosInfo(state).xI);
    }

    void calcPathErrorJacobianUtility(
        const State& state,
        Vector& pathError,
        Matrix& pathErrorJacobian) const;

    void applyCorrection(const State& state, const Vector& correction) const;

    size_t countActive(const State& s) const;

private:
    PosInfo& updPosInfo(const State& s) const;
    VelInfo& updVelInfo(const State& state) const;

    void calcPosInfo(const State& s, PosInfo& posInfo) const;
    void calcVelInfo(const State& s, VelInfo& velInfo) const;

    Vec3 findPrevPoint(const State& state, const CurveSegment& curve) const;

    Vec3 findNextPoint(const State& state, const CurveSegment& curve) const;

    const CurveSegment* findPrevActiveCurveSegment(
        const State& s,
        CurveSegmentIndex ix) const;
    const CurveSegment* findNextActiveCurveSegment(
        const State& s,
        CurveSegmentIndex ix) const;

    template <size_t N>
    void calcPathErrorVector(
        const State& state,
        const std::vector<LineSegment>& lines,
        std::array<CoordinateAxis, N> axes,
        Vector& pathError) const;

    template <size_t N>
    void calcPathErrorJacobian(
        const State& state,
        const std::vector<LineSegment>& lines,
        std::array<CoordinateAxis, N> axes,
        Matrix& J) const;

    // Make static or not?
    void calcLineSegments(
        const State& s,
        Vec3 p_O,
        Vec3 p_I,
        std::vector<LineSegment>& lines) const;

    double calcPathLength(
        const State& state,
        const std::vector<LineSegment>& lines) const;

    const CableSubsystem& getSubsystem() const
    {
        return *m_Subsystem;
    }

    CableSubsystem& updSubsystem()
    {
        return *m_Subsystem;
    }

    // Reference back to the subsystem.
    CableSubsystem* m_Subsystem; // TODO just a pointer?

    MobilizedBody m_OriginBody;
    Vec3 m_OriginPoint;

    MobilizedBody m_TerminationBody;
    Vec3 m_TerminationPoint;

    Array_<CurveSegment, CurveSegmentIndex> m_CurveSegments{};

    Real m_PathErrorBound = 1e-4;
    size_t m_PathMaxIter  = 50;

    // TOPOLOGY CACHE (set during realizeTopology())
    CacheEntryIndex m_PosInfoIx;
    CacheEntryIndex m_VelInfoIx;

    friend CurveSegment::Impl;
};

//==============================================================================
//                         SUBSYSTEM :: IMPL
//==============================================================================
class CableSubsystem::Impl : public Subsystem::Guts
{
public:
    struct SolverData
    {
        SolverData(int nActive)
        {
            static constexpr int Q = 4;
            static constexpr int C = 4;
            const int n            = nActive;

            lineSegments.resize(n + 1);
            pathErrorJacobian = Matrix(C * n, Q * n, 0.);
            pathCorrection    = Vector(Q * n, 0.);
            pathError         = Vector(C * n, 0.);
            mat               = Matrix(Q * n, Q * n, NaN);
            vec               = Vector(Q * n, NaN);
        }

        std::vector<CableSpan::LineSegment> lineSegments;

        Matrix pathErrorJacobian;
        Vector pathCorrection;
        Vector pathError;
        Matrix mat;
        // TODO Cholesky decomposition...
        FactorLU matInv;
        Vector vec;
    };

    struct CacheEntry
    {
        CacheEntry() = default;
        SolverData& updOrInsert(int nActive)
        {
            if (nActive <= 0) {
                throw std::runtime_error(
                    "Cannot produce solver data of zero dimension");
            }

            for (int i = m_Data.size(); i < nActive; ++i) {
                m_Data.emplace_back(i + 1);
            }

            return m_Data.at(nActive - 1);
        }

        std::vector<SolverData> m_Data;
    };

    Impl()
    {}
    ~Impl()
    {}
    Impl* cloneImpl() const override
    {
        return new Impl(*this);
    }

    int getNumPaths() const
    {
        return cables.size();
    }

    const CableSpan& getCablePath(CableSpanIndex index) const
    {
        return cables[index];
    }

    CableSpan& updCablePath(CableSpanIndex index)
    {
        return cables[index];
    }

    // Add a cable path to the list, bumping the reference count.
    CableSpanIndex adoptCablePath(CableSpan& path)
    {
        cables.push_back(path);
        return CableSpanIndex(cables.size() - 1);
    }

    // Return the MultibodySystem which owns this WrappingPathSubsystem.
    const MultibodySystem& getMultibodySystem() const
    {
        return MultibodySystem::downcast(getSystem());
    }

    // Return the SimbodyMatterSubsystem from which this WrappingPathSubsystem
    // gets the bodies to track.
    const SimbodyMatterSubsystem& getMatterSubsystem() const
    {
        return getMultibodySystem().getMatterSubsystem();
    }

    // Allocate state variables.
    int realizeSubsystemTopologyImpl(State& state) const override
    {
        // Briefly allow writing into the Topology cache; after this the
        // Topology cache is const.
        Impl* wThis = const_cast<Impl*>(this);

        wThis->realizeTopology(state);
        for (CableSpanIndex ix(0); ix < cables.size(); ++ix) {
            CableSpan& path = wThis->updCablePath(ix);
            path.updImpl().realizeTopology(state);
        }

        return 0;
    }

    void realizeTopology(State& state)
    {
        CacheEntry cache{};
        m_CacheIx = allocateCacheEntry(
            state,
            Stage::Instance,
            Stage::Infinity,
            new Value<CacheEntry>(cache));
    }

    CacheEntry& updCachedScratchboard(const State& state) const
    {
        return Value<CacheEntry>::updDowncast(updCacheEntry(state, m_CacheIx));
    }

    int calcDecorativeGeometryAndAppendImpl(
        const State& state,
        Stage stage,
        Array_<DecorativeGeometry>& decorations) const override
    {
        if (stage != Stage::Position)
            return 0;

        for (const CableSpan& cable : cables) {
            int returnValue = cable.getImpl().calcDecorativeGeometryAndAppend(
                state,
                stage,
                decorations);
            if (returnValue != 0) {
                return returnValue;
            }
        }
        return 0;
    }

    SimTK_DOWNCAST(Impl, Subsystem::Guts);

private:
    // TOPOLOGY STATE
    Array_<CableSpan, CableSpanIndex> cables;

    CacheEntryIndex m_CacheIx;
};

} // namespace SimTK

#endif
