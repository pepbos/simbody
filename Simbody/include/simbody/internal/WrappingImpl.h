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
    public:

    using FrenetFrame = ContactGeometry::FrenetFrame;
    using Variation   = ContactGeometry::GeodesicVariation;
    using Correction  = ContactGeometry::GeodesicCorrection;

    //==============================================================================
    //                      ???
    //==============================================================================
    // Represents the local surface wrapping problem.
    // Caches last computed geodesic as a warmstart.
    // Not exposed outside of simbody.
    // Not shared amongst different paths/spans or obstacles/CurveSegments.
    class LocalGeodesic
    {
        using FrenetFrame = ContactGeometry::FrenetFrame;
        using Variation   = ContactGeometry::GeodesicVariation;
        using Correction  = ContactGeometry::GeodesicCorrection;

    public:
        LocalGeodesic()                                    = default;
        ~LocalGeodesic()                                   = default;
        LocalGeodesic(LocalGeodesic&&) noexcept            = default;
        LocalGeodesic& operator=(LocalGeodesic&&) noexcept = default;
        LocalGeodesic(const LocalGeodesic&)                = delete;
        LocalGeodesic& operator=(const LocalGeodesic&)     = delete;

        LocalGeodesic(
            CableSubsystem subsystem,
            ContactGeometry geometry,
            Vec3 initPointGuess) :
            m_Subsystem(subsystem),
            m_Geometry(geometry), m_InitPointGuess(initPointGuess)
        {}

        // Allocate state variables and cache entries.
        void realizeTopology(State& state);
        void realizePosition(const State& state) const;

        // Some info that can be retrieved from cache.
        struct LocalGeodesicInfo
        {
            FrenetFrame KP{};
            FrenetFrame KQ{};

            Real length = NaN;

            Variation dKP{};
            Variation dKQ{};

            Status status = Status::Ok;
        };

        // Helper struct: Required data for shooting a new geodesic.
        struct GeodesicInitialConditions
        {
        private:
            GeodesicInitialConditions() = default;

        public:
            static GeodesicInitialConditions CreateCorrected(
                const FrenetFrame& KP,
                const Variation& dKP,
                Real l,
                const Correction& c);
            static GeodesicInitialConditions CreateFromGroundInSurfaceFrame(
                const Transform& X_GS,
                Vec3 x_G,
                Vec3 t_G,
                Real l);
            static GeodesicInitialConditions CreateZeroLengthGuess(
                Vec3 prev_QS,
                Vec3 xGuess_S);
            static GeodesicInitialConditions CreateAtTouchdown(
                Vec3 prev_QS,
                Vec3 next_PS,
                Vec3 trackingPointOnLine);

            Vec3 x{NaN, NaN, NaN};
            Vec3 t{NaN, NaN, NaN};
            Real l = NaN;
        };

        bool analyticFormAvailable() const
        {
            return m_Geometry.analyticFormAvailable();
        }

        const LocalGeodesicInfo& calcInitialGeodesic(
            State& s,
            const GeodesicInitialConditions& g0) const;

        // This will reevaluate the cached geodesic and status.
        // This will be called by the curve segment before updating the
        // position level cache variable of the CurveSegment.
        // TODO Unfortunate naming convention.
        const LocalGeodesicInfo& calcLocalGeodesicInfo(
            const State& s,
            Vec3 prev_QS,
            Vec3 next_PS) const; // TODO weird name

        // Apply the correction to the initial condition of the geodesic, and
        // shoot a new geodesic, updating the cache variable.
        void applyGeodesicCorrection(const State& s, const Correction& c) const;

        // Compute the path points of the current geodesic, and write them to
        // the buffer. These points are in local surface coordinates. Returns
        // the number of points written.
        void calcPathPoints(const State& state, std::vector<Vec3>& points)
            const;

        // Set the user defined point that controls the initial wrapping path.
        void setInitialPointGuess(Vec3 initPointGuess)
        {
            m_InitPointGuess = initPointGuess;
        }

        // Get the user defined point that controls the initial wrapping path.
        Vec3 getInitialPointGuess() const
        {
            return m_InitPointGuess;
        }

        Status getStatus(const State& s) const
        {
            return getCacheEntry(s).status;
        }

    private:
        // The cache entry: Curve in local surface coordinated.
        // This is an auto update discrete cache variable, which makes it
        // persist over integration steps.
        // TODO or should it be "normal" discrete?
        struct CacheEntry : LocalGeodesicInfo
        {
            Vec3 trackingPointOnLine{NaN, NaN, NaN};
            std::vector<FrenetFrame> frames{}; // Empty for analytic geoemetry
                                               // with no allocation overhead.
            double sHint = NaN;
        };

        const CacheEntry& getCacheEntry(const State& state) const
        {
            return Value<CacheEntry>::downcast(
                m_Subsystem.getDiscreteVarUpdateValue(state, m_CacheIx));
        }

        CacheEntry& updCacheEntry(const State& state) const
        {
            return Value<CacheEntry>::updDowncast(
                m_Subsystem.updDiscreteVarUpdateValue(state, m_CacheIx));
        }

        const CacheEntry& getPrevCacheEntry(const State& state) const
        {
            return Value<CacheEntry>::downcast(
                m_Subsystem.getDiscreteVariable(state, m_CacheIx));
        }

        CacheEntry& updPrevCacheEntry(State& state) const
        {
            return Value<CacheEntry>::updDowncast(
                m_Subsystem.updDiscreteVariable(state, m_CacheIx));
        }

        void calcCacheEntry(
            const Vec3& prev_QS,
            const Vec3& next_PS,
            CacheEntry& cache) const;

        void assertSurfaceBounds(const Vec3& prev_QS, const Vec3& next_PS)
            const;

        void calcTouchdownIfNeeded(
            const Vec3& prev_QS,
            const Vec3& next_PS,
            CacheEntry& cache) const;

        void calcLiftoffIfNeeded(
            const Vec3& prev_QS,
            const Vec3& next_PS,
            CacheEntry& cache) const;

        void shootNewGeodesic(
            const GeodesicInitialConditions& g0,
            CacheEntry& cache) const;

        //------------------------------------------------------------------------------
        CableSubsystem m_Subsystem;

        ContactGeometry m_Geometry;

        Vec3 m_InitPointGuess;

        size_t m_ProjectionMaxIter        = 10;
        Real m_ProjectionRequiredAccuracy = 1e-10;
        Real m_IntegratorAccuracy         = 1e-6;

        Real m_TouchdownAccuracy = 1e-3;
        size_t m_TouchdownIter   = 10;

        // TODO this must be a function argument such that the caller can
        // decide,
        size_t m_NumberOfAnalyticPoints = 10;

        DiscreteVariableIndex m_CacheIx;
    };

private:
    Impl() = default;

public:
    Impl(const Impl& source)            = default;
    Impl& operator=(const Impl& source) = default;
    ~Impl()                             = default;

    Impl(
            CableSpan path,
            CurveSegmentIndex ix,
            const MobilizedBody& mobod,
            const Transform& X_BS,
            ContactGeometry geometry,
            Vec3 initPointGuess);

    // Position level cache: Curve in ground frame.
    struct PosInfo
    {
        Transform X_GS{};

        FrenetFrame KP{};
        FrenetFrame KQ{};

        Variation dKP{};
        Variation dKQ{};

        Real length = NaN;

        // TODO add force & moment here?
        /* Vec3 unitForce {}; */
        /* Vec3 unitMoment {}; */
    };

    // Allocate state variables and cache entries.
    void realizeTopology(State& state);
    /* void realizeInstance(const State& state) const; */
    void realizePosition(const State& state) const;

    /* void realizeVelocity(const State& state) const; */
    /* void realizeAcceleration(const State& state) const; */
    /* void invalidateTopology(); */

    void realizeCablePosition(const State& state) const;

    void invalidatePositionLevelCache(const State& state) const
    {
        m_Subsystem.markCacheValueNotRealized(state, m_PosInfoIx);
    }

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

    const PosInfo& getPosInfo(const State& s) const
    {
        realizePosition(s);
        return Value<PosInfo>::downcast(
            m_Subsystem.getCacheEntry(s, m_PosInfoIx));
    }

    bool isActive(const State& s) const
    {
        return m_Surface.getStatus(s) == Status::Ok;
    }

    Status getStatus(const State& s) const
    {
        return m_Surface.getStatus(s);
    }

    void calcInitZeroLengthGeodesic(State& state, Vec3 prev_QS) const;

    void applyGeodesicCorrection(
        const State& state,
        const ContactGeometry::GeodesicCorrection& c) const;

    void calcPathPoints(const State& state, std::vector<Vec3>& points) const;

    /* const MobilizedBody& getMobilizedBody() const {return m_Mobod;} */
    void calcContactPointVelocitiesInGround(
        const State& s,
        Vec3& v_GP,
        Vec3& v_GQ) const
    {
        // TODO use builtin?
        /* v_GP = m_Mobod.findStationVelocityInGround(state, P_B); */

        // Get body kinematics in ground frame.
        const Vec3& x_BG = m_Mobod.getBodyOriginLocation(s);
        const Vec3& w_BG = m_Mobod.getBodyAngularVelocity(s);
        const Vec3& v_BG = m_Mobod.getBodyOriginVelocity(s);

        const PosInfo& pos = getPosInfo(s);

        // Relative contact point positions to body origin, expressed in ground
        // frame.
        Vec3 r_P = pos.KP.p() - x_BG;
        Vec3 r_Q = pos.KQ.p() - x_BG;

        // Compute contact point velocities in ground frame.
        v_GP = v_BG + w_BG % r_P;
        v_GQ = v_BG + w_BG % r_Q;
    }

    Transform calcSurfaceFrameInGround(const State& s) const
    {
        return m_Mobod.getBodyTransform(s).compose(m_Offset);
    }

    SpatialVec calcAppliedWrenchInGround(const State& s, Real tension) const;

    // TODO Force? or wrench?
    void applyBodyForce(
        const State& state,
        Real tension,
        Vector_<SpatialVec>& bodyForcesInG) const;

    // TODO allow for user to shoot his own geodesic.
    /* void calcGeodesic(State& state, Vec3 x, Vec3 t, Real l) const; */

private:
    PosInfo& updPosInfo(const State& state) const
    {
        return Value<PosInfo>::updDowncast(
            m_Subsystem.updCacheEntry(state, m_PosInfoIx));
    }

    void calcPosInfo(const State& state, PosInfo& posInfo) const;

    // TODO Required for accessing the cache variable?
    CableSubsystem m_Subsystem; // The subsystem this segment belongs to.
    CableSpan m_Path;           // The path this segment belongs to.
    CurveSegmentIndex m_Index;  // The index in its path.

    MobilizedBody m_Mobod;
    Transform m_Offset;

    LocalGeodesic m_Surface; // TODO rename

    // TOPOLOGY CACHE
    CacheEntryIndex m_PosInfoIx;
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
        CableSubsystem subsystem,
        MobilizedBody originBody,
        Vec3 originPoint,
        MobilizedBody terminationBody,
        Vec3 terminationPoint) :
        m_Subsystem(subsystem),
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

    CurveSegmentIndex adoptObstacle(CurveSegment& segment)
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
        m_Subsystem.invalidateSubsystemTopologyCache();
    }

    const PosInfo& getPosInfo(const State& state) const;
    const VelInfo& getVelInfo(const State& state) const;

    void calcInitZeroLengthGeodesic(State& s) const;

    void applyBodyForces(
        const State& state,
        Real tension,
        Vector_<SpatialVec>& bodyForcesInG) const;

private:
    PosInfo& updPosInfo(const State& s) const;
    VelInfo& updVelInfo(const State& state) const;

    void calcPosInfo(const State& s, PosInfo& posInfo) const;
    void calcVelInfo(const State& s, VelInfo& velInfo) const;

    size_t countActive(const State& s) const;

    Vec3 findPrevPoint(const State& state, CurveSegmentIndex ix) const;

    Vec3 findNextPoint(const State& state, CurveSegmentIndex ix) const;

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
        return m_Subsystem;
    }

    // Reference back to the subsystem.
    CableSubsystem m_Subsystem;

    MobilizedBody m_OriginBody;
    Vec3 m_OriginPoint;

    MobilizedBody m_TerminationBody;
    Vec3 m_TerminationPoint;

    Array_<CurveSegment, CurveSegmentIndex> m_CurveSegments{};

    Real m_PathErrorBound = 0.1;
    Real m_ObsErrorBound  = 0.1;
    size_t m_PathMaxIter  = 10;
    size_t m_ObsMaxIter   = 10;

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

        for (CableSpanIndex ix(0); ix < cables.size(); ++ix) {
            CableSpan& path = wThis->updCablePath(ix);
            path.updImpl().realizeTopology(state);
        }

        return 0;
    }

    int realizeSubsystemPositionImpl(const State& state) const override
    {
        for (CableSpanIndex ix(0); ix < cables.size(); ++ix) {
            const CableSpan& path = getCablePath(ix);
            path.getImpl().realizePosition(state);
        }
        return 0;
    }

    int realizeSubsystemVelocityImpl(const State& state) const override
    {
        for (CableSpanIndex ix(0); ix < cables.size(); ++ix) {
            const CableSpan& path = getCablePath(ix);
            path.getImpl().realizeVelocity(state);
        }
        return 0;
    }

    SimTK_DOWNCAST(Impl, Subsystem::Guts);

private:
    // TOPOLOGY STATE
    Array_<CableSpan, CableSpanIndex> cables;
};

} // namespace SimTK

#endif
