#ifndef SimTK_SIMBODY_WRAPPING_PATH_IMPL_H_
#define SimTK_SIMBODY_WRAPPING_PATH_IMPL_H_

#include "SimTKcommon/internal/State.h"
#include "SimTKmath.h"
#include "simbody/internal/Wrapping.h"
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/common.h"
#include "simmath/internal/ContactGeometry.h"

#include <functional>
#include <stdexcept>

namespace SimTK
{

SimTK_DEFINE_UNIQUE_INDEX_TYPE(PathSegmentIndex);
SimTK_DEFINE_UNIQUE_INDEX_TYPE(PathIndex);

//==============================================================================
//                                OBSTACLE IMPL
//==============================================================================
class WrapObstacle::Impl
{
    using FrenetFrame = ContactGeometry::FrenetFrame;
    using Variation = ContactGeometry::GeodesicVariation;
    using Correction = ContactGeometry::GeodesicCorrection;

    private:
    Impl()                                      = default;

    public:
    Impl(const Impl& source)            = default;
    Impl& operator=(const Impl& source) = default;
    ~Impl()                                     = default;

    /* Surface(const MobilizedBody& mobod, const Transform& X_BS, const ContactGeometry& geometry) */
    /*     : impl(std::make_shared<SurfaceImpl>(mobod, X_BS, geometry)) {} */
    /* Transform calcSurfaceToGroundTransform(const State& state) const {return impl->calcSurfaceToGroundTransform(state);} */

    Impl(
            WrappingPathSubsystem subsystem,
            const MobilizedBody& mobod,
            const Transform& X_BS,
            ContactGeometry geometry,
            Vec3 initPointGuess
        ) : 
        m_Subsystem(subsystem),
        m_Mobod(mobod),
        m_Offset(X_BS),
        m_Surface(subsystem, geometry, initPointGuess)
    {}

    enum class Status
    {
        Ok,
        Liftoff,
        Disabled,
    };

    // Ground frame solution.
    struct PosInfo
    {
        FrenetFrame KP {};
        FrenetFrame KQ {};

        Variation dKP {};
        Variation dKQ {};

        Real length = NaN;

        Status status = Status::Ok;
    };

    // Allocate state variables and cache entries.
    void realizeTopology(State& state);
    /* void realizeInstance(const State& state) const; */
    void realizePosition(const State& state) const;
    /* void realizeVelocity(const State& state) const; */
    /* void realizeAcceleration(const State& state) const; */
    /* void invalidateTopology(); */

    void invalidatePositionLevelCache(const State& state) const
    {
        m_Subsystem.markCacheValueNotRealized(state, m_PosInfoIx);
    }

    const PosInfo& getPosInfo(const State& state) const;

    bool isActive(const State& state) const;
    bool isDisabled(const State& state) const;

    PosInfo& updPosInfo(const State &state) const
    {
        return Value<PosInfo>::updDowncast(m_Subsystem.updCacheEntry(state, m_PosInfoIx));
    }

    size_t calcPathPoints(const State& state, std::vector<Vec3>& points) const;
    void applyGeodesicCorrection(const State& state, const ContactGeometry::GeodesicCorrection& c) const;
    void calcInitZeroLengthGeodesic(State& state, Vec3 prev_QS);

    private:
    void calcPosInfo(const State& state, PosInfo& posInfo) const;
    void calcGeodesic(State& state, Vec3 x, Vec3 t, Real l) const;

    public: // TODO public?
    //==============================================================================
    //                      SURFACE
    //==============================================================================
    // Represents the local surface wrapping problem.
    // Caches last computed geodesic as a warmstart.
    // Not exposed outside of simbody.
    // Not shared amongst different paths or obstacles.
    class Surface
    {
        using FrenetFrame = ContactGeometry::FrenetFrame;
        using Variation = ContactGeometry::GeodesicVariation;
        using Correction = ContactGeometry::GeodesicCorrection;

        public:
        Surface()                              = default;
        ~Surface() = default;
        Surface(Surface&&) noexcept            = default;
        Surface& operator=(Surface&&) noexcept = default;
        Surface(const Surface&)                = delete;
        Surface& operator=(const Surface&)     = delete;

        Surface(
                WrappingPathSubsystem subsystem,
                ContactGeometry geometry,
                Vec3 initPointGuess
               ) : 
            m_Subsystem(subsystem),
            m_Geometry(geometry),
            m_InitPointGuess(initPointGuess)
        {}

        // Allocate state variables and cache entries.
        void realizeTopology(State& state);
        void realizePosition(const State& state) const;

        struct LocalGeodesicInfo
        {
            FrenetFrame KP {};
            FrenetFrame KQ {};

            Real length = NaN;

            Variation dKP {};
            Variation dKQ {};

            Status status = Status::Ok;
        };

        struct GeodesicInitialConditions
        {
            private:
                GeodesicInitialConditions() = default;

            public:
            static GeodesicInitialConditions CreateCorrected(const FrenetFrame& KP, const Variation& dKP, Real l, const Correction& c);
            static GeodesicInitialConditions CreateInSurfaceFrame(const Transform& X_GS, Vec3 x_G, Vec3 t_G, Real l);
            static GeodesicInitialConditions CreateZeroLengthGuess(const Transform& X_GS, Vec3 prev_QS, Vec3 xGuess_S);
            static GeodesicInitialConditions CreateAtTouchdown(Vec3 prev_QS, Vec3 next_PS, Vec3 trackingPointOnLine);

            Vec3 x {NaN, NaN, NaN};
            Vec3 t {NaN, NaN, NaN};
            Real l = NaN;
        };

        const LocalGeodesicInfo& calcInitialGeodesic(State& s, const GeodesicInitialConditions& g0) const;
        // TODO weird name...
        const LocalGeodesicInfo& calcLocalGeodesic(const State& s, Vec3 prev_QS, Vec3 next_PS) const;
        void applyGeodesicCorrection(const State& s, const Correction& c) const;

        // 4. calcPathPoints
        size_t calcPathPoints(const State& state, std::vector<Vec3>& points) const;

        void setInitialPointGuess(Vec3 initPointGuess) {m_InitPointGuess = initPointGuess;}
        Vec3 getInitialPointGuess() const {return m_InitPointGuess;}

        private:

        struct CacheEntry : LocalGeodesicInfo
        {
            Vec3 trackingPointOnLine {NaN, NaN, NaN};
            std::vector<Vec3> points {};
            double sHint = NaN;
        };

        const CacheEntry& getCacheEntry(const State& state) const
        {
            return Value<CacheEntry>::downcast(
                    m_Subsystem.getDiscreteVarUpdateValue(state, m_CacheIx)
                    );
        }

        CacheEntry& updCacheEntry(const State& state) const
        {
            return Value<CacheEntry>::updDowncast(
                    m_Subsystem.updDiscreteVarUpdateValue(state, m_CacheIx)
                    );
        }

        const CacheEntry& getPrevCacheEntry(const State& state) const
        {
            return Value<CacheEntry>::downcast(
                    m_Subsystem.getDiscreteVariable(state, m_CacheIx)
                    );
        }

        CacheEntry& updPrevCacheEntry(State& state) const
        {
            return Value<CacheEntry>::updDowncast(
                    m_Subsystem.updDiscreteVariable(state, m_CacheIx)
                    );
        }

        void calcStatus(const Vec3& prev_QS, const Vec3& next_PS, CacheEntry& cache) const;
        void calcGeodesic(const GeodesicInitialConditions& g0, CacheEntry& cache) const;

        //------------------------------------------------------------------------------
        WrappingPathSubsystem m_Subsystem;

        ContactGeometry m_Geometry;

        Vec3 m_InitPointGuess;

        Real m_TouchdownAccuracy = 1e-3;
        size_t m_TouchdownIter = 10;

        DiscreteVariableIndex m_CacheIx;
    };

    private:

    // TODO Required for accessing the cache variable?
    WrappingPathSubsystem m_Subsystem; // The subsystem this segment belongs to.
    WrappingPath m_Path; // The path this segment belongs to.
    PathSegmentIndex m_PathSegmentIx; // The index in its path.

    MobilizedBody m_Mobod;
    Transform m_Offset;

    Surface m_Surface;

    // TOPOLOGY CACHE
    CacheEntryIndex       m_PosInfoIx;
};


//==============================================================================
//                         PATH :: IMPL
//==============================================================================
class WrappingPath::Impl {
    public:
        Impl(
                WrappingPathSubsystem subsystem,
                MobilizedBody originBody,
                Vec3 originPoint,
                MobilizedBody terminationBody,
                Vec3 terminationPoint):
            m_Subsystem(subsystem),
            m_OriginBody(originBody),
            m_OriginPoint(originPoint),
            m_TerminationBody(terminationBody),
            m_TerminationPoint(terminationPoint) {}

        struct PosInfo
        {
            Vec3 xO {NaN, NaN, NaN};
            Vec3 xI {NaN, NaN, NaN};

            Real l = NaN;

            std::vector<LineSegment> lines;

            Vector pathError {};
            Matrix pathMatrix {};
            Matrix pathErrorJacobian {};
            Vector pathCorrections {};

            size_t loopIter = 0;

            // TODO solver matrices
        };

        std::vector<WrapObstacle>& updObstacles() {return m_Obstacles;}
        const std::vector<WrapObstacle>& getObstacles() const {return m_Obstacles;}

        // Allocate state variables and cache entries.
        void realizeTopology(State& state);
        void realizePosition(const State& state) const;
        void realizeVelocity(const State& state) const;
        void realizeAcceleration(const State& state) const;
        void invalidateTopology()
        {   m_Subsystem.invalidateSubsystemTopologyCache(); }

        const PosInfo& getPosInfo(const State& state) const;

        void calcInitZeroLengthGeodesic(State& state, std::function<Vec3(size_t)> GetInitPointGuess) const;


    private:
        PosInfo& updPosInfo(const State& s) const;
        void calcPosInfo(const State& s, PosInfo& posInfo) const;

        size_t countActive(const State& s) const;

        Vec3 findPrevPoint(
                const State& state,
                PathSegmentIndex idx) const;

        Vec3 findNextPoint(
                const State& state,
                PathSegmentIndex idx) const;

        static Vec3 FindPrevPoint(
                const State& state,
                const Vec3& originPoint,
                const std::vector<WrapObstacle>& obs,
                size_t idx);

        static Vec3 FindNextPoint(
                const State& state,
                const Vec3& terminationPoint,
                const std::vector<WrapObstacle>& obs,
                size_t idx);

        WrapObstacle* findPrevActiveObstacle(const State& s, size_t obsIdx);
        WrapObstacle* findNextActiveObstacle(const State& s, size_t obsIdx);

        template<size_t N>
            static void calcPathErrorVector(
                    const State& state,
                    const std::vector<WrapObstacle>& obs,
                    const std::vector<LineSegment>& lines,
                    std::array<CoordinateAxis, N> axes,
                    Vector& pathError);

        template<size_t N>
            static void calcPathErrorJacobian(
                    const State& state,
                    const std::vector<WrapObstacle>& obs,
                    const std::vector<LineSegment>& lines,
                    std::array<CoordinateAxis, N> axes,
                    Matrix& J);

        static void calcLineSegments(
                const State& s,
                Vec3 p_O,
                Vec3 p_I,
                const std::vector<WrapObstacle>& obs,
                std::vector<LineSegment>& lines);

        static size_t CalcUpdatedObstaclePosInfo(
                const State& s,
                const Vec3& x_O,
                const Vec3& x_I,
                const std::vector<WrapObstacle>& obs,
                size_t maxIter,
                Real eps);

static double calcPathLength(
	const State& state,
    const std::vector<WrapObstacle>& obs,
    const std::vector<LineSegment>& lines);

        WrappingPathSubsystem m_Subsystem;

        MobilizedBody m_OriginBody;
        Vec3 m_OriginPoint;

        MobilizedBody m_TerminationBody;
        Vec3 m_TerminationPoint;

        std::vector<WrapObstacle> m_Obstacles {};

        Real m_PathErrorBound = 0.1;
        Real m_ObsErrorBound = 0.1;
        size_t m_PathMaxIter = 10;
        size_t m_ObsMaxIter = 10;

        // TOPOLOGY CACHE (set during realizeTopology())
        CacheEntryIndex       m_PosInfoIx;
        CacheEntryIndex       m_VelInfoIx;
        CacheEntryIndex       m_VizInfoIx;

        friend WrapObstacle::Impl;
};

//==============================================================================
//                         SUBSYSTEM :: IMPL
//==============================================================================
class WrappingPathSubsystem::Impl : public Subsystem::Guts {
    public:
        Impl() {}
        ~Impl() {}
        Impl* cloneImpl() const override 
        {   return new Impl(*this); }

        int getNumPaths() const {return paths.size();}

        const WrappingPath& getCablePath(WrappingPathIndex index) const 
        {   return paths[index]; }

        WrappingPath& updCablePath(WrappingPathIndex index) 
        {   return paths[index]; }

        // Add a cable path to the list, bumping the reference count.
        WrappingPathIndex adoptCablePath(WrappingPath& path) {
            paths.push_back(path);
            return WrappingPathIndex(paths.size()-1);
        }

        // Return the MultibodySystem which owns this WrappingPathSubsystem.
        const MultibodySystem& getMultibodySystem() const 
        {   return MultibodySystem::downcast(getSystem()); }

        // Return the SimbodyMatterSubsystem from which this WrappingPathSubsystem
        // gets the bodies to track.
        const SimbodyMatterSubsystem& getMatterSubsystem() const 
        {   return getMultibodySystem().getMatterSubsystem(); }

        // Allocate state variables.
        int realizeSubsystemTopologyImpl(State& state) const override {
            // Briefly allow writing into the Topology cache; after this the
            // Topology cache is const.
            Impl* wThis = const_cast<Impl*>(this);

            for (WrappingPathIndex ix(0); ix < paths.size(); ++ix) {
                WrappingPath& path = wThis->updCablePath(ix);
                path.updImpl().realizeTopology(state);
            }

            return 0;
        }

        int realizeSubsystemPositionImpl(const State& state) const override {
            for (WrappingPathIndex ix(0); ix < paths.size(); ++ix) {
                const WrappingPath& path = getCablePath(ix);
                path.getImpl().realizePosition(state);
            }
            return 0;
        }

        int realizeSubsystemVelocityImpl(const State& state) const override {
            for (WrappingPathIndex ix(0); ix < paths.size(); ++ix) {
                const WrappingPath& path = getCablePath(ix);
                path.getImpl().realizeVelocity(state);
            }
            return 0;
        }


        int realizeSubsystemAccelerationImpl(const State& state) const override {
            for (WrappingPathIndex ix(0); ix < paths.size(); ++ix) {
                const WrappingPath& path = getCablePath(ix);
                path.getImpl().realizeAcceleration(state);
            }
            return 0;
        }

        SimTK_DOWNCAST(Impl, Subsystem::Guts);

    private:
        // TOPOLOGY STATE
        Array_<WrappingPath, WrappingPathIndex> paths;
};

} // namespace SimTK

#endif
