#include "SimTKcommon/internal/CoordinateAxis.h"
#include "SimTKcommon/internal/ExceptionMacros.h"
#include "SimTKmath.h"
#include "Wrapping.h"
#include "WrappingImpl.h"
#include "simmath/internal/ContactGeometry.h"
#include <cstddef>
#include <stdexcept>

using namespace SimTK;

using GeodesicInfo = WrapObstacle::Impl::PosInfo;

//==============================================================================
//                                CONSTANTS
//==============================================================================
namespace {
    static const CoordinateAxis TangentAxis = ContactGeometry::TangentAxis;
    static const CoordinateAxis NormalAxis = ContactGeometry::NormalAxis;
    static const CoordinateAxis BinormalAxis = ContactGeometry::BinormalAxis;
    static const int GeodesicDOF = 4;
}

//==============================================================================
//                                SOLVER
//==============================================================================

namespace {
    struct SolverData;

    struct SolverData 
    {
        // A * x = b
        Matrix A;
        Vector x;
        Vector b;
    };

    class SolverDataCache
    {
        public:
        SolverDataCache(std::mutex& cacheMutex) :
            m_guard(cacheMutex) {}

        void setDataPtr(SolverData* ptr) {m_data = ptr;}

        SolverData& updData() {return *m_data;}

        private:
        std::lock_guard<std::mutex> m_guard;
        SolverData* m_data = nullptr;
    };

    SolverDataCache&& findDataCache(size_t nActive)
    {
        static constexpr int Q = 4;
        static constexpr int C = 4;
        static std::vector<SolverData> s_GlobalCache {};
        static std::mutex s_GlobalLock {};

        {
            std::lock_guard<std::mutex> lock(s_GlobalLock);

            for (size_t i = s_GlobalCache.size(); i < nActive - 1; ++i) {
                int n = i + 1;
                s_GlobalCache.emplace_back(
                        Matrix{C * n, Q * n, 0.},
                        Vector{Q * n, 0.},
                        Vector{C * n, 0.});
            }
        }

        return SolverDataCache(s_GlobalLock);
    }
} // namespace

//==============================================================================
//                                SOME MATH
//==============================================================================

namespace {
    using FrenetFrame = ContactGeometry::FrenetFrame;
    using Variation = ContactGeometry::GeodesicVariation;
    using Correction = ContactGeometry::GeodesicCorrection;
    using GeodesicInitialConditions = Surface::Impl::GeodesicInitialConditions;

    GeodesicInitialConditions calcCorrectedGeodesicInitConditions(const FrenetFrame& KP, const Variation& dKP, Real l, const Correction& c)
    {
        Surface::Impl::GeodesicInitialConditions g0;

        Vec3 v = dKP[1] * c;
        g0.x = KP.p() + v;

        Vec3 w = dKP[0] * c;
        const UnitVec3 t = KP.R().getAxisUnitVec(ContactGeometry::TangentAxis);
        g0.t = t + cross(w,t);

        Real dl = c[3];
        g0.l = l + dl;

        return g0;
    }
}

//==============================================================================
//                                SURFACE IMPL
//==============================================================================

//------------------------------------------------------------------------------
//                               REALIZE CACHE
//------------------------------------------------------------------------------
void Surface::Impl::realizeTopology(State &state)
{
	// Allocate an auto-update discrete variable for the last computed geodesic.
	LocalGeodesicInfo geodesic {};
	m_GeodesicInfoIx = m_Subsystem.allocateAutoUpdateDiscreteVariable(state, Stage::Velocity, new Value<LocalGeodesicInfo>(geodesic), Stage::Position);
}

void Surface::Impl::realizeInstance(const State &state) const
{
	throw std::runtime_error("NOTYETIMPLEMENTED");
}

void Surface::Impl::realizePosition(const State &state) const
{
	// Set the current local geodesic equal to the previous.
	if (!m_Subsystem.isDiscreteVarUpdateValueRealized(state, m_GeodesicInfoIx)) {
		updLocalGeodesicInfo(state) = getPrevLocalGeodesicInfo(state);
		m_Subsystem.markDiscreteVarUpdateValueRealized(state, m_GeodesicInfoIx);
	}
}

void Surface::Impl::realizeVelocity(const State &state) const
{
	throw std::runtime_error("NOTYETIMPLEMENTED");
}

void Surface::Impl::realizeAcceleration(const State &state) const
{
	throw std::runtime_error("NOTYETIMPLEMENTED");
}

//------------------------------------------------------------------------------
//                               PUBLIC METHODS
//------------------------------------------------------------------------------
void Surface::Impl::applyGeodesicCorrection(const State& state, const Correction& c) const
{
	// Get the previous geodesic.
	const LocalGeodesicInfo& g = getLocalGeodesicInfo(state);

    // Get corrected initial conditions.
    const GeodesicInitialConditions g0 = calcCorrectedGeodesicInitConditions(g.KP, g.dKP, g.length, c);

	// Shoot the new geodesic.
	calcLocalGeodesicInfo(g0.x, g0.t, g0.l, g.sHint, updLocalGeodesicInfo(state));
}

const Surface::WrappingStatus& Surface::Impl::calcWrappingStatus(
        const State& state,
        Vec3 prevPoint,
        Vec3 nextPoint, size_t maxIter, Real eps) const
{
    LocalGeodesicInfo& g = updLocalGeodesicInfo(state);

    if (g.disabled) {return g;}

    // Make sure that the previous point does not lie inside the surface.
    if (m_Geometry.calcSurfaceValue(prevPoint)) {
        // TODO use proper assert.
        throw std::runtime_error("Unable to wrap over surface: Preceding point lies inside the surface");
    }
    if (m_Geometry.calcSurfaceValue(nextPoint)) {
        // TODO use proper assert.
        throw std::runtime_error("Unable to wrap over surface: Next point lies inside the surface");
    }

    // Update the tracking point.
    ContactGeometry::PointOnLineResult result = m_Geometry.calcNearestPointOnLine(prevPoint, nextPoint, g.pointOnLine, maxIter, eps);
    g.pointOnLine = result.p; // TODO rename to point.

    // Attempt touchdown.
    if (g.liftoff) {
        const bool touchdown = result.isInsideSurface;
        g.liftoff = !touchdown;
    }

    // Detect liftoff.
    if (!g.liftoff) {
        g.liftoff = g.length == 0.;
        g.liftoff &= dot(prevPoint - g.KP.p(), g.KP.R().getAxisUnitVec(NormalAxis)) <= 0.;
        g.liftoff &= dot(nextPoint - g.KP.p(), g.KP.R().getAxisUnitVec(NormalAxis)) <= 0.;
    }

    return g;
}

size_t Surface::Impl::calcPathPoints(const State& state, std::vector<Vec3>& points) const
{
    size_t count = 0;
    const LocalGeodesicInfo& g = getLocalGeodesicInfo(state);
    for (Vec3 p: g.points) {
        points.push_back(p);
    }
    return count;
}

void Surface::Impl::setInitialPointGuess(Vec3 pointGuess) {
	throw std::runtime_error("NOTYETIMPLEMENTED");
}

Vec3 Surface::Impl::getInitialPointGuess() const {
	throw std::runtime_error("NOTYETIMPLEMENTED");
}

//------------------------------------------------------------------------------
//                               PRIVATE METHODS
//------------------------------------------------------------------------------
void Surface::Impl::calcLocalGeodesicInfo(Vec3 x, Vec3 t, Real l, Real sHint,
        LocalGeodesicInfo& geodesic) const
{
	// Compute geodesic start boundary frame and variation.
	m_Geometry.calcNearestFrenetFrameFast(x, t, geodesic.KP);
	m_Geometry.calcGeodesicStartFrameVariation(geodesic.KP, geodesic.dKP);

    // Compute geodesic end boundary frame amd variation (shoot new geodesic).
	m_Geometry.calcGeodesicEndFrameVariationImplicitly(
			geodesic.KP.p(),
			geodesic.KP.R().getAxisUnitVec(ContactGeometry::TangentAxis),
			l,
			sHint,
			geodesic.KQ,
			geodesic.dKQ,
			geodesic.points);

	// TODO  update step size.
	// TODO  update line tracking?
	throw std::runtime_error("NOTYETIMPLEMENTED");
}

//==============================================================================
//                               OBSTACLE
//==============================================================================

void WrapObstacle::Impl::realizeTopology(State &state)
{
	// Allocate position level cache.
	PosInfo posInfo {};
	m_PosInfoIx = m_Subsystem.allocateCacheEntry(state, Stage::Position, new Value<PosInfo>(posInfo));
}

void WrapObstacle::Impl::realizePosition(const State &state) const
{
	if (!m_Subsystem.isCacheValueRealized(state, m_PosInfoIx)) {
		calcPosInfo(state, updPosInfo(state));
		m_Subsystem.markCacheValueRealized(state, m_PosInfoIx);
	}
}

const WrapObstacle::Impl::PosInfo& WrapObstacle::Impl::getPosInfo(const State &state) const
{
    realizePosition(state);
    return Value<PosInfo>::downcast(m_Subsystem.getCacheEntry(state, m_PosInfoIx));
}

void WrapObstacle::Impl::applyGeodesicCorrection(const State& state, const WrapObstacle::Impl::Correction& c) const
{
    // Apply correction to curve.
    m_Surface.getImpl().applyGeodesicCorrection(state, c);

	// Invalidate position level cache.
	m_Subsystem.markCacheValueNotRealized(state, m_PosInfoIx);
}

void WrapObstacle::Impl::calcPosInfo(
		const State& state,
		PosInfo& posInfo) const
{
    if(isDisabled(state))
    {
        return;
    }

	// Get tramsform from local surface frame to ground.
	Transform X_GS = m_Mobod.getBodyTransform(state).compose(m_Offset);

    // Get the path points before and after this segment.
    Vec3 prev_G = m_Path.getImpl().findPrevPoint(state, m_PathSegmentIx);
    Vec3 next_G = m_Path.getImpl().findNextPoint(state, m_PathSegmentIx);

    // Transform the prev and next path points to the surface frame.
    Vec3 prev_S = X_GS.shiftBaseStationToFrame(prev_G);
    Vec3 next_S = X_GS.shiftBaseStationToFrame(next_G);

    // Detect liftoff, touchdown and potential invalid configurations.
    // TODO this doesnt follow the regular invalidation scheme...
    m_Surface.getImpl().calcUpdatedStatus(state, prev_S, next_S);

	// Grab the last geodesic that was computed.
	const Surface::LocalGeodesic& geodesic_S = m_Surface.getImpl().getGeodesic(state);

    // Store the the local geodesic in ground frame.
    posInfo.KP = X_GS.compose(geodesic_S.KP);
    posInfo.KQ = X_GS.compose(geodesic_S.KQ);

    posInfo.dKP[0] = X_GS.R() * geodesic_S.dKP[0];
    posInfo.dKP[1] = X_GS.R() * geodesic_S.dKP[1];

    posInfo.dKQ[0] = X_GS.R() * geodesic_S.dKQ[0];
    posInfo.dKQ[1] = X_GS.R() * geodesic_S.dKQ[1];

    // TODO use SpatialVec for variation.
    /* for (size_t i = 0; i < GeodesicDOF; ++i) { */
    /*     posInfo.dKP[i][0] = X_GS.xformFrameVecToBase(geodesic_S.dKP[i][0]) */
    /*     posInfo.dKP[i][1] = X_GS.xformFrameVecToBase(geodesic_S.dKP[i][1]) */

    /*     posInfo.dKQ[i][0] = X_GS.xformFrameVecToBase(geodesic_S.dKQ[i][0]) */
    /*     posInfo.dKQ[i][1] = X_GS.xformFrameVecToBase(geodesic_S.dKQ[i][1]) */
    /* } */

    posInfo.length = geodesic_S.length;

    throw std::runtime_error("NOTYETIMPLEMENTED: Check transformation order");
}

void WrapObstacle::Impl::calcGeodesic(State& state, Vec3 x, Vec3 t, Real l) const
{
    m_Surface.getImpl().calcGeodesic(state, x, t, l);
	m_Subsystem.markCacheValueNotRealized(state, m_PosInfoIx);
}

//==============================================================================
//                                PATH HELPERS
//==============================================================================

namespace
{
    using GeodesicJacobian = Vec4;
    using LineSegment = WrappingPath::LineSegment;
    using PointVariation = ContactGeometry::GeodesicPointVariation;
    using FrameVariation = ContactGeometry::GeodesicFrameVariation;
    using Variation = ContactGeometry::GeodesicVariation;
    using Correction = ContactGeometry::GeodesicCorrection;
/* using LocalGeodesic = WrapObstacle::LocalGeodesicInfo; */

int countActive(const State& s, const std::vector<WrapObstacle>& obstacles)
{
    int n = 0;
    for (const WrapObstacle& o : obstacles) {
        if (o.isActive(s)) {
            ++n;
        }
    }
    return n;
}

const WrapObstacle* FindPrevActiveObstacle(
	const State& state,
    const std::vector<WrapObstacle>& obs,
    size_t idx)
{
    for (ptrdiff_t i = idx - 1; i > 0; --i) {
        // Find the active segment before the current.
        if (obs.at(i).isActive(state)) {
            return &obs.at(i);
        }
    }
    return nullptr;
}

const WrapObstacle* FindNextActiveObstacle(
	const State& state,
    const std::vector<WrapObstacle>& obs,
    size_t idx)
{
    // Find the active segment after the current.
    for (ptrdiff_t i = idx + 1; i < obs.size(); ++i) {
        if (obs.at(i).isActive(state)) {
            return &obs.at(i);
        }
    }
    return nullptr;
}

static const int N_PATH_CONSTRAINTS = 4;

// TODO this is awkward
void addBlock(
        const Vec4& values,
    Matrix& block)
{
    for (int i = 0; i < Vec4::size(); ++i) {
        block[0][i] = values[i];
    }
}

void addDirectionJacobian(
    const LineSegment& e,
	const UnitVec3& axis,
    const PointVariation& dx,
    Matrix& J,
    bool invert = false)
{
    Vec3 y = axis - e.d * dot(e.d,axis);
    y /= e.l * (invert ? 1. : -1);
    addBlock(~dx * y, J);
}

double calcPathError(const LineSegment& e, const Rotation& R, CoordinateAxis axis)
{
    return dot(e.d, R.getAxisUnitVec(axis));
}

GeodesicJacobian addPathErrorJacobian(
    const LineSegment& e,
	const UnitVec3& axis,
    const Variation& dK,
	Matrix& J,
    bool invertV = false)
{
	addDirectionJacobian(e, axis, dK[1], J, invertV);
	addBlock(~dK[0] * cross(axis,e.d), J);
}

}

//==============================================================================
//                                PATH IMPL
//==============================================================================

void WrappingPath::Impl::realizeTopology(State &state)
{
	PosInfo posInfo {};
	m_PosInfoIx = m_Subsystem.allocateCacheEntry(state, Stage::Position, new Value<PosInfo>(posInfo));
}

void WrappingPath::Impl::realizePosition(const State &state) const
{
	if (m_Subsystem.isCacheValueRealized(state, m_PosInfoIx)) {return;}
	calcPosInfo(state, updPosInfo(state));
	m_Subsystem.markCacheValueRealized(state, m_PosInfoIx);
}

const WrappingPath::Impl::PosInfo& WrappingPath::Impl::getPosInfo(const State &state) const
{
	realizePosition(state);
    return Value<PosInfo>::downcast(m_Subsystem.getCacheEntry(state, m_PosInfoIx));
}

WrappingPath::Impl::PosInfo& WrappingPath::Impl::updPosInfo(const State &state) const
{
    return Value<PosInfo>::updDowncast(m_Subsystem.updCacheEntry(state, m_PosInfoIx));
}

void WrappingPath::Impl::calcInitZeroLengthGeodesic(State& state, std::function<Vec3(size_t)> GetInitPointGuess) const
{
	Vec3 prev = m_OriginBody.getBodyTransform(state).shiftFrameStationToBase(m_OriginPoint);
	size_t i = 0;
	for (const WrapObstacle& obstacle: m_Obstacles) {
	    Vec3 xGuess = GetInitPointGuess(++i);
	    obstacle.getImpl().calcGeodesic(state, xGuess, xGuess - prev, 0.);

		prev = obstacle.getImpl().getPosInfo(state).KQ.p();
	}
}

template<size_t N>
void WrappingPath::Impl::calcPathErrorVector(
	const State& state,
    const std::vector<WrapObstacle>& obs,
    const std::vector<LineSegment>& lines,
    std::array<CoordinateAxis, N> axes,
    Vector& pathError)
{
    size_t lineIx = 0;
    ptrdiff_t row  = -1;

    for (const WrapObstacle& o : obs) {
        if (!o.isActive(state)) {
            continue;
        }

        const GeodesicInfo& g = o.getImpl().getPosInfo(state);
        for (CoordinateAxis axis: axes) {
            pathError(++row)  = calcPathError(lines.at(lineIx), g.KP.R(), axis);
        }
        ++lineIx;
        for (CoordinateAxis axis: axes) {
            pathError(++row)  = calcPathError(lines.at(lineIx), g.KQ.R(), axis);
        }
    }
}

template<size_t N>
void WrappingPath::Impl::calcPathErrorJacobian(
	const State& state,
    const std::vector<WrapObstacle>& obs,
    const std::vector<LineSegment>& lines,
    std::array<CoordinateAxis, N> axes,
    Matrix& J)
{
    constexpr size_t Nq = GeodesicDOF;
    const size_t n     = countActive(state, obs);

    SimTK_ASSERT(J.rows() == n * N, "Invalid number of rows in jacobian matrix");
    SimTK_ASSERT(J.cols() == n * Nq, "Invalid number of columns in jacobian matrix");

    size_t row = 0;
    size_t col = 0;
    for (size_t i = 0; i < n; ++i) {
        if (!obs.at(i).isActive(state)) {
            continue;
        }
        const GeodesicInfo& g = obs.at(i).getImpl().getPosInfo(state);

        const LineSegment& l_P = lines.at(i);
        const LineSegment& l_Q = lines.at(i + 1);

        const WrapObstacle* prev = FindPrevActiveObstacle(state, obs, i);
        const WrapObstacle* next = FindNextActiveObstacle(state, obs, i);

        for (CoordinateAxis axis: axes) {
            const UnitVec3 a_P = g.KP.R().getAxisUnitVec(axis);
            const Variation& dK_P = g.dKP;

			addPathErrorJacobian(l_P, a_P, dK_P, J.block(row, col, 1, Nq));

			if (prev) {
				const Variation& prev_dK_Q = prev->getImpl().getPosInfo(state).dKQ;
				addDirectionJacobian(l_P, a_P, prev_dK_Q[1], J.block(row, col-Nq, 1, Nq), true);
			}
			++row;
		}

        for (CoordinateAxis axis: axes) {
            const UnitVec3 a_Q = g.KQ.R().getAxisUnitVec(axis);
            const Variation& dK_Q = g.dKQ;

			addPathErrorJacobian(l_Q, a_Q, dK_Q, J.block(row, col, 1, Nq), true);

			if (next) {
				const Variation& next_dK_P = next->getImpl().getPosInfo(state).dKP;
				addDirectionJacobian(l_Q, a_Q, next_dK_P[1], J.block(row, col+Nq, 1, Nq));
			}
			++row;
		}

        col += Nq;
    };
}

double WrappingPath::Impl::calcPathLength(
	const State& state,
    const std::vector<WrapObstacle>& obs,
    const std::vector<LineSegment>& lines)
{
    double lTot = 0.;
    for (const LineSegment& line : lines) {
        // TODO spell out as length.
        lTot += line.l;
    }

    for (const WrapObstacle& obstacle : obs) {
        if (!obstacle.isActive(state))
		{
            continue;
		}
        lTot += obstacle.getImpl().getPosInfo(state).length;
    }
    return lTot;
}

void WrappingPath::Impl::calcLineSegments(
	const State& state,
    Vec3 p_O,
    Vec3 p_I,
    const std::vector<WrapObstacle>& obs,
    std::vector<LineSegment>& lines)
{
    const size_t n = obs.size();
    lines.reserve(n + 1);
    lines.clear();

    Vec3 lineStart = std::move(p_O);
    for (const WrapObstacle& o : obs) {
        if (!o.isActive(state)) {
            continue;
        }

        const GeodesicInfo& g = o.getImpl().getPosInfo(state);
        const Vec3 lineEnd = g.KP.p();
        lines.emplace_back(lineStart, lineEnd);

        lineStart = g.KQ.p();
    }
    lines.emplace_back(lineStart, p_I);
}

Vec3 WrappingPath::Impl::FindPrevPoint(
        const State& state,
        const Vec3& originPoint,
    const std::vector<WrapObstacle>& obs,
    size_t idx)
{
    const WrapObstacle* prev = FindPrevActiveObstacle(state, obs, idx);
    return prev ? prev->getImpl().getPosInfo(state).KQ.p() : originPoint;
}

Vec3 WrappingPath::Impl::FindNextPoint(
        const State& state,
        const Vec3& terminationPoint,
    const std::vector<WrapObstacle>& obs,
    size_t idx)
{
    const WrapObstacle* next = FindNextActiveObstacle(state, obs, idx);
    return next ? next->getImpl().getPosInfo(state).KP.p() : terminationPoint;
}

void WrappingPath::Impl::calcPosInfo(const State& s, PosInfo& posInfo) const
{
	// Path origin and termination point.
	const Vec3 x_O = m_OriginBody.getBodyTransform(s).shiftFrameStationToBase(m_OriginPoint);
	const Vec3 x_I = m_TerminationBody.getBodyTransform(s).shiftFrameStationToBase(m_TerminationPoint);

	const std::array<CoordinateAxis, 2> axes {NormalAxis, BinormalAxis};

    for (posInfo.loopIter = 0; posInfo.loopIter < m_PathMaxIter; ++posInfo.loopIter) {
        // Compute the straight-line segments.
        calcLineSegments(s, x_O, x_I, m_Obstacles, posInfo.lines);

        // Evaluate path error, and stop when converged.
        calcPathErrorVector<2>(s, m_Obstacles, posInfo.lines, axes, posInfo.pathError);
        const Real maxPathError = posInfo.pathError.normInf();
        if (maxPathError < m_PathErrorBound) {
            return;
        }

        // Evaluate the path error jacobian.
        calcPathErrorJacobian<2>(s, m_Obstacles, posInfo.lines, axes, posInfo.pathErrorJacobian);

        // Compute path corrections.
        calcPathCorrections(s, m_Obstacles, posInfo.pathError, posInfo.pathErrorJacobian, posInfo.pathMatrix, posInfo.pathCorrections);

        // Apply path corrections.
		throw std::runtime_error("NOTYETIMPLEMENTED");
		const Correction* corrIt = nullptr; // TODO
        for (const WrapObstacle& obstacle : m_Obstacles) {
            if (!obstacle.isActive(s)) {
                continue;
            }
            obstacle.getImpl().applyGeodesicCorrection(s, *corrIt);
            ++corrIt;
        }

        // Path has changed: invalidate each segment's cache.
        for (const WrapObstacle& obstacle : m_Obstacles) {
            obstacle.getImpl().invalidatePositionLevelCache(s);
        }
    }

	throw std::runtime_error("Failed to converge");
}

//==============================================================================
//                                WRAP OBSTACLE IMPL
//==============================================================================

namespace
{

void calcPathCorrections(const State& s, const std::vector<WrapObstacle>& obs, const Vector& pathError, const Matrix& pathErrorJacobian, Matrix& pathMatrix, Vector& pathCorrections);

}

//==============================================================================
//            SUBSYSTEM
//==============================================================================

bool WrappingPathSubsystem::isInstanceOf(const Subsystem& s) {
    return Impl::isA(s.getSubsystemGuts());
}

const WrappingPathSubsystem& WrappingPathSubsystem::
downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return static_cast<const WrappingPathSubsystem&>(s);
}
WrappingPathSubsystem& WrappingPathSubsystem::
updDowncast(Subsystem& s) {
    assert(isInstanceOf(s));
    return static_cast<WrappingPathSubsystem&>(s);
}

const WrappingPathSubsystem::Impl& WrappingPathSubsystem::
getImpl() const {
    return SimTK_DYNAMIC_CAST_DEBUG<const Impl&>(getSubsystemGuts());
}
WrappingPathSubsystem::Impl& WrappingPathSubsystem::
updImpl() {
    return SimTK_DYNAMIC_CAST_DEBUG<Impl&>(updSubsystemGuts());
}

// Create Subsystem but don't associate it with any System. This isn't much use
// except for making std::vectors, which require a default constructor to be 
// available.
WrappingPathSubsystem::WrappingPathSubsystem() 
{   adoptSubsystemGuts(new Impl()); }

WrappingPathSubsystem::WrappingPathSubsystem(MultibodySystem& mbs) 
{   adoptSubsystemGuts(new Impl());
    mbs.adoptSubsystem(*this); } // steal ownership

int WrappingPathSubsystem::getNumPaths() const
{   return getImpl().getNumPaths(); }

const WrappingPath& WrappingPathSubsystem::
getPath(WrappingPathIndex cableIx) const
{   return getImpl().getCablePath(cableIx); }

WrappingPath& WrappingPathSubsystem::
updPath(WrappingPathIndex cableIx)
{   return updImpl().updCablePath(cableIx); }
