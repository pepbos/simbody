#include "SimTKcommon/internal/CoordinateAxis.h"
#include "SimTKcommon/internal/ExceptionMacros.h"
#include "SimTKmath.h"
#include "Wrapping.h"
#include "WrappingImpl.h"
#include "simmath/internal/ContactGeometry.h"
#include <cstddef>
#include <stdexcept>

using namespace SimTK;

using Correction = ContactGeometry::GeodesicCorrection;
using FrenetFrame = ContactGeometry::FrenetFrame;
using GeodesicInitialConditions = WrapObstacle::Impl::Surface::GeodesicInitialConditions;
using LocalGeodesicInfo = WrapObstacle::Impl::Surface::LocalGeodesicInfo;
using GeodesicInfo = WrapObstacle::Impl::GeodesicInfo;
using Variation = ContactGeometry::GeodesicVariation;

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
    class SolverDataCache;
    SolverDataCache& findDataCache(size_t nActive);

    struct SolverData 
    {
        std::vector<LineSegment> lineSegments;

        Matrix pathErrorJacobian;
        Vector pathCorrection;
        Vector pathError;
        Matrix mat;
        // TODO Cholesky decomposition...
        FactorLU matInv;
        Vector vec;
    };

    class SolverDataCache
    {
        private:
        SolverDataCache(size_t n) {
            static constexpr int Q = 4;
            static constexpr int C = 4;

            m_Data.lineSegments.resize(n+1);
            m_Data.pathErrorJacobian = Matrix(C * n, Q * n, 0.);
            m_Data.pathCorrection    = Vector(Q * n, 0.);
            m_Data.pathError         = Vector(C * n, 0.);
            m_Data.mat               = Matrix(Q*n, Q*n, NaN);
            m_Data.vec               = Vector(Q*n, NaN);
        }

        public:
        void lock() {
            m_Mutex.lock();
        }

        void unlock() {
            m_Mutex.unlock();
        }

        SolverData& updData() {return m_Data;}

        private:
        SolverData m_Data;
        std::mutex m_Mutex {};

        friend SolverDataCache& findDataCache(size_t nActive);
    };

    SolverDataCache& findDataCache(size_t nActive)
    {
        static std::vector<SolverDataCache> s_GlobalCache {};
        static std::mutex s_GlobalLock {};

        {
            std::lock_guard<std::mutex> lock(s_GlobalLock);

            for (size_t i = s_GlobalCache.size(); i < nActive - 1; ++i) {
                int n = i + 1;
                s_GlobalCache.emplace_back(nActive);
            }
        }

        return s_GlobalCache.at(nActive-1);
    }

    const Correction* calcPathCorrections(SolverData& data) {
        Real w = data.pathError.normInf();

        data.mat = data.pathErrorJacobian.transpose() * data.pathErrorJacobian;
        for (int i = 0; i < data.mat.nrow(); ++i) {
            data.mat[i][i] += w;
        }
        data.matInv = data.mat;
        data.vec = data.pathErrorJacobian.transpose() * data.pathError;
        data.matInv.solve(data.vec, data.pathCorrection);

        static_assert(
                sizeof(Correction) == sizeof(double) * GeodesicDOF,
                "Invalid size of corrections vector");
        SimTK_ASSERT(data.pathCorrection.size() * sizeof(double) == n * sizeof(Correction),
                "Invalid size of path corrections vector");
        return reinterpret_cast<const Correction*>(&data.pathCorrection[0]);
    }
} // namespace

//==============================================================================
//                          GEODESIC INITIAL CONDITIONS
//==============================================================================

GeodesicInitialConditions GeodesicInitialConditions::CreateCorrected(const FrenetFrame& KP, const Variation& dKP, Real l, const Correction& c)
{
    GeodesicInitialConditions g0;

    Vec3 v = dKP[1] * c;
    g0.x = KP.p() + v;

    Vec3 w = dKP[0] * c;
    const UnitVec3 t = KP.R().getAxisUnitVec(ContactGeometry::TangentAxis);
    g0.t = t + cross(w,t);

    Real dl = c[3];
    g0.l = l + dl;

    return g0;
}

GeodesicInitialConditions GeodesicInitialConditions::CreateInSurfaceFrame(const Transform& X_GS, Vec3 x_G, Vec3 t_G, Real l)
{
    GeodesicInitialConditions g0;

    g0.x = X_GS.shiftBaseStationToFrame(x_G);
    g0.t = X_GS.xformBaseVecToFrame(t_G);
    g0.l = l;

    return g0;
}

GeodesicInitialConditions GeodesicInitialConditions::CreateZeroLengthGuess(const Transform& X_GS, Vec3 prev_QS, Vec3 xGuess_S)
{
    GeodesicInitialConditions g0;

    g0.x = xGuess_S;
    g0.t = g0.x - X_GS.shiftBaseStationToFrame(prev_QS);
    g0.l = 0.;

    return g0;
}

GeodesicInitialConditions GeodesicInitialConditions::CreateAtTouchdown(Vec3 prev_QS, Vec3 next_PS, Vec3 trackingPointOnLine)
{
    GeodesicInitialConditions g0;

    g0.x = trackingPointOnLine;
    g0.t = next_PS - prev_QS;
    g0.l = 0.;

    return g0;
}

//==============================================================================
//                                SURFACE IMPL
//==============================================================================
using Surface = WrapObstacle::Impl::Surface;

//------------------------------------------------------------------------------
//                               REALIZE CACHE
//------------------------------------------------------------------------------
void Surface::realizeTopology(State &s)
{
	// Allocate an auto-update discrete variable for the last computed geodesic.
	CacheEntry cache {};
	m_CacheIx = m_Subsystem.allocateAutoUpdateDiscreteVariable(s, Stage::Velocity, new Value<CacheEntry>(cache), Stage::Position);
}

void Surface::realizePosition(const State &s) const
{
	if (!m_Subsystem.isDiscreteVarUpdateValueRealized(s, m_CacheIx)) {
		updCacheEntry(s) = getPrevCacheEntry(s);
		m_Subsystem.markDiscreteVarUpdateValueRealized(s, m_CacheIx);
	}
}

//------------------------------------------------------------------------------
//                               PUBLIC METHODS
//------------------------------------------------------------------------------
const LocalGeodesicInfo& Surface::calcInitialGeodesic(State& s, const GeodesicInitialConditions& g0) const
{
    // TODO is this correct?
    CacheEntry& cache = updPrevCacheEntry(s);
    calcGeodesic(g0, cache);
    updCacheEntry(s) = cache;
    m_Subsystem.markDiscreteVarUpdateValueRealized(s, m_CacheIx);
}

const LocalGeodesicInfo& Surface::calcLocalGeodesic(const State& s, Vec3 prev_QS, Vec3 next_PS) const
{
    realizePosition(s);
    CacheEntry& cache = updCacheEntry(s);
    calcStatus(prev_QS, next_PS, cache);
    return cache;
}

void Surface::applyGeodesicCorrection(const State& s, const Correction& c) const
{
    realizePosition(s);

	// Get the previous geodesic.
	const CacheEntry& g = getCacheEntry(s);

    // Get corrected initial conditions.
    const GeodesicInitialConditions g0 = GeodesicInitialConditions::CreateCorrected(g.KP, g.dKP, g.length, c);

	// Shoot the new geodesic.
	calcGeodesic(g0, updCacheEntry(s));
}

void Surface::calcGeodesic(const GeodesicInitialConditions& g0, CacheEntry& cache) const
{
	// Compute geodesic start boundary frame and variation.
	m_Geometry.calcNearestFrenetFrameFast(g0.x, g0.t, cache.KP);
	m_Geometry.calcGeodesicStartFrameVariation(cache.KP, cache.dKP);

    // Compute geodesic end boundary frame amd variation (shoot new geodesic).
	m_Geometry.calcGeodesicEndFrameVariationImplicitly(
			cache.KP.p(),
			cache.KP.R().getAxisUnitVec(ContactGeometry::TangentAxis),
			g0.l,
			cache.sHint,
			cache.KQ,
			cache.dKQ,
			cache.points);

	// TODO  update step size.
	// TODO  update line tracking?
	throw std::runtime_error("NOTYETIMPLEMENTED");
}

void Surface::calcStatus(const Vec3& prev_QS, const Vec3& next_PS, CacheEntry& cache) const
{
    LocalGeodesicInfo& g = cache;

    if (g.status == Status::Disabled) {return;}

    // Make sure that the previous point does not lie inside the surface.
    if (m_Geometry.calcSurfaceValue(prev_QS)) {
        // TODO use proper assert.
        throw std::runtime_error("Unable to wrap over surface: Preceding point lies inside the surface");
    }
    if (m_Geometry.calcSurfaceValue(next_PS)) {
        // TODO use proper assert.
        throw std::runtime_error("Unable to wrap over surface: Next point lies inside the surface");
    }

    Vec3& pTrack = cache.trackingPointOnLine;

    // Detect touchdown.
    bool detectedTouchdown = g.status == Status::Liftoff;
    detectedTouchdown &= m_Geometry.calcNearestPointOnLine(prev_QS, next_PS, pTrack, m_TouchdownIter, m_TouchdownAccuracy);
    if (detectedTouchdown) {
        g.status = Status::Ok;
        GeodesicInitialConditions g0 = GeodesicInitialConditions::CreateAtTouchdown(prev_QS, next_PS, cache.trackingPointOnLine);
        calcGeodesic(g0, cache);
    }

    // Detect liftoff.
    bool detectedLiftoff = g.status == Status::Ok;
    detectedLiftoff &= g.length == 0.;
    detectedLiftoff &= dot(prev_QS - g.KP.p(), g.KP.R().getAxisUnitVec(NormalAxis)) > 0.;
    detectedLiftoff &= dot(next_PS - g.KP.p(), g.KP.R().getAxisUnitVec(NormalAxis)) > 0.;
    if (detectedLiftoff) {
        g.status = Status::Liftoff;
        pTrack = g.KP.p();
    }
}

size_t Surface::calcPathPoints(const State& s, std::vector<Vec3>& points) const
{
    realizePosition(s);

    size_t count = 0;
    const CacheEntry& cache = getCacheEntry(s);
    for (Vec3 p: cache.points) {
        points.push_back(p);
        ++count;
    }
	throw std::runtime_error("NOTYETIMPLEMENTED for analytic");
    return count;
}

//------------------------------------------------------------------------------
//                               PRIVATE METHODS
//------------------------------------------------------------------------------

//==============================================================================
//                               OBSTACLE
//==============================================================================

bool WrapObstacle::Impl::isActive(const State& s) const
{
    return getPosInfo(s).status == Status::Ok;
}

void WrapObstacle::Impl::realizeTopology(State &s)
{
	// Allocate position level cache.
	PosInfo posInfo {};
	m_PosInfoIx = m_Subsystem.allocateCacheEntry(s, Stage::Position, new Value<PosInfo>(posInfo));
}

void WrapObstacle::Impl::realizePosition(const State &s) const
{
	if (!m_Subsystem.isCacheValueRealized(s, m_PosInfoIx)) {
		calcPosInfo(s, updPosInfo(s));
		m_Subsystem.markCacheValueRealized(s, m_PosInfoIx);
	}
}

const WrapObstacle::Impl::PosInfo& WrapObstacle::Impl::getPosInfo(const State &s) const
{
    realizePosition(s);
    return Value<PosInfo>::downcast(m_Subsystem.getCacheEntry(s, m_PosInfoIx));
}

void WrapObstacle::Impl::applyGeodesicCorrection(const State& s, const WrapObstacle::Impl::Correction& c) const
{
    // Apply correction to curve.
    m_Surface.getImpl().applyGeodesicCorrection(s, c);

	// Invalidate position level cache.
	m_Subsystem.markCacheValueNotRealized(s, m_PosInfoIx);
}

void WrapObstacle::Impl::calcPosInfo(
		const State& s,
		PosInfo& posInfo) const
{
    if(isDisabled(s))
    {
        return;
    }

	// Get tramsform from local surface frame to ground.
	Transform X_GS = m_Mobod.getBodyTransform(s).compose(m_Offset);

    // Get the path points before and after this segment.
    Vec3 prev_G = m_Path.getImpl().findPrevPoint(s, m_PathSegmentIx);
    Vec3 next_G = m_Path.getImpl().findNextPoint(s, m_PathSegmentIx);

    // Transform the prev and next path points to the surface frame.
    Vec3 prev_S = X_GS.shiftBaseStationToFrame(prev_G);
    Vec3 next_S = X_GS.shiftBaseStationToFrame(next_G);

    // Detect liftoff, touchdown and potential invalid configurations.
    // TODO this doesnt follow the regular invalidation scheme...
    m_Surface.getImpl().calcUpdatedStatus(s, prev_S, next_S);

	// Grab the last geodesic that was computed.
	const Surface::LocalGeodesic& geodesic_S = m_Surface.getImpl().getGeodesic(s);

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

void WrapObstacle::Impl::calcGeodesic(State& s, Vec3 x, Vec3 t, Real l) const
{
    m_Surface.getImpl().calcGeodesic(s, x, t, l);
	m_Subsystem.markCacheValueNotRealized(s, m_PosInfoIx);
}

//==============================================================================
//                                PATH HELPERS
//==============================================================================

namespace
{
    using GeodesicJacobian = Vec4;
    using PointVariation = ContactGeometry::GeodesicPointVariation;
    using FrameVariation = ContactGeometry::GeodesicFrameVariation;
    using Variation = ContactGeometry::GeodesicVariation;
    using Correction = ContactGeometry::GeodesicCorrection;
/* using LocalGeodesic = WrapObstacle::LocalGeodesicInfo; */

const WrapObstacle* FindPrevActiveObstacle(
	const State& s,
    const std::vector<WrapObstacle>& obs,
    size_t idx)
{
    for (ptrdiff_t i = idx - 1; i > 0; --i) {
        // Find the active segment before the current.
        if (obs.at(i).isActive(s)) {
            return &obs.at(i);
        }
    }
    return nullptr;
}

const WrapObstacle* FindNextActiveObstacle(
	const State& s,
    const std::vector<WrapObstacle>& obs,
    size_t idx)
{
    // Find the active segment after the current.
    for (ptrdiff_t i = idx + 1; i < obs.size(); ++i) {
        if (obs.at(i).isActive(s)) {
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

void WrappingPath::Impl::realizeTopology(State &s)
{
	PosInfo posInfo {};
	m_PosInfoIx = m_Subsystem.allocateCacheEntry(s, Stage::Position, new Value<PosInfo>(posInfo));
}

void WrappingPath::Impl::realizePosition(const State &s) const
{
	if (m_Subsystem.isCacheValueRealized(s, m_PosInfoIx)) {return;}
	calcPosInfo(s, updPosInfo(s));
	m_Subsystem.markCacheValueRealized(s, m_PosInfoIx);
}

const WrappingPath::Impl::PosInfo& WrappingPath::Impl::getPosInfo(const State &s) const
{
	realizePosition(s);
    return Value<PosInfo>::downcast(m_Subsystem.getCacheEntry(s, m_PosInfoIx));
}

WrappingPath::Impl::PosInfo& WrappingPath::Impl::updPosInfo(const State &s) const
{
    return Value<PosInfo>::updDowncast(m_Subsystem.updCacheEntry(s, m_PosInfoIx));
}

void WrappingPath::Impl::calcInitZeroLengthGeodesic(State& s, std::function<Vec3(size_t)> GetInitPointGuess) const
{
	Vec3 prev = m_OriginBody.getBodyTransform(s).shiftFrameStationToBase(m_OriginPoint);
	size_t i = 0;
	for (const WrapObstacle& obstacle: m_Obstacles) {
	    Vec3 xGuess = GetInitPointGuess(++i);
	    obstacle.getImpl().calcGeodesic(s, xGuess, xGuess - prev, 0.);

		prev = obstacle.getImpl().getPosInfo(s).KQ.p();
	}
}

template<size_t N>
void WrappingPath::Impl::calcPathErrorVector(
	const State& s,
    const std::vector<WrapObstacle>& obs,
    const std::vector<LineSegment>& lines,
    std::array<CoordinateAxis, N> axes,
    Vector& pathError)
{
    size_t lineIx = 0;
    ptrdiff_t row  = -1;

    for (const WrapObstacle& o : obs) {
        if (!o.isActive(s)) {
            continue;
        }

        const GeodesicInfo& g = o.getImpl().getPosInfo(s);
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
	const State& s,
    const std::vector<WrapObstacle>& obs,
    const std::vector<LineSegment>& lines,
    std::array<CoordinateAxis, N> axes,
    Matrix& J)
{
    constexpr size_t Nq = GeodesicDOF;
    const size_t n     = countActive(s, obs);

    SimTK_ASSERT(J.rows() == n * N, "Invalid number of rows in jacobian matrix");
    SimTK_ASSERT(J.cols() == n * Nq, "Invalid number of columns in jacobian matrix");

    size_t row = 0;
    size_t col = 0;
    for (size_t i = 0; i < n; ++i) {
        if (!obs.at(i).isActive(s)) {
            continue;
        }
        const GeodesicInfo& g = obs.at(i).getImpl().getPosInfo(s);

        const LineSegment& l_P = lines.at(i);
        const LineSegment& l_Q = lines.at(i + 1);

        const WrapObstacle* prev = FindPrevActiveObstacle(s, obs, i);
        const WrapObstacle* next = FindNextActiveObstacle(s, obs, i);

        for (CoordinateAxis axis: axes) {
            const UnitVec3 a_P = g.KP.R().getAxisUnitVec(axis);
            const Variation& dK_P = g.dKP;

			addPathErrorJacobian(l_P, a_P, dK_P, J.block(row, col, 1, Nq));

			if (prev) {
				const Variation& prev_dK_Q = prev->getImpl().getPosInfo(s).dKQ;
				addDirectionJacobian(l_P, a_P, prev_dK_Q[1], J.block(row, col-Nq, 1, Nq), true);
			}
			++row;
		}

        for (CoordinateAxis axis: axes) {
            const UnitVec3 a_Q = g.KQ.R().getAxisUnitVec(axis);
            const Variation& dK_Q = g.dKQ;

			addPathErrorJacobian(l_Q, a_Q, dK_Q, J.block(row, col, 1, Nq), true);

			if (next) {
				const Variation& next_dK_P = next->getImpl().getPosInfo(s).dKP;
				addDirectionJacobian(l_Q, a_Q, next_dK_P[1], J.block(row, col+Nq, 1, Nq));
			}
			++row;
		}

        col += Nq;
    };
}

double WrappingPath::Impl::calcPathLength(
	const State& s,
    const std::vector<WrapObstacle>& obs,
    const std::vector<LineSegment>& lines)
{
    double lTot = 0.;
    for (const LineSegment& line : lines) {
        // TODO spell out as length.
        lTot += line.l;
    }

    for (const WrapObstacle& obstacle : obs) {
        if (!obstacle.isActive(s))
		{
            continue;
		}
        lTot += obstacle.getImpl().getPosInfo(s).length;
    }
    return lTot;
}

void WrappingPath::Impl::calcLineSegments(
	const State& s,
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
        if (!o.isActive(s)) {
            continue;
        }

        const GeodesicInfo& g = o.getImpl().getPosInfo(s);
        const Vec3 lineEnd = g.KP.p();
        lines.emplace_back(lineStart, lineEnd);

        lineStart = g.KQ.p();
    }
    lines.emplace_back(lineStart, p_I);
}

Vec3 WrappingPath::Impl::FindPrevPoint(
        const State& s,
        const Vec3& originPoint,
    const std::vector<WrapObstacle>& obs,
    size_t idx)
{
    const WrapObstacle* prev = FindPrevActiveObstacle(s, obs, idx);
    return prev ? prev->getImpl().getPosInfo(s).KQ.p() : originPoint;
}

Vec3 WrappingPath::Impl::FindNextPoint(
        const State& s,
        const Vec3& terminationPoint,
    const std::vector<WrapObstacle>& obs,
    size_t idx)
{
    const WrapObstacle* next = FindNextActiveObstacle(s, obs, idx);
    return next ? next->getImpl().getPosInfo(s).KP.p() : terminationPoint;
}

size_t WrappingPath::Impl::countActive(const State& s) const
{
    size_t count = 0;
    for (const WrapObstacle& o : m_Obstacles) {
        if (o.getImpl().isActive(s)) {
            ++count;
        }
    }
    return count;
}

void WrappingPath::Impl::calcPosInfo(const State& s, PosInfo& posInfo) const
{
	// Path origin and termination point.
	const Vec3 x_O = m_OriginBody.getBodyTransform(s).shiftFrameStationToBase(m_OriginPoint);
	const Vec3 x_I = m_TerminationBody.getBodyTransform(s).shiftFrameStationToBase(m_TerminationPoint);

	const std::array<CoordinateAxis, 2> axes {NormalAxis, BinormalAxis};

    for (posInfo.loopIter = 0; posInfo.loopIter < m_PathMaxIter; ++posInfo.loopIter) {
        const size_t nActive = countActive(s);

        // Grab the shared data cache for computing the matrices, and lock it.
        SolverDataCache& data = findDataCache(nActive);
        data.lock();

        // Compute the straight-line segments.
        calcLineSegments(s, x_O, x_I, m_Obstacles, data.updData().lineSegments);

        // Evaluate path error, and stop when converged.
        calcPathErrorVector<2>(s, m_Obstacles, posInfo.lines, axes, data.updData().pathError);
        const Real maxPathError = posInfo.pathError.normInf();
        if (maxPathError < m_PathErrorBound) {
            return;
        }

        // Evaluate the path error jacobian.
        calcPathErrorJacobian<2>(s, m_Obstacles, posInfo.lines, axes, posInfo.pathErrorJacobian);

        // Compute path corrections.
		const Correction* corrIt = calcPathCorrections(data.updData());

        // Apply path corrections.
        for (const WrapObstacle& obstacle : m_Obstacles) {
            if (!obstacle.isActive(s)) {
                continue;
            }
            obstacle.getImpl().applyGeodesicCorrection(s, *corrIt);
            ++corrIt;
        }

        // Release the lock on the shared data.
        data.unlock();

        // Path has changed: invalidate each segment's cache.
        for (const WrapObstacle& obstacle : m_Obstacles) {
            obstacle.getImpl().invalidatePositionLevelCache(s);
        }
    }

	throw std::runtime_error("Failed to converge");
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
