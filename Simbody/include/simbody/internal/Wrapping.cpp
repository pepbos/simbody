#include "SimTKcommon/internal/CoordinateAxis.h"
#include "SimTKmath.h"
#include "Wrapping.h"
#include "WrappingImpl.h"
#include "simmath/internal/ContactGeometry.h"
#include <stdexcept>

using namespace SimTK;

//==============================================================================
//                                SOME MATH
//==============================================================================
namespace
{

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
	// Transform the local geodesic to ground frame.
	Transform X_GS = m_Mobod.getBodyTransform(state).compose(m_Offset);
	const Surface::LocalGeodesic& geodesic_S = m_Surface.getImpl().getGeodesic(state);

    posInfo.KP = X_GS.compose(geodesic_S.KP);
    posInfo.KQ = X_GS.compose(geodesic_S.KQ);

    posInfo.dKP[0] = X_GS.R() * geodesic_S.dKP[0];
    posInfo.dKP[1] = X_GS.R() * geodesic_S.dKP[1];

    posInfo.dKQ[0] = X_GS.R() * geodesic_S.dKQ[0];
    posInfo.dKQ[1] = X_GS.R() * geodesic_S.dKQ[1];

    posInfo.length = geodesic_S.length;

    throw std::runtime_error("NOTYETIMPLEMENTED: Check transformation order");
}

void WrapObstacle::Impl::calcPosInfo(
		const State& state,
		Vec3 prev,
		Vec3 next,
		size_t maxIter, Real eps) const
{
    m_Surface.getImpl().calcWrappingStatus(state, prev, next, maxIter, eps);

    if(isDisabled(state))
    {
        return;
    }

    calcPosInfo(state, updPosInfo(state));
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
    using Geodesic = WrapObstacle::Impl::PosInfo;
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

}

//==============================================================================
//                                PATH IMPL
//==============================================================================

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

//==============================================================================
//                                WRAP OBSTACLE IMPL
//==============================================================================

namespace
{

static const int N_PATH_CONSTRAINTS = 4;

using ActiveLambda = std::function<
    void(size_t prev, size_t next, size_t current, bool isActive)>;

void CallCurrentWithPrevAndNext(
    size_t prev,
    size_t next,
    size_t current,
    bool isActive,
    ActiveLambda& f)
{
    f(prev, next, current, isActive);
}

template <typename... FUNCS>
void CallCurrentWithPrevAndNext(
    size_t prev,
    size_t next,
    size_t current,
    bool isActive,
    ActiveLambda& f,
    FUNCS&&... fs)
{
    f(prev, next, current, isActive);
    CallCurrentWithPrevAndNext(
        prev,
        next,
        current,
        isActive,
        std::forward<FUNCS>(fs)...);
}

template <typename... FUNCS>
void MapWithPrevAndNext(
	const State& s,
    const std::vector<WrapObstacle>& obs,
    ActiveLambda& f,
    FUNCS&&... fs)
{
    const ptrdiff_t n = obs.size();
    ptrdiff_t next    = 0;
    ptrdiff_t prev    = -1;

    for (ptrdiff_t i = 0; i < n; ++i) {
        // Find the active segment before the current.
        if (i > 0) {
            if (obs.at(i - 1).isActive(s)) {
                prev = i - 1;
            }
        }

        // Find the active segment after the current.
        if (next <= i) {
            for (; ++next < n;) {
                const WrapObstacle& o = obs.at(next);
                if (o.isActive(s)) {
                    break;
                }
            }
        }

        CallCurrentWithPrevAndNext(
            prev < 0 ? n : prev,
            next,
            i,
            obs.at(i).isActive(s),
            f,
            std::forward<FUNCS>(fs)...);
    }
}

void addDirectionJacobian(
    const LineSegment& e,
	const UnitVec3& axis,
    const PointVariation& dx,
    MatrixView J,
    bool invert = false)
{
    Vec3 y = axis - e.d * dot(e.d,axis);
    y /= e.l * (invert ? 1. : -1);
    J += ~dx * y;
}

double calcPathError(const LineSegment& e, const Rotation& R, CoordinateAxis axis)
{
    return dot(e.d, R.getAxisUnitVec(axis));
}

GeodesicJacobian addPathErrorJacobian(
    const LineSegment& e,
	const UnitVec3& axis,
    const Variation& dK,
	MatrixView J,
    bool invertV = false)
{
	addDirectionJacobian(e, axis, dK[1], J, invertV);
    J += cross(axis,e.d).transpose() * dK[0];
}

void calcPathError(
	const State& state,
    const std::vector<WrapObstacle>& obs,
    const std::vector<LineSegment>& lines,
    Vector& pathError)
{
    size_t i      = 0;
    ptrdiff_t row = -1;
    for (const WrapObstacle& o : obs) {
        const Geodesic& g = o.getGeodesic(state);
        pathError(++row)  = calcPathError(lines.at(i), g.KP.R(), NormalAxis);
        pathError(++row)  = calcPathError(lines.at(i), g.KP.R(), BinormalAxis);
        ++i;
        pathError(++row)  = calcPathError(lines.at(i), g.KQ.R(), NormalAxis);
        pathError(++row)  = calcPathError(lines.at(i), g.KQ.R(), BinormalAxis);
    }
}

void writePathErrorJaobianBlock(
		const Transform& K,
		const Variation& dK,
		const LineSegment& l,
		MatrixView J)
{
	const UnitVec3 nP = K.R().getAxisUnitVec(NormalAxis);
	const UnitVec3 bP = K.R().getAxisUnitVec(BinormalAxis);
	
}

void calcPathErrorJacobian(
	const State& s,
    const std::vector<WrapObstacle>& obs,
    const std::vector<LineSegment>& lines,
    Matrix& J)
{
    constexpr size_t Q = GeodesicDOF;
    constexpr size_t C = N_PATH_CONSTRAINTS;
    const size_t n     = countActive(s, obs);

    J.resize(n * C, n * Q);

    size_t row = 0;
    size_t col = 0;
    Vec4 block { NaN, NaN, NaN, NaN,};

    ActiveLambda f = [&](size_t prevIx, size_t nextIx, size_t i, bool isActive) {
        if (!isActive) {
            return;
        }

        const Geodesic& g = obs.at(i).getGeodesic(s);

		{
			const LineSegment& l_P = lines.at(i);

			const UnitVec3 nP = g.KP.R().getAxisUnitVec(NormalAxis);
			const UnitVec3 bP = g.KP.R().getAxisUnitVec(BinormalAxis);

			addPathErrorJacobian(l_P, nP, g.dKP, J.block(row, col, 1, Q));
			addPathErrorJacobian(l_P, bP, g.dKP, J.block(row + 1, col, 1, Q));

			if (prevIx != n) {
				const Geodesic& prev = obs.at(prevIx).getGeodesic(s);

				addDirectionJacobian(l_P, nP, prev.dKP[1], J.block(row, col-Q, 1, Q), true);
				addDirectionJacobian(l_P, bP, prev.dKP[1], J.block(row + 1, col-Q, 1, Q), true);
			}
			row += 2;
		}

		// TODO looked like this
		/* calcPathErrorJacobian(l_Q, a_Q, g.v_Q, g.w_Q, true); */
		// dx_P  dR_P
		// addDirectionJacobian(l_Q, n_Q, prev.dx_P, J.block())
		{
			const LineSegment& l_Q = lines.at(i + 1);

			const UnitVec3 nQ = g.KQ.R().getAxisUnitVec(NormalAxis);
			const UnitVec3 bQ = g.KQ.R().getAxisUnitVec(BinormalAxis);

			addPathErrorJacobian(l_Q, nQ, g.dKQ, J.block(row, col, 1, Q), true);
			addPathErrorJacobian(l_Q, bQ, g.dKQ, J.block(row + 1, col, 1, Q), true);

			if (nextIx != n) {
				/* pathErrorJacobian.block<1, Q>(row, col + Q) = */
				/* 	calcDirectionJacobian(l_Q, a_Q, obs.at(next).getGeodesic().v_P); */

				const Geodesic& next = obs.at(nextIx).getGeodesic(s);

				addDirectionJacobian(l_Q, nQ, next.dKP[1], J.block(row, col+Q, 1, Q));
				addDirectionJacobian(l_Q, bQ, next.dKP[1], J.block(row + 1, col+Q, 1, Q));
			}
			++row;
		}

        col += Q;
    };

    MapWithPrevAndNext(s, obs, f);
}

void calcPathCorrections(const State& s, const std::vector<WrapObstacle>& obs, const Vector& pathError, const Matrix& pathErrorJacobian, Matrix& pathMatrix, Vector& pathCorrections);

double calcPathLength(
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
        lTot += obstacle.getGeodesic(s).length;
    }
    return lTot;
}

void calcLineSegments(
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
    for (size_t i = 0; i < n; ++i) {
        if (!obs.at(i).isActive(s)) {
            continue;
        }

        const Geodesic& g = obs.at(i).getGeodesic(s);
        const Vec3 lineEnd   = g.KP.p();
        lines.emplace_back(lineStart, lineEnd);

        lineStart = g.KQ.p();
    }
    lines.emplace_back(lineStart, p_I);
}

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


//==============================================================================
//                               WRAPPING PATH
//==============================================================================

void WrappingPath::Impl::realizeTopology(State &state)
{
	PosInfo posInfo {};
	m_PosInfoIx = m_Subsystem.allocateCacheEntry(state, Stage::Position, new Value<PosInfo>(posInfo));
}

void WrappingPath::Impl::realizePosition(const State &state) const
{
	if (m_Subsystem.isCacheValueRealized(state, m_PosInfoIx)) {return;}
	calcPosInfo(updPosInfo(state));
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

bool WrapObstacle::isPointBelowSurface(const State& state, Vec3 point) const
{
	const Transform& X_BS = getSurfaceToBaseTransform(state);
	return m_Surface.getGeometry().calcSurfaceValue(X_BS.shiftBaseStationToFrame(point)) < 0.;
}

void WrappingPath::Impl::calcPosInfo(const State& s, PosInfo& posInfo, bool preventLiftOff = false) const
{
	// Path oigin and termination points.
	Vec3 x_O = m_OriginBody.getBodyTransform(s).shiftFrameStationToBase(m_OriginPoint);
	Vec3 x_I = m_TerminationBody.getBodyTransform(s).shiftFrameStationToBase(m_TerminationPoint);

	// Helper for detecting if start/end points lie inside the surfacce.
	std::function<void(Vec3, const WrapObstacle&, Vec3)> DetectInsideSurfaceError =
		[&](Vec3 x_O, const WrapObstacle& obstacle, Vec3 x_I) { // TODO also supply the status.
			// Check that previous point does not lie inside the surface.
			if (obstacle.isPointBelowSurface(s, x_O)) {
				throw std::runtime_error("Start point lies inside the surface");
			}

			// Check that next point does not lie inside the surface.
			if (obstacle.isPointBelowSurface(s, x_I)) {
				throw std::runtime_error("End point lies inside the surface");
			}

			if (preventLiftOff) {
				return;
			}

			// TODO detect touhdown
			// TODO detect liftoff
			throw std::runtime_error("NOTYETIMPLEMENTED");
		};

    for (posInfo.loopIter = 0; posInfo.loopIter < m_MaxIter; ++posInfo.loopIter) {

        // Detect touchdown & liftoff.
		// doForEachObjectWithPrevAndNextPoint(UpdateLiftoffAndTouchdown)
		callCurrentWithPrevAndNext(s, DetectInsideSurfaceError);

        // Compute the line segments.
        calcLineSegments(s, x_O, x_I, m_Obstacles, posInfo.lines);

        calcPathError(s, m_Obstacles, posInfo.lines, posInfo.pathError);

        const Real maxPathError = posInfo.pathError.normInf();

        // Evaluate path error, and stop when converged.
        if (maxPathError < m_PathErrorBound) {
            return;
        }

        // Evaluate the path error jacobian.
        calcPathErrorJacobian(s, m_Obstacles, posInfo.lines, posInfo.pathErrorJacobian);

        // Compute path corrections.
        calcPathCorrections(s, m_Obstacles, posInfo.pathError, posInfo.pathErrorJacobian, posInfo.pathMatrix, posInfo.pathCorrections);

        // Apply path corrections.
		throw std::runtime_error("NOTYETIMPLEMENTED");
		const Correction* corrIt = nullptr; // TODO
        for (const WrapObstacle& obstacle : m_Obstacles) {
            if (!obstacle.isActive(s)) {
                continue;
            }
            obstacle.applyGeodesicCorrection(s, *corrIt);
            ++corrIt;
        }
    }

	throw std::runtime_error("Failed to converge");
}
