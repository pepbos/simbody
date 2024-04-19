#include "SimTKmath.h"
#include "Wrapping.h"
#include "simmath/internal/ContactGeometry.h"
#include <stdexcept>

using namespace SimTK;

//==============================================================================
//                                SOME MATH
//==============================================================================
namespace
{

using GeodesicJacobian = Vec4;
using LineSegment = WrappingPathImpl::LineSegment;
using PointVariation = ContactGeometry::GeodesicPointVariation;
using FrameVariation = ContactGeometry::GeodesicFrameVariation;
using Variation = ContactGeometry::GeodesicVariation;

using Geodesic = WrapObstacle::PosInfo;
using LocalGeodesic = WrapObstacle::LocalGeodesicInfo;

static const CoordinateAxis TangentAxis = ContactGeometry::TangentAxis;
static const CoordinateAxis NormalAxis = ContactGeometry::NormalAxis;
static const CoordinateAxis BinormalAxis = ContactGeometry::BinormalAxis;
static const int GeodesicDOF = 4;
static const int N_PATH_CONSTRAINTS = 4;

using Status = WrapObstacle::Status;

	Vec3 calcCorrectedGeodesicStartPoint(
			const ContactGeometry::GeodesicCorrection& c,
			const ContactGeometry::GeodesicVariation& dKP,
			const ContactGeometry::FrenetFrame& KP
			)
	{
		Vec3 v = dKP[1] * c;
		return KP.p() + v;
	}

	Vec3 calcCorrectedGeodesicStartTangent(
			const ContactGeometry::GeodesicCorrection& c,
			const ContactGeometry::GeodesicVariation& dKP,
			const ContactGeometry::FrenetFrame& KP
			)
	{
		Vec3 w = dKP[0] * c;
		const UnitVec3 t = KP.R().getAxisUnitVec(ContactGeometry::TangentAxis);
		return t + cross(w, t);
	}

	double calcCorrectedGeodesicLength(
			const ContactGeometry::GeodesicCorrection& c,
			Real length)
	{
		return length + c[3];
	}

	void xformSurfaceGeodesicToBase(
			const WrapObstacle::LocalGeodesicInfo& geodesic_S,
			WrapObstacle::PosInfo& geodesic_B,
			Transform& X_BS) {

		// TODO check transformation order.
		geodesic_B.KP = X_BS.compose(geodesic_S.KP);
		geodesic_B.KQ = X_BS.compose(geodesic_S.KQ);

		geodesic_B.dKP[0] = X_BS.R() * geodesic_S.dKP[0];
		geodesic_B.dKP[1] = X_BS.R() * geodesic_S.dKP[1];

		geodesic_B.dKQ[0] = X_BS.R() * geodesic_S.dKQ[0];
		geodesic_B.dKQ[1] = X_BS.R() * geodesic_S.dKQ[1];

		geodesic_B.length = geodesic_S.length;

		throw std::runtime_error("NOTYETIMPLEMENTED: Check transformation order");
	}

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

const WrapObstacle::LocalGeodesicInfo& WrapObstacle::calcInitZeroLengthGeodesicGuess(State& s, Vec3 xPrev) const
{
	Vec3 x = getInitialPointGuess();
	Vec3 t = (x - xPrev);

	// Shoot a zero-length geodesic as initial guess.
	calcLocalGeodesic(x, t, 0., 0., updPrevLocalGeodesicInfo(s));
}

void WrappingPath::Impl::calcInitZeroLengthGeodesicSegmentsGuess(State& s) const
{
	Vec3 prev = m_OriginBody.getBodyTransform(s).shiftFrameStationToBase(m_OriginPoint);
	for (const WrapObstacle& obstacle: m_Obstacles) {
		prev = obstacle.calcInitZeroLengthGeodesicGuess(s, prev).KQ.p();
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

void WrappingPath::Impl::calcPosInfo(PosInfo& posInfo) const
{
	throw std::runtime_error("NOTYETIMPLEMENTED");
}

//==============================================================================
//                               OBSTACLE
//==============================================================================

void WrapObstacle::realizeTopology(State &state)
{
	// Allocate an auto-update discrete variable for the last computed geodesic.
	LocalGeodesicInfo warmStartInfo {};
	m_WarmStartInfoDIx = m_Subsystem.allocateAutoUpdateDiscreteVariable(state, Stage::Velocity, new Value<LocalGeodesicInfo>(warmStartInfo), Stage::Position);

	// Allocate position level cache.
	PosInfo posInfo {};
	m_PosInfoIx = m_Subsystem.allocateCacheEntry(state, Stage::Position, new Value<PosInfo>(posInfo));
}

void WrapObstacle::realizePosition(const State &state) const
{
	// Set the current local geodesic equal to the previous.
	if (!m_Subsystem.isDiscreteVarUpdateValueRealized(state, m_WarmStartInfoDIx)) {
		updLocalGeodesicInfo(state) = getPrevLocalGeodesicInfo(state);
		m_Subsystem.markDiscreteVarUpdateValueRealized(state, m_WarmStartInfoDIx);
	}
}

const WrapObstacle::PosInfo& WrapObstacle::getGeodesic(const State &state) const
{
	if (!m_Subsystem.isCacheValueRealized(state, m_PosInfoIx)) {
		calcPosInfo(state, getLocalGeodesicInfo(state), updPosInfo(state));
		m_Subsystem.markCacheValueRealized(state, m_PosInfoIx);
	}
    return Value<PosInfo>::downcast(m_Subsystem.getCacheEntry(state, m_PosInfoIx));
}

void WrapObstacle::applyGeodesicCorrection(const State& state, const WrapObstacle::Correction& c) const
{
	// Get prev geodesic.
	const LocalGeodesicInfo& geodesic = getLocalGeodesicInfo(state);

	// Apply geodesic correction.
	Vec3 x = calcCorrectedGeodesicStartPoint(c, geodesic.dKP, geodesic.KP);
	Vec3 t = calcCorrectedGeodesicStartTangent(c, geodesic.dKP, geodesic.KP);
	Real l = calcCorrectedGeodesicLength(c, geodesic.length);
	Real sHint = geodesic.sHint;

	// Shoot the new geodesic.
	calcLocalGeodesic(x, t, l, sHint, updLocalGeodesicInfo(state));

	// Invalidate position level cache.
	m_Subsystem.markCacheValueNotRealized(state, m_PosInfoIx);
}

void WrapObstacle::calcLocalGeodesic( Vec3 x, Vec3 t, Real l, Real sHint,
		WrapObstacle::LocalGeodesicInfo& geodesic) const
{
	const ContactGeometry& geometry = m_Surface.getGeometry();

	// Compute geodesic start boundary frame.
	geometry.calcNearestFrenetFrameFast(x, t, geodesic.KP);
	geometry.calcGeodesicStartFrameVariation(geodesic.KP, geodesic.dKP);

	// Compute geodesic end boundary frame (shoot new geodesic).
	geometry.calcGeodesicEndFrameVariationImplicitly(
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

void WrapObstacle::calcPosInfo(
		const State& state,
		const WrapObstacle::LocalGeodesicInfo& localGeodesic,
		PosInfo& posInfo) const
{
	// Transform the local geodesic to ground frame.
	Transform X_BS = m_Surface.calcSurfaceToGroundTransform(state);
	xformSurfaceGeodesicToBase(localGeodesic, posInfo, X_BS);
}
