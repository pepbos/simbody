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

GeodesicJacobian calcDirectionJacobian(
    const LineSegment& e,
    UnitVec3 axis,
    const PointVariation& v)
{
    Vec3 y = axis - e.d * dot(e.d,axis);
    y /= e.l;
    return ~v * y;
}

GeodesicJacobian calcPathErrorJacobian(
    const LineSegment& line,
    UnitVec3 axis,
    const PointVariation& v,
    const FrameVariation& w,
    bool invertV = false)
{
    GeodesicJacobian jacobian =
        calcDirectionJacobian(line, axis, v) * (invertV ? -1. : 1.);
    jacobian += cross(axis,line.d).transpose() * w;
    return jacobian;
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

const WrapObstacle::PosInfo& WrapObstacle::getPosInfo(const State &state) const
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
