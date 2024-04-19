#include "SimTKmath.h"
#include "Wrapping.h"

using namespace SimTK;

//==============================================================================
//            WRAPPING PATH
//==============================================================================

void WrappingPath::Impl::realizeTopology(State &state)
{
	PosInfo posInfo {};
	m_PosInfoIx = m_Subsystem->allocateCacheEntry(state, Stage::Position, new Value<PosInfo>(posInfo));

	VelInfo velInfo {};
	m_VelInfoIx = m_Subsystem->allocateCacheEntry(state, Stage::Velocity, new Value<VelInfo>(velInfo));

	VizInfo vizInfo {};
	m_VizInfoIx = m_Subsystem->allocateCacheEntry(state, Stage::Position, new Value<VizInfo>(vizInfo));
}

void WrappingPath::Impl::realizePosition(const State &state) const
{
	if (m_Subsystem->isCacheValueRealized(state, m_PosInfoIx)) {return;}
	calcPosInfo(updPosInfo(state));
	m_Subsystem->markCacheValueRealized(state, m_PosInfoIx);
}

const WrappingPath::Impl::PosInfo& WrappingPath::Impl::getPosInfo(const State &state) const
{
	realizePosition(state);
    return Value<PosInfo>::downcast(m_Subsystem->getCacheEntry(state, m_PosInfoIx));
}

WrappingPath::Impl::PosInfo& WrappingPath::Impl::updPosInfo(const State &state) const
{
    return Value<PosInfo>::updDowncast(m_Subsystem->updCacheEntry(state, m_PosInfoIx));
}
