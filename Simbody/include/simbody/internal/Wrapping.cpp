#include "SimTKmath.h"
#include "Wrapping.h"
#include <stdexcept>

using namespace SimTK;

//==============================================================================
//            WRAPPING PATH
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
