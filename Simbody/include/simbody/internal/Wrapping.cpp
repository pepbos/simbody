#include "SimTKmath.h"
#include "Wrapping.h"
#include <stdexcept>

using namespace SimTK;

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
