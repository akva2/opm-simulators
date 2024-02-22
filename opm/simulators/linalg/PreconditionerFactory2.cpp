#include <config.h>

#include <opm/simulators/linalg/PreconditionerFactory_impl.hpp>

namespace Opm {

INSTANCE_PF(double,2)

#if FLOW_INSTANCE_FLOAT
INSTANCE_PF(float,2)
#endif

}
