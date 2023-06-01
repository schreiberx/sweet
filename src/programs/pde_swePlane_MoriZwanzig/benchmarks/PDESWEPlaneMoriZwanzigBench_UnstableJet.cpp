#include "PDESWEPlaneMoriZwanzigBench_UnstableJet.hpp"

double PDESWEPlaneMoriZwanzigBench_UnstableJet::u_scale = 1.;

void PDESWEPlaneMoriZwanzigBench_UnstableJet::set_u_scale()
{
	this->u_scale = 50. * std::sqrt(shackPDESWEPlane->h0 / 1e4 * shackPDESWEPlane->gravitation / 9.80616);
	std::cout << "U SCALE " << this->u_scale << std::endl;
}

