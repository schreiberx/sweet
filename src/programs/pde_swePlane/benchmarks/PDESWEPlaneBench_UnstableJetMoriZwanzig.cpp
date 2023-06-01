#include "PDESWEPlaneBench_UnstableJetMoriZwanzig.hpp"

double PDESWEPlaneBench_UnstableJetMoriZwanzig::u_scale = 1.;

void PDESWEPlaneBench_UnstableJetMoriZwanzig::set_u_scale()
{
	////this->u_scale = std::sqrt(shackPDESWEPlane->h0 / 1e4);
	this->u_scale = 50. * std::sqrt(shackPDESWEPlane->h0 / 1e4 * shackPDESWEPlane->gravitation / 9.80616);
	std::cout << "U SCALE " << this->u_scale << std::endl;
}

