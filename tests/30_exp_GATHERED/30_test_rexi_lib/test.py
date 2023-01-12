#! /usr/bin/env python3

import mule_local.rexi.tests.test_cirexi_coeffs

import mule_local.rexi.tests.test_cirexi_coeffs_rescale_beta

import mule_local.rexi.tests.test_function_identities

import mule_local.rexi.tests.test_rexi_approximations


exec_program('mule.benchmark.cleanup_all', catch_output=False)
