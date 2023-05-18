#include <mpi.h>
#include "../interface/SphereDataVars.hpp"
#include "SphereDataCtxSDC.hpp"
#include "ceval.hpp"

extern "C"
{
    void cecho_error(sweet::SphereData_Spectral *sd,
                     int step)
    {
        // not implemented
    }

    void cecho_residual(SphereDataCtxSDC *i_ctx,
                        double i_norm,
                        int i_current_proc)
    {
        // get the residual vector
        std::vector<std::vector<double>> &residuals = i_ctx->get_residuals();

        // save the residual
        residuals[i_current_proc].push_back(i_norm);
    }

    void cecho_output_invariants(SphereDataCtxSDC *i_ctx,
                                 SphereDataVars *i_Y,
                                 int i_current_proc,
                                 int i_current_step,
                                 int i_current_iter,
                                 int i_nnodes,
                                 int i_niters)
    {
        const sweet::SphereData_Spectral &phi_pert_Y = i_Y->get_phi_pert();
        const sweet::SphereData_Spectral &vrt_Y = i_Y->get_vrt();
        const sweet::SphereData_Spectral &div_Y = i_Y->get_div();

        // get the sweet::SphereOperators object from context
        sweet::SphereOperators *sphereOperators = i_ctx->get_sphere_operators();

        // compute the invariants
        i_ctx->diagnostics.update_phi_vrt_div_2_mass_energy_enstrophy(
            sphereOperators,
            phi_pert_Y,
            vrt_Y,
            div_Y,
            i_ctx->shackSphereDataOps->sphere_radius,
            i_ctx->shackPDESWESphere->gravitation);

        std::cout << "[MULE] libpfasst.mass_s" << std::setfill('0') << std::setw(5) << i_current_step;
        std::cout << " = " << std::setprecision(20) << i_ctx->diagnostics.total_mass << std::endl;
        std::cout << "[MULE] libpfasst.energy_s" << std::setfill('0') << std::setw(5) << i_current_step;
        std::cout << " = " << std::setprecision(20) << i_ctx->diagnostics.total_energy << std::endl;
        std::cout << "[MULE] libpfasst.potential_enstrophy_s" << std::setfill('0') << std::setw(5) << i_current_step;
        std::cout << " = " << std::setprecision(20) << i_ctx->diagnostics.total_potential_enstrophy << std::endl;

        // save the invariants for plotting at the end
        i_ctx->save_physical_invariants(i_current_step);
    }

    void cecho_output_solution(SphereDataCtxSDC *i_ctx,
                               SphereDataVars *i_Y,
                               int i_current_proc,
                               int i_current_step,
                               int i_current_iter,
                               int i_nnodes,
                               int i_niters)
    {
        // get the pointer to the Simulation Variables object
        // sweet::ShackDictionary *shackDict = i_ctx->get_simulation_variables();

        // update timecontrol information
        i_ctx->shackTimestepControl->current_timestep_nr = i_current_step + 1;
        auto current_dt = i_ctx->shackTimestepControl->current_timestep_size;
        i_ctx->shackTimestepControl->current_simulation_time = (i_current_step + 1) * current_dt;

        // check if we should write output
        if (!timestep_check_output(i_ctx, i_current_iter, i_niters))
        {
            return;
        }

        // update when to write output the next time
        i_ctx->shackIOData->output_next_sim_seconds += i_ctx->shackIOData->output_each_sim_seconds;

        const sweet::SphereData_Spectral &phi_pert_Y = i_Y->get_phi_pert();
        const sweet::SphereData_Spectral &vrt_Y = i_Y->get_vrt();
        const sweet::SphereData_Spectral &div_Y = i_Y->get_div();

        // write the data to file
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank != 0)
        {
            return;
        }
        std::string filename = "prog_phi_pert";
        write_file(i_ctx, phi_pert_Y, filename.c_str());

        filename = "prog_vrt";
        write_file(i_ctx, vrt_Y, filename.c_str());

        filename = "prog_div";
        write_file(i_ctx, div_Y, filename.c_str());
    }
}
