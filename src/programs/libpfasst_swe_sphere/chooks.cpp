#include <mpi.h>
#include "SphereDataVars.hpp"
#include "SphereDataCtx.hpp"
#include "ceval.hpp"

extern "C"
{

    void cecho_error(SphereData_Spectral* sd,
                     int step)
    {
      // not implemented
    }

    void cecho_residual(SphereDataCtx *i_ctx,
                        double i_norm,
                        int i_current_proc)
    {
        // get the residual vector
        std::vector<std::vector<double> >& residuals = i_ctx->get_residuals();

        // save the residual
        residuals[i_current_proc].push_back(i_norm);
    }

    void cecho_output_invariants(SphereDataCtx *i_ctx,
                                 SphereDataVars *i_Y,
                                 int i_current_proc,
                                 int i_current_step,
                                 int i_current_iter,
                                 int i_nnodes,
                                 int i_niters
                                 )
    {
        const SphereData_Spectral& phi_Y  = i_Y->get_phi();
        const SphereData_Spectral& vort_Y = i_Y->get_vort();
        const SphereData_Spectral& div_Y  = i_Y->get_div();

        // get the current space-time level
        const int level = i_Y->get_level();

        // get the simulation variables
        SimulationVariables* simVars         = i_ctx->get_simulation_variables();

        // get the SphereDiagnostics object from context
        SphereHelpers_Diagnostics* sphereDiagnostics = i_ctx->get_sphere_diagnostics();

        // get the SphereOperators object from context
        SphereOperators_SphereData* sphereOperators     = i_ctx->get_sphere_operators(level);

        // compute the invariants
        sphereDiagnostics->update_phi_vrt_div_2_mass_energy_enstrophy(
                                       *sphereOperators,
                                       phi_Y,
                                       vort_Y,
                                       div_Y,
                                       *simVars
                                       );

        std::cout << std::setprecision(20)
              << "mass = " << simVars->diag.total_mass
              << " energy = " << simVars->diag.total_energy
              << " potential_enstrophy = " << simVars->diag.total_potential_enstrophy
              << std::endl;

        // save the invariants for plotting at the end
        i_ctx->save_physical_invariants(i_current_step);
    }

    void cecho_output_jump(SphereDataCtx *i_ctx,
                           SphereDataVars *i_Y,
                           int i_current_proc,
                           int i_current_step,
                           int i_current_iter,
                           int i_nnodes,
                           int i_niters
                           )
    {
        const SphereData_Spectral& phi_Y  = i_Y->get_phi();
        const SphereData_Spectral& vort_Y = i_Y->get_vort();
        const SphereData_Spectral& div_Y  = i_Y->get_div();

        // get the pointer to the Simulation Variables object
        SimulationVariables* simVars = i_ctx->get_simulation_variables();
        simVars->timecontrol.current_timestep_nr = i_current_step + 1;
        auto current_dt = simVars->timecontrol.current_timestep_size;
        simVars->timecontrol.current_simulation_time = (i_current_step + 1) * current_dt;

        // write the data to file
        std::string filename = "prog_jump_phi_current_proc_"+std::to_string(i_current_proc)
                                    +"_current_iter_"+std::to_string(i_current_iter)
                                    +"_nnodes_"      +std::to_string(i_nnodes)
                                    +"_niters_"      +std::to_string(i_niters);
        write_file(*i_ctx, phi_Y, filename.c_str());

        filename = "prog_jump_vort_current_proc_"+std::to_string(i_current_proc)
                        +"_current_iter_"+std::to_string(i_current_iter)
                        +"_nnodes_"      +std::to_string(i_nnodes)
                        +"_niters_"      +std::to_string(i_niters);
        write_file(*i_ctx, vort_Y, filename.c_str());

        filename = "prog_jump_div_current_proc_"+std::to_string(i_current_proc)
                        +"_current_iter_"+std::to_string(i_current_iter)
                        +"_nnodes_"      +std::to_string(i_nnodes)
                        +"_niters_"      +std::to_string(i_niters);
        write_file(*i_ctx, div_Y, filename.c_str());
    }



    void cecho_output_solution(SphereDataCtx *i_ctx,
                               SphereDataVars *i_Y,
                               int i_current_proc,
                               int i_current_step,
                               int i_current_iter,
                               int i_nnodes,
                               int i_niters
                               )
    {
        // get the pointer to the Simulation Variables object
        SimulationVariables* simVars = i_ctx->get_simulation_variables();

        // update timecontrol information
        simVars->timecontrol.current_timestep_nr = i_current_step + 1;
        auto current_dt = simVars->timecontrol.current_timestep_size;
        simVars->timecontrol.current_simulation_time = (i_current_step + 1) * current_dt;

        // check if output is necessary
        const auto output_dt = simVars->iodata.output_each_sim_seconds;
        if (output_dt < 0) {
            // no output required
            return;
        }
        if (output_dt != 0) {
            // == 0 means output at every step
            if (std::fmod(simVars->timecontrol.current_simulation_time, output_dt) != 0) {
                // we haven't reached the next output time yet -> no output required
                return;
            }
        }

        const SphereData_Spectral& phi_Y  = i_Y->get_phi();
        const SphereData_Spectral& vort_Y = i_Y->get_vort();
        const SphereData_Spectral& div_Y  = i_Y->get_div();

        // write the data to file
        std::string filename = "prog_phi_current_proc_"+std::to_string(i_current_proc)
                                  +"_current_iter_"+std::to_string(i_current_iter)
                                  +"_nnodes_"      +std::to_string(i_nnodes)
                                  +"_niters_"      +std::to_string(i_niters);
        write_file(*i_ctx, phi_Y, filename.c_str());

        filename = "prog_vort_current_proc_"+std::to_string(i_current_proc)
                      +"_current_iter_"+std::to_string(i_current_iter)
                      +"_nnodes_"      +std::to_string(i_nnodes)
                      +"_niters_"      +std::to_string(i_niters);
        write_file(*i_ctx, vort_Y, filename.c_str());

        filename = "prog_div_current_proc_"+std::to_string(i_current_proc)
                      +"_current_iter_"+std::to_string(i_current_iter)
                      +"_nnodes_"      +std::to_string(i_nnodes)
                      +"_niters_"      +std::to_string(i_niters);
        write_file(*i_ctx, div_Y, filename.c_str());

    }

}
