#include "PlaneDataVars.hpp"
#include "PlaneDataCtx.hpp"
#include "ceval.hpp"

extern "C"
{
  void cecho_error(PlaneData* sd, 
		   int step)
  {
    // not implemented
  }

  void cecho_residual(PlaneDataCtx *i_ctx,
		      double i_norm,
		      int i_current_proc,
		      int i_current_step,
		      int i_current_iter,
		      int i_nnodes,
		      int i_niters)
  {
    // get the residual vector
    std::vector<std::vector<double> >& residuals = i_ctx->get_residuals();

    // save the residual
    residuals[i_current_proc].push_back(i_norm);
  }

  void cecho_output_invariants(
			       PlaneDataCtx *i_ctx,
			       PlaneDataVars *i_Y,
			       int i_current_proc,
			       int i_current_step,
			       int i_current_iter,
			       int i_nnodes,
			       int i_niters
			       )
  {
    const PlaneData& h_Y    = i_Y->get_h();
    const PlaneData& vort_Y = i_Y->get_vort();
    const PlaneData& div_Y  = i_Y->get_div();

    // get the current space-time level
    const int level = i_Y->get_level();

    // get the simulation variables
    SimulationVariables* simVars        = i_ctx->get_simulation_variables();

    // get the PlaneDiagnostics object from context
    PlaneDiagnostics* planeDiagnostics = i_ctx->get_plane_diagnostics();

    // get the PlaneOperators object from context
    PlaneOperators* planeOperators     = i_ctx->get_plane_operators(level);

    // compute the invariants
    // planeDiagnostics->update_h_vort_div_2_mass_energy_enstrophy(
    // 								*planeOperators,
    // 								h_Y,
    // 								vort_Y,
    // 								div_Y,
    // 								*simVars
    // 								);

    std::cout << std::setprecision(20) 
	      << "mass = " << simVars->diag.total_mass
	      << " energy = " << simVars->diag.total_energy
	      << " potential_enstrophy = " << simVars->diag.total_potential_enstrophy
	      << std::endl;

    // save the invariants for plotting at the end
    i_ctx->save_physical_invariants(i_current_step);

  }

  void cecho_output_jump(
			 PlaneDataCtx *i_ctx,
			 PlaneDataVars *i_Y,
			 int i_current_proc,
			 int i_current_step,
			 int i_current_iter,
			 int i_nnodes,
			 int i_niters
			 )
  {
    const PlaneData& h_Y    = i_Y->get_h();
    const PlaneData& vort_Y = i_Y->get_vort();
    const PlaneData& div_Y  = i_Y->get_div();

    // get the SimulationVariables object from context
    SimulationVariables* simVars(i_ctx->get_simulation_variables());

    // write the data to file
    std::string filename = "prog_jump_h_current_proc_"+std::to_string(i_current_proc)
                                +"_current_step_"+std::to_string(i_current_step)
                                +"_current_iter_"+std::to_string(i_current_iter)
                                +"_nnodes_"      +std::to_string(i_nnodes)
                                +"_niters_"      +std::to_string(i_niters);
    write_file(*i_ctx, h_Y, filename.c_str());

    filename = "prog_jump_vort_current_proc_"+std::to_string(i_current_proc)
                    +"_current_step_"+std::to_string(i_current_step)
                    +"_current_iter_"+std::to_string(i_current_iter)
                    +"_nnodes_"      +std::to_string(i_nnodes)
                    +"_niters_"      +std::to_string(i_niters);
    write_file(*i_ctx, vort_Y, filename.c_str());

    filename = "prog_jump_div_current_proc_"+std::to_string(i_current_proc)
                    +"_current_step_"+std::to_string(i_current_step)
                    +"_current_iter_"+std::to_string(i_current_iter)
                    +"_nnodes_"      +std::to_string(i_nnodes)
                    +"_niters_"      +std::to_string(i_niters);
    write_file(*i_ctx, div_Y, filename.c_str());

  }



void cecho_output_solution(
			     PlaneDataCtx *i_ctx,
			     PlaneDataVars *i_Y,
			     int i_current_proc,
			     int i_current_step,
			     int i_current_iter,
			     int i_nnodes,
			     int i_niters
			     )
  {
    const PlaneData& h_Y    = i_Y->get_h();
    const PlaneData& vort_Y = i_Y->get_vort();
    const PlaneData& div_Y  = i_Y->get_div();

    // get the SimulationVariables object from context
    SimulationVariables* simVars(i_ctx->get_simulation_variables());

    // write the data to file
    std::string filename = "prog_h_current_proc_"+std::to_string(i_current_proc)
                                +"_current_step_"+std::to_string(i_current_step)
                                +"_current_iter_"+std::to_string(i_current_iter)
                                +"_nnodes_"      +std::to_string(i_nnodes)
                                +"_niters_"      +std::to_string(i_niters);
    write_file(*i_ctx, h_Y, filename.c_str());

    filename = "prog_vort_current_proc_"+std::to_string(i_current_proc)
                    +"_current_step_"+std::to_string(i_current_step)
                    +"_current_iter_"+std::to_string(i_current_iter)
                    +"_nnodes_"      +std::to_string(i_nnodes)
                    +"_niters_"      +std::to_string(i_niters);
    write_file(*i_ctx, vort_Y, filename.c_str());

    filename = "prog_div_current_proc_"+std::to_string(i_current_proc)
                    +"_current_step_"+std::to_string(i_current_step)
                    +"_current_iter_"+std::to_string(i_current_iter)
                    +"_nnodes_"      +std::to_string(i_nnodes)
                    +"_niters_"      +std::to_string(i_niters);
    write_file(*i_ctx, div_Y, filename.c_str());

  }

}
