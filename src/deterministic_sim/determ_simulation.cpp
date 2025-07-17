#include "openmc/deterministic_sim/deterministic.h"

#include "openmc/deterministic_sim/deterministic.h"
#include "openmc/eigenvalue.h"
#include "openmc/geometry.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/output.h"
#include "openmc/plot.h"
#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/random_ray/random_ray.h"
#include "openmc/simulation.h"
#include "openmc/source.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/tally_scoring.h"
#include "openmc/timer.h"
#include "openmc/weight_windows.h"

namespace openmc {

//==============================================================================
// Non-member functions
//==============================================================================

// logic has to be to connect tally objects with source(s) and integrate over
// them source regions need

int openmc_run_determ()
{
  int err = 0;

  openmc::simulation::time_total.start();

  // intitialize structures
  openmc_simulation_init();

  openmc_sim_det_init();

  return err;
}

// Enforces restrictions on inputs in deterministic mode.  While there are
// many features that don't make sense in deterministic mode, and are therefore
// unsupported, we limit our testing/enforcement operations only to inputs
// that may cause erroneous/misleading output or crashes from the solver.
void validate_determ_input()
{

  for (auto& tally : model::tallies) {

    // Validate score types
    for (auto score_bin : tally->scores_) {
      switch (score_bin) {
      case SCORE_FLUX:
      case SCORE_TOTAL:
      case SCORE_FISSION:
      case SCORE_NU_FISSION:
      case SCORE_EVENTS:
        break;
      default:
        fatal_error(
          "Invalid score specified. Only flux, total, fission, nu-fission, and "
          "event scores are supported in random ray mode.");
      }
    }

    // Validate filter types
    for (auto f : tally->filters()) {
      auto& filter = *model::tally_filters[f];

      switch (filter.type()) {
      case FilterType::CELL:
      case FilterType::CELL_INSTANCE:
      case FilterType::DISTRIBCELL:
      case FilterType::ENERGY:
      case FilterType::MATERIAL:
      case FilterType::MESH:
      case FilterType::UNIVERSE:
      case FilterType::PARTICLE:
        break;
      default:
        fatal_error("Invalid filter specified. Only cell, cell_instance, "
                    "distribcell, energy, material, mesh, and universe filters "
                    "are supported in random ray mode.");
      }
    }
  }
}

void openmc_sim_det_init()
{
  std::cout << "Nothing to do here yet" << std::endl;
}

void DetermSimulation::dummy_initialize_targets()
{
  // this is test_data
  Position ll {-10, -10, 0};
  cont int N {10};
  Position ur {10, 10, 0};
  Position dd = (ur - ll) = / N;

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      targets_.push_back(
        SourceSite(ll.x + dd.x * i, ll.x + dd.y * j, ll.z + dd.z * 0));
    }
  }
}

void DetermSimulation::initialize_targets()
{
  dummy_initialize_targets();
}

void DetermSimulations::initialize_sources()
{
  dummy_initialize_soruces();
}

void DetermSimulation::dummy_initialize_sources()
{
  // this is test_data
  SourceSite s0, s1;
  s0.x = -2.134;
  s0.y = -2.134;
  s1.x = 2.134;
  s1.y = 2.134;
  sources_.push_back(s0);
  sources_.push_back(s1);
}

void DetermSimulation::line_source_element()
{
  // return a
}

DetermSimulation::simulate()
{

  while (simulation::current_batch < settings::n_batches) {
    initialize_batch();
    initialize_generation();

#pragma omp parallel
    for (i = 0; i < settings::n_particles; ii) {
      DetermRay ray(i);
      ray.transport_history_based_single_ray();
    } // rays in batch

  } // batches
}

DetermRay::DetermRay(uint64_t ray_id, SourceSite ss, SourceSite tgt)
{
  // find direction
  Position uu = tgt.u() - ss.u();
  // normalize
  uu = uu / u.dot(u);
}

void DetermRay::transport_history_based_single_ray() {} // namespace openmc
