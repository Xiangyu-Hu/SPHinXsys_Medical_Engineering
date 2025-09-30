#include "bidirectional_boundary_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
BufferIndicationCK::BufferIndicationCK(AlignedBoxByCell &aligned_box_part)
    : BaseLocalDynamics<AlignedBoxByCell>(aligned_box_part),
      part_id_(aligned_box_part.getPartID()),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_buffer_indicator_(particles_->registerStateVariable<int>("BufferIndicator"))
{
    particles_->addEvolvingVariable<int>("BufferIndicator");
}
//=================================================================================================//
BufferOutflowIndication::BufferOutflowIndication(AlignedBoxByCell &aligned_box_part)
    : BaseLocalDynamics<AlignedBoxByCell>(aligned_box_part),
      part_id_(aligned_box_part.getPartID()),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      sv_total_real_particles_(particles_->svTotalRealParticles()),
      dv_buffer_indicator_(
          particles_->registerStateVariable<int>("BufferIndicator")),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_life_status_(particles_->registerStateVariable<int>("LifeStatus", 0)) {}
//=================================================================================================//
BufferOutflowIndication::UpdateKernel::
    IsDeletable::IsDeletable(int part_id, AlignedBox *aligned_box,
                             Vecd *pos, int *buffer_particle_indicator)
    : part_id_(part_id), aligned_box_(aligned_box), pos_(pos),
      buffer_indicator_(buffer_particle_indicator) {}
//=================================================================================================//
OutflowParticleDeletion::OutflowParticleDeletion(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      remove_real_particle_method_(particles_),
      dv_life_status_(particles_->getVariableByName<int>("LifeStatus")) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
