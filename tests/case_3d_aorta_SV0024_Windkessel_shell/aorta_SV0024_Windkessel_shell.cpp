/**
 * @file 	aorta_SV0024_Windkessel_shell.cpp
 * @brief   Aorta with shell, imposed velocity inlet and windkessel outlet.
 */

#include "sphinxsys.h"
#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "particle_generation_and_detection.h"
#include "windkessel_bc.h"
#include "hemodynamic_indices.h"

using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
std::string full_path_to_file = "./input/aorta_0154_0001.stl";
BoundingBox system_domain_bounds(Vecd(-6.0E-2, -3.0E-2, -1.0E-2), Vecd(6.0E-2, 8.0E-2, 20.0E-2));
//----------------------------------------------------------------------
//	Global parameters on the fluid and solid properties.
//----------------------------------------------------------------------
Real rho0_f = 1060.0;                   
Real mu_f = 0.0035;
Real U_f = 3.0;
Real c_f = 10.0*U_f;
Real scaling = 1.0E-2;
Real dp_0 = 0.06 * scaling;
Vecd translation(0.0, 0.0, 0.0);
Real shell_resolution = dp_0;
Real thickness = 0.25 * scaling;
Real rho0_s = 1000;                /** Normalized density. */
Real Youngs_modulus = 7.5e5;    /** Normalized Youngs Modulus. */
Real poisson = 0.49;               /** Poisson ratio. */
Real physical_viscosity = 2000;
//----------------------------------------------------------------------
//	Define buffer locations.
//----------------------------------------------------------------------
struct RotationResult
{
    Vec3d axis;
    Real angle;
};

RotationResult RotationCalculator(Vecd target_normal, Vecd standard_direction)
{
    target_normal.normalize();

    Vec3d axis = standard_direction.cross(target_normal);
    Real angle = std::acos(standard_direction.dot(target_normal));

    if (axis.norm() < 1e-6)
    {
        if (standard_direction.dot(target_normal) < 0)
        {
            axis = Vec3d(1, 0, 0);
            angle = M_PI;
        }
        else
        {
            axis = Vec3d(0, 0, 1);
            angle = 0;
        }
    }
    else
    {
        axis.normalize();
    }

    return {axis, angle};
}

Vecd standard_direction(1, 0, 0);

// inlet (-1.015, 4.519, 2.719), (0.100, 0.167, 0.981), 1.366
Real A_in = 5.9825 * scaling * scaling;
Real radius_inlet = 1.366 * scaling;
Vec3d inlet_half = Vec3d(2.0 * dp_0, 1.6 * scaling, 1.6 * scaling);
Vec3d inlet_vector(0.100, 0.167, 0.981);
Vec3d inlet_normal = inlet_vector.normalized();
Vec3d inlet_center = Vec3d(-1.015, 4.519, 2.719) * scaling;
Vec3d inlet_cut_translation = inlet_center - inlet_normal * 1.5 * dp_0;
Vec3d inlet_blood_cut_translation = inlet_center;
Vec3d inlet_buffer_translation = inlet_center + inlet_normal * 4.0 * dp_0;
RotationResult inlet_rotation_result = RotationCalculator(inlet_normal, standard_direction);
Rotation3d inlet_emitter_rotation(inlet_rotation_result.angle, inlet_rotation_result.axis);
Rotation3d inlet_disposer_rotation(inlet_rotation_result.angle + Pi, inlet_rotation_result.axis);

// outlet1 (-0.947, 4.322, 10.110), (0.607, 0.373, 0.702), 0.280
Real A_out1 = 0.2798 * scaling * scaling;
Vec3d outlet_1_half = Vec3d(2.0 * dp_0, 0.5 * scaling, 0.5 * scaling);
Vec3d outlet_1_vector(0.607, 0.373, 0.702);
Vec3d outlet_1_normal = outlet_1_vector.normalized();
Vec3d outlet_1_center = Vec3d(-0.947, 4.322, 10.110) * scaling;
Vec3d outlet_1_cut_translation = outlet_1_center + outlet_1_normal * 1.5 * dp_0;
Vec3d outlet_1_blood_cut_translation = outlet_1_center;
Vec3d outlet_1_buffer_translation = outlet_1_center - outlet_1_normal * 4.0 * dp_0;
RotationResult outlet_1_rotation_result = RotationCalculator(outlet_1_normal, standard_direction);
Rotation3d outlet_1_disposer_rotation(outlet_1_rotation_result.angle, outlet_1_rotation_result.axis);
Rotation3d outlet_1_emitter_rotation(outlet_1_rotation_result.angle + Pi, outlet_1_rotation_result.axis);

// outlet2 (-2.650, 3.069, 10.911), (-0.099, 0.048, 0.994), 0.253
Real A_out2 = 0.2043 * scaling * scaling;
Vec3d outlet_2_half = Vec3d(2.0 * dp_0, 1.0 * scaling, 1.0 * scaling);
Vec3d outlet_2_vector(-0.099, 0.048, 0.994);
Vec3d outlet_2_normal = outlet_2_vector.normalized();
Vec3d outlet_2_center = Vec3d(-2.650, 3.069, 10.911) * scaling;
Vec3d outlet_2_cut_translation = outlet_2_center + outlet_2_normal * 1.5 * dp_0;
Vec3d outlet_2_blood_cut_translation = outlet_2_center;
Vec3d outlet_2_buffer_translation = outlet_2_center - outlet_2_normal * 4.0 * dp_0;
RotationResult outlet_2_rotation_result = RotationCalculator(outlet_2_normal, standard_direction);
Rotation3d outlet_2_disposer_rotation(outlet_2_rotation_result.angle, outlet_2_rotation_result.axis);
Rotation3d outlet_2_emitter_rotation(outlet_2_rotation_result.angle + Pi, outlet_2_rotation_result.axis);

// outlet3 (-2.946, 1.789, 10.007), (-0.146, -0.174, 0.974), 0.356
Real A_out3 = 0.3721 * scaling * scaling;
Real radius_outlet3 = 0.356 * scaling;
Vec3d outlet_3_half = Vec3d(2.0 * dp_0, 0.4 * scaling, 0.4 * scaling);
Vec3d outlet_3_vector(-0.146, -0.174, 0.974);
Vec3d outlet_3_normal = outlet_3_vector.normalized();
Vec3d outlet_3_center = Vec3d(-2.946, 1.789, 10.007) * scaling;
Vec3d outlet_3_cut_translation = outlet_3_center + outlet_3_normal * 1.5 * dp_0;
Vec3d outlet_3_blood_cut_translation = outlet_3_center;
Vec3d outlet_3_buffer_translation = outlet_3_center - outlet_3_normal * 4.0 * dp_0;
RotationResult outlet_3_rotation_result = RotationCalculator(outlet_3_normal, standard_direction);
Rotation3d outlet_3_disposer_rotation(outlet_3_rotation_result.angle, outlet_3_rotation_result.axis);
Rotation3d outlet_3_emitter_rotation(outlet_3_rotation_result.angle + Pi, outlet_3_rotation_result.axis);

// outlet4 (-1.052, 1.152, 9.669), (0.568, 0.428, 0.703), 0.395
Real A_out4 = 0.4567 * scaling * scaling;
Vec3d outlet_4_half = Vec3d(2.0 * dp_0, 1.0 * scaling, 1.0 * scaling);
Vec3d outlet_4_vector(0.568, 0.428, 0.703);
Vec3d outlet_4_normal = outlet_4_vector.normalized();
Vec3d outlet_4_center = Vec3d(-1.052, 1.152, 9.669) * scaling;
Vec3d outlet_4_cut_translation = outlet_4_center + outlet_4_normal * 1.5 * dp_0;
Vec3d outlet_4_blood_cut_translation = outlet_4_center;
Vec3d outlet_4_buffer_translation = outlet_4_center - outlet_4_normal * 4.0 * dp_0;
RotationResult outlet_4_rotation_result = RotationCalculator(outlet_4_normal, standard_direction);
Rotation3d outlet_4_disposer_rotation(outlet_4_rotation_result.angle, outlet_4_rotation_result.axis);
Rotation3d outlet_4_emitter_rotation(outlet_4_rotation_result.angle + Pi, outlet_4_rotation_result.axis);

// outlet5 (-1.589, -0.797, 0.247), (-0.033, 0.073, -0.997), 0.921
Real A_out5 = 2.6756 * scaling * scaling;
Vec3d outlet_5_half = Vec3d(2.0 * dp_0, 1.5 * scaling, 1.5 * scaling);
Vec3d outlet_5_vector(-0.033, 0.073, -0.997);
Vec3d outlet_5_normal = outlet_5_vector.normalized();
Vec3d outlet_5_center = Vec3d(-1.589, -0.797, 0.247) * scaling;
Vec3d outlet_5_cut_translation = outlet_5_center + outlet_5_normal * 1.5 * dp_0;
Vec3d outlet_5_blood_cut_translation = outlet_5_center;
Vec3d outlet_5_buffer_translation = outlet_5_center - outlet_5_normal * 4.0 * dp_0;
RotationResult outlet_5_rotation_result = RotationCalculator(outlet_5_normal, standard_direction);
Rotation3d outlet_5_disposer_rotation(outlet_5_rotation_result.angle, outlet_5_rotation_result.axis);
Rotation3d outlet_5_emitter_rotation(outlet_5_rotation_result.angle + Pi, outlet_5_rotation_result.axis);

//----------------------------------------------------------------------
//	Define the imported model.
//----------------------------------------------------------------------
namespace SPH
{
class ShellShape : public ComplexShape
{
public:
    explicit ShellShape(const std::string &shape_name) : ComplexShape(shape_name),
        mesh_shape_(new TriangleMeshShapeSTL(full_path_to_file, translation, scaling))
    {
        add<ExtrudeShape<TriangleMeshShapeSTL>>(shell_resolution/2, full_path_to_file, translation, scaling);
    }

    TriangleMeshShapeSTL* getMeshShape() const
    {
        return mesh_shape_.get();
    }

private:
    std::unique_ptr<TriangleMeshShapeSTL> mesh_shape_;
};

class WaterBlock : public ComplexShape
{
public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(full_path_to_file, translation, scaling);
    }
};
//----------------------------------------------------------------------
//	Shell particle generation.
//----------------------------------------------------------------------
class FromSTLFile;
template <>
class ParticleGenerator<SurfaceParticles, FromSTLFile> : public ParticleGenerator<SurfaceParticles>
{
    Real mesh_total_area_;
    Real particle_spacing_;
    const Real thickness_;
    Real avg_particle_volume_;
    size_t planned_number_of_particles_;

    TriangleMeshShapeSTL* mesh_shape_;
    Shape &initial_shape_;

public:
    explicit ParticleGenerator(SPHBody& sph_body, SurfaceParticles &surface_particles, TriangleMeshShapeSTL* mesh_shape, Real shell_thickness) 
        : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles),
        mesh_total_area_(0),
        particle_spacing_(sph_body.getSPHAdaptation().ReferenceSpacing()),
        thickness_(shell_thickness),
        avg_particle_volume_(pow(particle_spacing_, Dimensions - 1) * thickness_),
        planned_number_of_particles_(0),
        mesh_shape_(mesh_shape), initial_shape_(sph_body.getInitialShape()) 
    {
        if (!mesh_shape_)
        {
            std::cerr << "Error: Mesh shape is not set!" << std::endl;
            return;
        }

        if (!initial_shape_.isValid())
        {
            std::cout << "\n BaseParticleGeneratorLattice Error: initial_shape_ is invalid." << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            throw;
        }
    }

    virtual void prepareGeometricData() override
    {
        
        // Preload vertex positions
        StdVec<std::array<Real, 3>> vertex_positions;
        auto &V = mesh_shape_->getVertices();
        int num_vertices = static_cast<int>(V.size());
        vertex_positions.reserve(num_vertices);
        for (int i = 0; i < num_vertices; i++)
        {
            const auto &p = V[i];
            vertex_positions.push_back({Real(p[0]), Real(p[1]), Real(p[2])});
        }

        // Preload face
        StdVec<std::array<int, 3>> faces;
        auto &F = mesh_shape_->getFaces();
        int num_faces = static_cast<int>(F.size());
        std::cout << "num_faces calculation = " << num_faces << std::endl;
        faces.reserve(num_faces);
        for (int i = 0; i < num_faces; i++)
        {
            faces.push_back({ F[i][0], F[i][1], F[i][2] });
        }

        // Calculate total volume
        std::vector<Real> face_areas(num_faces);
        for (int i = 0; i < num_faces; ++i)
        {
            Vec3d vertices[3];
            for (int j = 0; j < 3; ++j)
            {
                const auto& pos = vertex_positions[faces[i][j]];
                vertices[j] = Vec3d(pos[0], pos[1], pos[2]);
            }

            Real each_area = calculateEachFaceArea(vertices);
            face_areas[i] = each_area;
            mesh_total_area_ += each_area;
        }

        Real number_of_particles = mesh_total_area_ * thickness_ / avg_particle_volume_ + 0.5;
        planned_number_of_particles_ = int(number_of_particles);
        std::cout << "planned_number_of_particles calculation = " << planned_number_of_particles_ << std::endl;

        // initialize a uniform distribution between 0 (inclusive) and 1 (exclusive)
        std::mt19937_64 rng;
        std::uniform_real_distribution<Real> unif(0, 1);

        // Calculate the interval based on the number of particles.
        Real interval = planned_number_of_particles_ / (num_faces + TinyReal);  // if planned_number_of_particles_ >= num_faces, every face will generate particles
        if (interval <= 0)
            interval = 1; // It has to be lager than 0.

        for (int i = 0; i < num_faces; ++i)
        {
            Vec3d vertices[3];
            for (int j = 0; j < 3; ++j)
            {
                const auto& pos = vertex_positions[faces[i][j]];
                vertices[j] = Vec3d(pos[0], pos[1], pos[2]);
            }

            Real random_real = unif(rng);
            if (random_real <= interval && base_particles_.TotalRealParticles() < planned_number_of_particles_)
            {
                // Generate particle at the center of this triangle face
                // generateParticleAtFaceCenter(vertices);
                
                // Generate particles on this triangle face, unequal
                int particles_per_face = std::max(1, int(planned_number_of_particles_ * (face_areas[i] / mesh_total_area_)));
                generateParticlesOnFace(vertices, particles_per_face);
            }
        }
    }

private:
    Real calculateEachFaceArea(const Vec3d vertices[3])
    {
        Vec3d edge1 = vertices[1] - vertices[0];
        Vec3d edge2 = vertices[2] - vertices[0];
        Real area = 0.5 * edge1.cross(edge2).norm();
        return area;
    }

    void generateParticleAtFaceCenter(const Vec3d vertices[3])
    {
        Vec3d face_center = (vertices[0] + vertices[1] + vertices[2]) / 3.0;

        addPositionAndVolumetricMeasure(face_center, avg_particle_volume_/thickness_);
        addSurfaceProperties(initial_shape_.findNormalDirection(face_center), thickness_);
    }

    void generateParticlesOnFace(const Vec3d vertices[3], int particles_per_face)
    {
        for (int k = 0; k < particles_per_face; ++k)
        {
            Real u = static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX);
            Real v = static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX);

            if (u + v > 1.0) {
                u = 1.0 - u;
                v = 1.0 - v;
            }

            Vec3d particle_position = (1 - u - v) * vertices[0] + u * vertices[1] + v * vertices[2];
            Vec3d normal = initial_shape_.findNormalDirection(particle_position);
            Vec3d offset_position = particle_position + 0.5 * particle_spacing_ * normal;

            addPositionAndVolumetricMeasure(offset_position, avg_particle_volume_ / thickness_);
            addSurfaceProperties(normal, thickness_);
        }
    }
};

struct InflowVelocity
{
    Real u_ave, interval_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ave(0.0), interval_(0.66) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = velocity;
        int n = static_cast<int>(current_time / interval_);
        Real t_in_cycle = current_time - n * interval_;

        u_ave = 5.0487;
        Real a[8] = {4.5287, -4.3509, -5.8551, -1.5063, 1.2800, 0.9012, 0.0855, -0.0480};
        Real b[8] = {-8.0420, -6.2637, 0.7465, 3.5239, 1.6283, -0.1306, -0.2738, -0.0449};

        Real w = 2 * Pi / 1.0;
        for (size_t i = 0; i < 8; i++)
        {
            u_ave = u_ave + a[i] * cos(w * (i + 1) * t_in_cycle) + b[i] * sin(w * (i + 1) * t_in_cycle);
        }
        target_velocity[0] = 2.0 * u_ave * (1.0 - (position[1] * position[1] + position[2] * position[2]) / radius_inlet / radius_inlet);
        target_velocity[1] = 0.0;
        target_velocity[2] = 0.0;    

        return target_velocity;
    }
};
//----------------------------------------------------------------------
//	Boundary constrain
//----------------------------------------------------------------------
class BoundaryGeometry : public BodyPartByParticle
{
  public:
    BoundaryGeometry(SPHBody &body, const std::string &body_part_name) : BodyPartByParticle(body)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };
    virtual ~BoundaryGeometry() {};

  private:
    bool tagManually(size_t index_i)
    {
        const Vecd &pos = base_particles_.ParticlePositions()[index_i];

        Vecd relative_position_inlet = pos - inlet_center;
        Vecd relative_position_outlet_01 = pos - outlet_1_center;
        Vecd relative_position_outlet_02 = pos - outlet_2_center;
        Vecd relative_position_outlet_03 = pos - outlet_3_center;
        Vecd relative_position_outlet_04 = pos - outlet_4_center;
        Vecd relative_position_outlet_05 = pos - outlet_5_center;

        Real proj_inlet = relative_position_inlet.dot(inlet_normal);
        Real proj_out_01 = relative_position_outlet_01.dot(-outlet_1_normal);
        Real proj_out_02 = relative_position_outlet_02.dot(-outlet_2_normal);
        Real proj_out_03 = relative_position_outlet_03.dot(-outlet_3_normal);
        Real proj_out_04 = relative_position_outlet_04.dot(-outlet_4_normal);
        Real proj_out_05 = relative_position_outlet_05.dot(-outlet_5_normal);

        bool is_inlet =
            (ABS(proj_inlet)  < 4 * dp_0) && (pos[1] > 0.) && (pos[2] < 0.033);
        bool is_outlet_01 =
            (ABS(proj_out_01) < 4 * dp_0) && (pos[1] > 0.038) && (pos[2] > 0.05);
        bool is_outlet_02 =
            (ABS(proj_out_02) < 4 * dp_0) && (pos[2] > 0.106);
        bool is_outlet_03 =
            (ABS(proj_out_03) < 4 * dp_0) && (pos[1] > 0.014) && (pos[1] < 0.023) && (pos[0] < -0.02);
        bool is_outlet_04 =
            (ABS(proj_out_04) < 4 * dp_0) && (pos[1] < 0.0147) && (pos[0] > -0.015) && (pos[2] > 0.);
        bool is_outlet_05 =
            (ABS(proj_out_05) < 4 * dp_0) && (pos[2] < 0.01);

        bool tagged = is_inlet || is_outlet_01 || is_outlet_02
                                         || is_outlet_03 || is_outlet_04 || is_outlet_05;

        return tagged;
    };
};
} // namespace SPH

//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, dp_0);
    sph_system.setRunParticleRelaxation(false); // Tag for run particle relaxation for body-fitted distribution
    sph_system.setReloadParticles(true);       // Tag for computation with save particles distribution
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(0.5);
    if (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
    {
        water_block.generateParticlesWithReserve<BaseParticles, Reload>(in_outlet_particle_buffer, water_block.getName());
    }
    else
    {
        water_block.defineBodyLevelSetShape()->correctLevelSetSign();
        water_block.generateParticles<BaseParticles, Lattice>();
    }
    ShellShape body_from_mesh("BodyFromMesh");
    TriangleMeshShapeSTL* mesh_shape = body_from_mesh.getMeshShape();
    SolidBody shell_body(sph_system, makeShared<ShellShape>("ShellBody"));
    shell_body.defineAdaptation<SPHAdaptation>(1.15, dp_0/shell_resolution);
    shell_body.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    if (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
    {
        shell_body.generateParticles<SurfaceParticles, Reload>(shell_body.getName());
    }
    else
    {
        shell_body.defineBodyLevelSetShape(2.0)->correctLevelSetSign();
        shell_body.generateParticles<SurfaceParticles, FromSTLFile>(mesh_shape, thickness);
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation shell_inner(shell_body);
    ContactRelationFromShellToFluid water_shell_contact(water_block, {&shell_body}, {false});
    ContactRelationFromFluidToShell shell_water_contact(shell_body, {&water_block}, {false});
    ShellInnerRelationWithContactKernel shell_curvature_inner(shell_body, water_block);
    ComplexRelation water_block_complex(water_block_inner, {&water_shell_contact});

    //----------------------------------------------------------------------
    //	SPH Particle relaxation section
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        InnerRelation blood_inner(water_block);

        AlignedBoxByCell inlet_detection_box(shell_body, AlignedBox(xAxis, Transform(Rotation3d(inlet_emitter_rotation), Vec3d(inlet_cut_translation)), inlet_half));
        AlignedBoxByCell outlet01_detection_box(shell_body, AlignedBox(xAxis, Transform(Rotation3d(outlet_1_emitter_rotation), Vec3d(outlet_1_cut_translation)), outlet_1_half));
        AlignedBoxByCell outlet02_detection_box(shell_body, AlignedBox(xAxis, Transform(Rotation3d(outlet_2_emitter_rotation), Vec3d(outlet_2_cut_translation)), outlet_2_half));
        AlignedBoxByCell outlet03_detection_box(shell_body, AlignedBox(xAxis, Transform(Rotation3d(outlet_3_emitter_rotation), Vec3d(outlet_3_cut_translation)), outlet_3_half));
        //AlignedBoxByCell outlet04_detection_cylinder(shell_body, AlignedBox(xAxis, Transform(Rotation3d(outlet_4_emitter_rotation), Vec3d(outlet_4_cut_translation)), outlet_4_half[1], outlet_4_half[0]));
        AlignedBoxByCell outlet04_detection_cylinder(shell_body, AlignedBox(xAxis, Transform(Rotation3d(outlet_4_emitter_rotation), Vec3d(outlet_4_cut_translation)), outlet_4_half));
        AlignedBoxByCell outlet05_detection_box(shell_body, AlignedBox(xAxis, Transform(Rotation3d(outlet_5_emitter_rotation), Vec3d(outlet_5_cut_translation)), outlet_5_half));

        /** cut blood */
        AlignedBoxByCell inlet_blood_detection_box(water_block, AlignedBox(xAxis, Transform(Rotation3d(inlet_emitter_rotation), Vec3d(inlet_blood_cut_translation)), inlet_half));
        AlignedBoxByCell outlet01_blood_detection_box(water_block, AlignedBox(xAxis, Transform(Rotation3d(outlet_1_emitter_rotation), Vec3d(outlet_1_blood_cut_translation)), outlet_1_half));
        AlignedBoxByCell outlet02_blood_detection_box(water_block, AlignedBox(xAxis, Transform(Rotation3d(outlet_2_emitter_rotation), Vec3d(outlet_2_blood_cut_translation)), outlet_2_half));
        AlignedBoxByCell outlet03_blood_detection_box(water_block, AlignedBox(xAxis, Transform(Rotation3d(outlet_3_emitter_rotation), Vec3d(outlet_3_blood_cut_translation)), outlet_3_half));
        AlignedBoxByCell outlet04_blood_detection_box(water_block, AlignedBox(xAxis, Transform(Rotation3d(outlet_4_emitter_rotation), Vec3d(outlet_4_blood_cut_translation)), outlet_4_half));
        AlignedBoxByCell outlet05_blood_detection_box(water_block, AlignedBox(xAxis, Transform(Rotation3d(outlet_5_emitter_rotation), Vec3d(outlet_5_blood_cut_translation)), outlet_5_half));

        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        /** A  Physics relaxation step. */
        SurfaceRelaxationStep relaxation_step_inner(shell_inner);
        ShellNormalDirectionPrediction shell_normal_prediction(shell_inner, shell_resolution * 1.0);

        //RelaxationStepInner relaxation_step_inner_blood(blood_inner);
        RelaxationStepComplex relaxation_step_blood(blood_inner, water_shell_contact);

        // here, need a class to switch particles in aligned box to ghost particles (not real particles)
        SimpleDynamics<ParticlesInAlignedBoxDetectionByCell> inlet_particles_detection(inlet_detection_box);
        SimpleDynamics<ParticlesInAlignedBoxDetectionByCell> outlet01_particles_detection(outlet01_detection_box);
        SimpleDynamics<ParticlesInAlignedBoxDetectionByCell> outlet02_particles_detection(outlet02_detection_box);
        SimpleDynamics<ParticlesInAlignedBoxDetectionByCell> outlet03_particles_detection(outlet03_detection_box);
        SimpleDynamics<ParticlesInAlignedBoxDetectionByCell> outlet04_particles_detection(outlet04_detection_cylinder);
        SimpleDynamics<ParticlesInAlignedBoxDetectionByCell> outlet05_particles_detection(outlet05_detection_box);

        SimpleDynamics<ParticlesInAlignedBoxDetectionByCell> inlet_blood_particles_detection(inlet_blood_detection_box);
        SimpleDynamics<ParticlesInAlignedBoxDetectionByCell> outlet01_blood_particles_detection(outlet01_blood_detection_box);
        SimpleDynamics<ParticlesInAlignedBoxDetectionByCell> outlet02_blood_particles_detection(outlet02_blood_detection_box);
        SimpleDynamics<ParticlesInAlignedBoxDetectionByCell> outlet03_blood_particles_detection(outlet03_blood_detection_box);
        SimpleDynamics<ParticlesInAlignedBoxDetectionByCell> outlet04_blood_particles_detection(outlet04_blood_detection_box);
        SimpleDynamics<ParticlesInAlignedBoxDetectionByCell> outlet05_blood_particles_detection(outlet05_blood_detection_box);

        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_shell_to_vtp({shell_body});
        write_shell_to_vtp.addToWrite<Vecd>(shell_body, "NormalDirection");
        BodyStatesRecordingToVtp write_blood_to_vtp({water_block});
        BodyStatesRecordingToVtp write_all_bodies_to_vtp({sph_system});
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files({ &shell_body, &water_block });

        ParticleSorting particle_sorting_shell(shell_body);
        ParticleSorting particle_sorting_water(water_block);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        relaxation_step_inner.getOnSurfaceBounding().exec();
        relaxation_step_blood.SurfaceBounding().exec();
        write_shell_to_vtp.writeToFile(0.0);
        write_blood_to_vtp.writeToFile(0.0);
        shell_body.updateCellLinkedList();
        water_block.updateCellLinkedList();
        //----------------------------------------------------------------------
        //	Particle relaxation time stepping start here.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 5000)
        {
            relaxation_step_inner.exec();
            relaxation_step_blood.exec();
            ite_p += 1;
            if (ite_p % 500 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the imported model N = " << ite_p << "\n";

                if (ite_p % 1000 == 0)
                {
                    write_shell_to_vtp.writeToFile(ite_p);
                    write_blood_to_vtp.writeToFile(ite_p);
                }
            }
        }
        std::cout << "The physics relaxation process of imported model finish !" << std::endl;

        shell_normal_prediction.smoothing_normal_exec();

        inlet_particles_detection.exec();
        particle_sorting_shell.exec();
        shell_body.updateCellLinkedList();
        outlet01_particles_detection.exec();
        particle_sorting_shell.exec();
        shell_body.updateCellLinkedList();
        outlet02_particles_detection.exec();
        particle_sorting_shell.exec();
        shell_body.updateCellLinkedList();
        outlet03_particles_detection.exec();
        particle_sorting_shell.exec();
        shell_body.updateCellLinkedList();
        outlet04_particles_detection.exec();
        particle_sorting_shell.exec();
        shell_body.updateCellLinkedList();
        outlet05_particles_detection.exec();
        particle_sorting_shell.exec();
        shell_body.updateCellLinkedList();

        inlet_blood_particles_detection.exec();
        particle_sorting_water.exec();
        water_block.updateCellLinkedList();
        outlet01_blood_particles_detection.exec();
        particle_sorting_water.exec();
        water_block.updateCellLinkedList();
        outlet02_blood_particles_detection.exec();
        particle_sorting_water.exec();
        water_block.updateCellLinkedList();
        outlet03_blood_particles_detection.exec();
        particle_sorting_water.exec();
        water_block.updateCellLinkedList();
        outlet04_blood_particles_detection.exec();
        particle_sorting_water.exec();
        water_block.updateCellLinkedList();
        outlet05_blood_particles_detection.exec();
        particle_sorting_water.exec();
        water_block.updateCellLinkedList();

        write_all_bodies_to_vtp.writeToFile(ite_p);
        write_particle_reload_files.writeToFile(0);

        return 0;
    }

    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------

    //----------------------------------------------------------------------
    //	Solid dynamics
    //----------------------------------------------------------------------
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> shell_corrected_configuration(shell_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> shell_stress_relaxation_first(shell_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> shell_stress_relaxation_second(shell_inner);
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> shell_time_step_size(shell_body);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> shell_average_curvature(shell_curvature_inner);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> shell_update_normal(shell_body);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>
        shell_velocity_damping(0.5, shell_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>
        shell_rotation_damping(0.5, shell_inner, "AngularVelocity", physical_viscosity);

    /** Exert constrain on shell. */
    BoundaryGeometry boundary_geometry(shell_body, "BoundaryGeometry");
    SimpleDynamics<FixBodyPartConstraint> constrain_holder(boundary_geometry);
    
    //----------------------------------------------------------------------
    //	Fluid dynamics
    //----------------------------------------------------------------------
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> kernel_correction_complex(water_block_inner, water_shell_contact);
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_shell_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex>
        boundary_indicator(water_block_inner, water_shell_contact);
    /** Pressure relaxation algorithm without Riemann solver for viscous flows. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionWithWallRiemann> pressure_relaxation(water_block_inner, water_shell_contact);
    /** Pressure relaxation algorithm by using position verlet time stepping. */
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_shell_contact);
    /* Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_f);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);
    /** Computing viscous acceleration. */
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(water_block_inner, water_shell_contact);
    /** Impose transport velocity. */
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>>
        transport_velocity_correction(water_block_inner, water_shell_contact);
    
    // add buffers and boundary conditions
    AlignedBoxByCell inlet_emitter(water_block, AlignedBox(xAxis, Transform(Rotation3d(inlet_emitter_rotation), Vec3d(inlet_buffer_translation)), inlet_half));
    fluid_dynamics::BidirectionalBuffer<fluid_dynamics::NonPrescribedPressure> inflow_injection(inlet_emitter, in_outlet_particle_buffer);
    AlignedBoxByCell outflow_emitter_1(water_block, AlignedBox(xAxis, Transform(Rotation3d(outlet_1_emitter_rotation), Vec3d(outlet_1_buffer_translation)), outlet_1_half));
    fluid_dynamics::WindkesselOutletBidirectionalBuffer outflow_injection_1(outflow_emitter_1, in_outlet_particle_buffer);
    AlignedBoxByCell outflow_emitter_2(water_block, AlignedBox(xAxis, Transform(Rotation3d(outlet_2_emitter_rotation), Vec3d(outlet_2_buffer_translation)), outlet_2_half));
    fluid_dynamics::WindkesselOutletBidirectionalBuffer outflow_injection_2(outflow_emitter_2, in_outlet_particle_buffer);
    AlignedBoxByCell outflow_emitter_3(water_block, AlignedBox(xAxis, Transform(Rotation3d(outlet_3_emitter_rotation), Vec3d(outlet_3_buffer_translation)), outlet_3_half));
    fluid_dynamics::WindkesselOutletBidirectionalBuffer outflow_injection_3(outflow_emitter_3, in_outlet_particle_buffer);
    AlignedBoxByCell outflow_emitter_4(water_block, AlignedBox(xAxis, Transform(Rotation3d(outlet_4_emitter_rotation), Vec3d(outlet_4_buffer_translation)), outlet_4_half));
    fluid_dynamics::WindkesselOutletBidirectionalBuffer outflow_injection_4(outflow_emitter_4, in_outlet_particle_buffer);
    AlignedBoxByCell outflow_emitter_5(water_block, AlignedBox(xAxis, Transform(Rotation3d(outlet_5_emitter_rotation), Vec3d(outlet_5_buffer_translation)), outlet_5_half));
    fluid_dynamics::WindkesselOutletBidirectionalBuffer outflow_injection_5(outflow_emitter_5, in_outlet_particle_buffer);

    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_shell_contact);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> emitter_buffer_inflow_condition(inlet_emitter);
    SimpleDynamics<fluid_dynamics::WindkesselBoundaryCondition> outflow_pressure_condition1(outflow_emitter_1);
    SimpleDynamics<fluid_dynamics::WindkesselBoundaryCondition> outflow_pressure_condition2(outflow_emitter_2);
    SimpleDynamics<fluid_dynamics::WindkesselBoundaryCondition> outflow_pressure_condition3(outflow_emitter_3);
    SimpleDynamics<fluid_dynamics::WindkesselBoundaryCondition> outflow_pressure_condition4(outflow_emitter_4);
    SimpleDynamics<fluid_dynamics::WindkesselBoundaryCondition> outflow_pressure_condition5(outflow_emitter_5);

    ReduceDynamics<fluid_dynamics::SectionTransientFlowRate> compute_inlet_transient_flow_rate(inlet_emitter, A_in);
    //----------------------------------------------------------------------
    //	FSI
    //----------------------------------------------------------------------
    InteractionWithUpdate<solid_dynamics::WallShearStress> viscous_force_from_fluid(shell_water_contact);
    SimpleDynamics<solid_dynamics::HemodynamicIndiceCalculation> hemodynamic_indice_calculation(shell_body, 0.66);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_on_shell(shell_water_contact);
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(shell_body);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Pressure");
    body_states_recording.addToWrite<Real>(water_block, "DensityChangeRate");
    body_states_recording.addToWrite<Real>(water_block, "Density");
    body_states_recording.addToWrite<int>(water_block, "BufferIndicator");
    body_states_recording.addToWrite<Vecd>(shell_body, "NormalDirection");
    body_states_recording.addToWrite<Matd>(shell_body, "MidSurfaceCauchyStress");
    body_states_recording.addDerivedVariableRecording<SimpleDynamics<Displacement>>(shell_body);
    body_states_recording.addToWrite<Real>(shell_body, "Average1stPrincipleCurvature");
    body_states_recording.addToWrite<Real>(shell_body, "Average2ndPrincipleCurvature");
    body_states_recording.addToWrite<Vecd>(shell_body, "WallShearStress");
    body_states_recording.addToWrite<Real>(shell_body, "TimeAveragedWallShearStress");
    body_states_recording.addToWrite<Real>(shell_body, "OscillatoryShearIndex");

    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists(); 
    sph_system.initializeSystemConfigurations();
    shell_corrected_configuration.exec();
    shell_average_curvature.exec();
    constrain_holder.exec();
    water_block_complex.updateConfiguration();
    shell_water_contact.updateConfiguration();
    boundary_indicator.exec();
    inflow_injection.tag_buffer_particles.exec();
    outflow_injection_1.tag_buffer_particles.exec();
    outflow_injection_2.tag_buffer_particles.exec();
    outflow_injection_3.tag_buffer_particles.exec();
    outflow_injection_4.tag_buffer_particles.exec();
    outflow_injection_5.tag_buffer_particles.exec();
    
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0.0;
    int screen_output_interval = 10;
    Real end_time = 3.3;   /**< End time. */
    Real Output_Time = 0.01; /**< Time stamps for output of body states. */
    Real dt = 0.0;          /**< Default acoustic time step sizes. */
    Real dt_s = 0.0; /**< Default acoustic time step sizes for solid. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    Real accumulated_time = 0.006;
    int updateP_n = 0;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(0);
    //----------------------------------------------------------------------
    //	Setup Windkessel model parameters.
    //----------------------------------------------------------------------
    outflow_pressure_condition1.getTargetPressure()->setWindkesselParams(7.13E07, 8.26E-10, 1.20E09, accumulated_time);
    outflow_pressure_condition2.getTargetPressure()->setWindkesselParams(7.13E07, 8.26E-10, 1.20E09, accumulated_time);
    outflow_pressure_condition3.getTargetPressure()->setWindkesselParams(6.02E07, 9.79E-10, 1.01E09, accumulated_time);
    outflow_pressure_condition4.getTargetPressure()->setWindkesselParams(6.89E07, 8.55E-10, 1.16E09, accumulated_time);
    outflow_pressure_condition5.getTargetPressure()->setWindkesselParams(9.80E06, 6.02E-09, 1.65E08, accumulated_time);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < Output_Time)
        {  
            time_instance = TickCount::now();
            Real Dt = get_fluid_advection_time_step_size.exec();          
            update_fluid_density.exec();
            kernel_correction_complex.exec();
            viscous_acceleration.exec();
            transport_velocity_correction.exec();

            /** FSI for viscous force. */
            viscous_force_from_fluid.exec();
            hemodynamic_indice_calculation.exec(Dt);

            interval_computing_time_step += TickCount::now() - time_instance;
            /** Dynamics including pressure relaxation. */
            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);

                pressure_relaxation.exec(dt);
                /** FSI for pressure force. */
                pressure_force_on_shell.exec();

                kernel_summation.exec();

                emitter_buffer_inflow_condition.exec();  

                // windkessel model implementation
                if (physical_time >= updateP_n * accumulated_time)
                {
                    outflow_pressure_condition1.getTargetPressure()->updateNextPressure();
                    outflow_pressure_condition2.getTargetPressure()->updateNextPressure();
                    outflow_pressure_condition3.getTargetPressure()->updateNextPressure();
                    outflow_pressure_condition4.getTargetPressure()->updateNextPressure();
                    outflow_pressure_condition5.getTargetPressure()->updateNextPressure();
                   
                    ++updateP_n;
                }

                outflow_pressure_condition1.exec(dt); 
                outflow_pressure_condition2.exec(dt);
                outflow_pressure_condition3.exec(dt); 
                outflow_pressure_condition4.exec(dt); 
                outflow_pressure_condition5.exec(dt);   

                density_relaxation.exec(dt);

                 Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    dt_s = shell_time_step_size.exec();
                    if (dt - dt_s_sum < dt_s)
                        dt_s = dt - dt_s_sum;
                    shell_stress_relaxation_first.exec(dt_s);

                    constrain_holder.exec(dt_s);
                    shell_velocity_damping.exec(dt_s);
                    shell_rotation_damping.exec(dt_s);
                    constrain_holder.exec(dt_s);

                    shell_stress_relaxation_second.exec(dt_s);
                    dt_s_sum += dt_s;
                }
                average_velocity_and_acceleration.update_averages_.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9)
                          << "N=" << number_of_iterations
                          << "	Time = " << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "	dt_s = " << dt_s << "\n";
            }
            number_of_iterations++;

            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();

            /** Water block configuration and periodic condition. */
            inflow_injection.injection.exec();
            outflow_injection_1.injection.exec();
            outflow_injection_2.injection.exec();
            outflow_injection_3.injection.exec();
            outflow_injection_4.injection.exec();
            outflow_injection_5.injection.exec();

            inflow_injection.deletion.exec();
            outflow_injection_1.deletion.exec();
            outflow_injection_2.deletion.exec();
            outflow_injection_3.deletion.exec();
            outflow_injection_4.deletion.exec();
            outflow_injection_5.deletion.exec();

            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }

            water_block.updateCellLinkedList();
            shell_update_normal.exec();
            shell_body.updateCellLinkedList();
            shell_curvature_inner.updateConfiguration();
            shell_average_curvature.exec();
            shell_water_contact.updateConfiguration();
            water_block_complex.updateConfiguration();

            interval_updating_configuration += TickCount::now() - time_instance;
            boundary_indicator.exec();

            inflow_injection.tag_buffer_particles.exec();
            outflow_injection_1.tag_buffer_particles.exec();
            outflow_injection_2.tag_buffer_particles.exec();
            outflow_injection_3.tag_buffer_particles.exec();
            outflow_injection_4.tag_buffer_particles.exec();
            outflow_injection_5.tag_buffer_particles.exec();
        }
        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();
        compute_inlet_transient_flow_rate.exec();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    
    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
              << interval_computing_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    return 0;
}
