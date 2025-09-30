/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file mesh_iterators.h
 * @brief This is for the base functions for mesh iterator.
 * There are two types of functions: one is for static ranged
 * which are defined by template parameters,
 * the other for dynamics ranges which are input parameters.
 * @author  Xiangyu Hu
 */

#ifndef MESH_ITERATORS_H
#define MESH_ITERATORS_H

#include "base_data_package.h"

#include "execution_policy.h"

namespace SPH
{
/** iteration with void (non_value_returning) function. 2D case. */
template <int lower0, int upper0,
          int lower1, int upper1, typename FunctionOnEach>
inline void mesh_for_each2d(const FunctionOnEach &function);
template <int lower, int upper, typename FunctionOnEach>
inline void mesh_for_each2d(const FunctionOnEach &function)
{
    mesh_for_each2d<lower, upper, lower, upper, FunctionOnEach>(function);
};

template <typename FunctionOnEach>
inline void mesh_for_each_neighbor2d(int depth, const FunctionOnEach &function);
template <typename FunctionOnEach>
inline void mesh_for_each_neighbor3d(int depth, const FunctionOnEach &function);

/** iteration with boolean return function. 2D case. */
template <int lower0, int upper0,
          int lower1, int upper1, typename CheckOnEach>
inline Array2i mesh_find_if2d(const CheckOnEach &function);
template <int lower, int upper, typename CheckOnEach>
inline Array2i mesh_find_if2d(const CheckOnEach &function)
{
    return mesh_find_if2d<lower, upper, lower, upper, CheckOnEach>(function);
};
template <int lower, int upper, typename CheckOnEach>
inline bool mesh_any_of2d(const CheckOnEach &function)
{
    return (mesh_find_if2d<lower, upper, lower, upper, CheckOnEach>(
                function) != Array2i(upper, upper))
        .all();
};

/** iteration with void (non_value_returning) function. 3D case. */
template <int lower0, int upper0,
          int lower1, int upper1,
          int lower2, int upper2, typename FunctionOnEach>
inline void mesh_for_each3d(const FunctionOnEach &function);
template <int lower, int upper, typename FunctionOnEach>
inline void mesh_for_each3d(const FunctionOnEach &function)
{
    mesh_for_each3d<lower, upper, lower, upper, lower, upper, FunctionOnEach>(function);
};

/** iteration with boolean return function.  3D case. */
template <int lower0, int upper0,
          int lower1, int upper1,
          int lower2, int upper2, typename CheckOnEach>
inline Array3i mesh_find_if3d(const CheckOnEach &function);
template <int lower, int upper, typename CheckOnEach>
inline Array3i mesh_find_if3d(const CheckOnEach &function)
{
    return mesh_find_if3d<lower, upper, lower, upper, lower, upper, CheckOnEach>(function);
};
template <int lower, int upper, typename CheckOnEach>
inline bool mesh_any_of3d(const CheckOnEach &function)
{
    return (mesh_find_if3d<lower, upper, lower, upper, lower, upper, CheckOnEach>(
                function) != Array3i(upper, upper, upper))
        .all();
};

template <typename FunctionOnEach>
void mesh_for_each(const Arrayi &lower, const Arrayi &upper, const FunctionOnEach &function);
template <typename FunctionOnEach>
void mesh_for_column_major(const Arrayi &lower, const Arrayi &upper, const FunctionOnEach &function);
template <typename FunctionOnEach>
Arrayi mesh_find_if(const Arrayi &lower, const Arrayi &upper, const FunctionOnEach &function);
template <typename FunctionOnEach>
bool mesh_any_of(const Arrayi &lower, const Arrayi &upper, const FunctionOnEach &function)
{
    return mesh_find_if(lower, upper, function).matrix() != upper.matrix();
};

using MeshRange = std::pair<Arrayi, Arrayi>;
/** Iterator on the mesh by looping index. sequential computing. */
template <typename LocalFunction, typename... Args>
void mesh_for(const MeshRange &mesh_range, const LocalFunction &local_function, Args &&...args);
/** Iterator on the mesh by looping index. parallel computing. */
template <typename LocalFunction, typename... Args>
void mesh_parallel_for(const MeshRange &mesh_range, const LocalFunction &local_function, Args &&...args);

/** Iterator on the mesh by looping index. parallel computing. */
template <typename LocalFunction, typename... Args>
void mesh_for(const execution::SequencedPolicy &seq, const MeshRange &mesh_range,
              const LocalFunction &local_function, Args &&...args)
{
    mesh_for(mesh_range, local_function, std::forward<Args>(args)...);
};

template <typename LocalFunction, typename... Args>
void mesh_for(const execution::ParallelPolicy &par_host, const MeshRange &mesh_range,
              const LocalFunction &local_function, Args &&...args)
{
    mesh_parallel_for(mesh_range, local_function, std::forward<Args>(args)...);
};

template <typename FunctionOnData>
void package_for(const execution::SequencedPolicy &seq, UnsignedInt start_index,
                 UnsignedInt num_grid_pkgs, const FunctionOnData &function)
{
    for (size_t i = start_index; i != num_grid_pkgs; ++i)
        function(i);
}

template <typename FunctionOnData>
void package_for(const execution::ParallelPolicy &par_host, UnsignedInt start_index,
                 UnsignedInt num_grid_pkgs, const FunctionOnData &function)
{
    parallel_for(IndexRange(start_index, num_grid_pkgs), [&](const IndexRange &r)
                 {
                    for (size_t i = r.begin(); i != r.end(); ++i)
                    {
                        function(i);
                    } }, ap);
}

template <typename FunctionOnData>
void package_for(const execution::ParallelDevicePolicy &par_device,
                 UnsignedInt start_index, UnsignedInt num_grid_pkgs,
                 const FunctionOnData &function);
} // namespace SPH
#endif // MESH_ITERATORS_H