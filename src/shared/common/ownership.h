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
 * @file 	ownership.h
 * @brief 	Here, the classes related to the ownership of objects are defined
 * 			by implementing Resource Acquisition Is Initialization (RAII).
 *          The basic idea in SPHinXsys to use defined ownership.
 * 			Unique pointers are used as much as possible and shared pointer only used
 * 			when a default smart pointer is generated as argument in the constructors.
 *          In this case, for objects with statically defined lifetime, references will be used.
 *          For those with dynamical lifetime, we will use smart-pointer keepers to handel.
 *          Basically, new a raw pointer is forbidden through out the code except for
 * 			raw matrixes, such as the cell linked list mesh, that will be deleted in destructor,
 *          especially, when new features are introduced.
 *          Generally, the classes here are private member and
 * 			the smart pointer is used only when the memory is allocated.
 *          After that, we take reference or raw pointer out and use them as observer.
 * @author	Chi Zhang and Xiangyu Hu
 */
#ifndef OWNERSHIP_H
#define OWNERSHIP_H

#include "base_data_type.h"

#include <string_view>

template <typename T>
constexpr std::string_view type_name()
{
#if defined(__clang__)
    constexpr std::string_view p = __PRETTY_FUNCTION__;
    constexpr std::string_view prefix = "std::string_view type_name() [T = ";
    constexpr std::string_view suffix = "]";
#elif defined(__GNUC__)
    constexpr std::string_view p = __PRETTY_FUNCTION__;
    constexpr std::string_view prefix = "constexpr std::string_view type_name() [with T = ";
    constexpr std::string_view suffix = "; std::string_view = std::basic_string_view<char>]";
#elif defined(_MSC_VER)
    constexpr std::string_view p = __FUNCSIG__;
    constexpr std::string_view prefix = "class std::basic_string_view<char,struct std::char_traits<char> > __cdecl type_name<";
    constexpr std::string_view suffix = ">(void)";
#else
#error "Unsupported compiler"
#endif
    const std::size_t start = p.find(prefix) + prefix.size();
    const std::size_t end = p.rfind(suffix);
    return p.substr(start, end - start);
}

namespace SPH
{
template <class CastingType, class OwnerType, class CastedType>
CastingType *DynamicCast(OwnerType *owner, CastedType *casted)
{
    CastingType *tmp = dynamic_cast<CastingType *>(casted);
    if (tmp == nullptr)
    {
        std::cout << "\n Error: pointer DynamicCasting " << type_name<CastedType>() << " leads to nullptr! \n";
        std::cout << "\n This error locates in " << type_name<OwnerType>() << '\n';
        exit(1);
    }
    return tmp;
}

template <class CastingType, class OwnerType, class CastedType>
CastingType &DynamicCast(OwnerType *owner, CastedType &casted)
{
    CastingType *tmp = dynamic_cast<CastingType *>(&casted);
    if (tmp == nullptr)
    {
        std::cout << "\n Error: reference DynamicCasting " << type_name<CastedType>() << " leads to nullptr! \n";
        std::cout << "\n This error locates in " << type_name<OwnerType>() << '\n';
        exit(1);
    }
    return *tmp;
}

template <class T>
using UniquePtr = std::unique_ptr<T>;

template <class T, typename... Args>
UniquePtr<T> makeUnique(Args &&...args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

/**
 * @class UniquePtrKeeper
 * @brief A wrapper to provide an ownership for a new derived object
 * which previous often generated by new a raw pointer.
 */
template <class BaseType>
class UniquePtrKeeper
{
  public:
    /** output the observer as derived pointer */
    template <class DerivedType, typename... Args>
    DerivedType *createPtr(Args &&...args)
    {
        ptr_member_.reset(new DerivedType(std::forward<Args>(args)...));
        return static_cast<DerivedType *>(ptr_member_.get());
    };

    /** output the observer as derived reference */
    template <class DerivedType, typename... Args>
    DerivedType &createRef(Args &&...args)
    {
        ptr_member_.reset(new DerivedType(std::forward<Args>(args)...));
        return *static_cast<DerivedType *>(ptr_member_.get());
    };

    /** output the observer as pointer */
    BaseType *movePtr(UniquePtr<BaseType> moved_unique_ptr)
    {
        ptr_member_ = std::move(moved_unique_ptr);
        return ptr_member_.get();
    };

    BaseType *getPtr()
    {
        return ptr_member_.get();
    };

  private:
    UniquePtr<BaseType> ptr_member_;
};

/**
 * @class UniquePtrsKeeper
 * @brief A wrapper to provide an ownership for
 * a vector of base class pointers which point to derived objects.
 * It should be a private member.
 */
template <class BaseType>
class UniquePtrsKeeper
{
  public:
    /** used to create a new derived object in the vector
     * and output its pointer as observer */
    template <class DerivedType, typename... Args>
    DerivedType *createPtr(Args &&...args)
    {
        ptr_keepers_.push_back(UniquePtrKeeper<BaseType>());
        BaseType *observer =
            ptr_keepers_.back().template createPtr<DerivedType>(std::forward<Args>(args)...);
        return static_cast<DerivedType *>(observer);
    };

    UniquePtrKeeper<BaseType> &operator[](size_t index)
    {
        if (index < ptr_keepers_.size())
        {
            return ptr_keepers_[index];
        }
        std::cout << "\n Error in UniquePtrsKeeper : UniquePtr index is out of bound! \n";
        exit(1);
    }

    size_t size() const
    {
        return ptr_keepers_.size();
    }

  private:
    std::vector<UniquePtrKeeper<BaseType>> ptr_keepers_;
};

template <class T>
using SharedPtr = std::shared_ptr<T>;

template <class T, typename... Args>
SharedPtr<T> makeShared(Args &&...args)
{
    return std::make_shared<T>(std::forward<Args>(args)...);
}

/**
 * @class SharedPtrKeeper
 * @brief A wrapper to provide an shared ownership for a new derived object
 * which previous often generated by new a raw pointer.
 */
template <class BaseType>
class SharedPtrKeeper
{
  public:
    /** output the observer as pointer */
    template <class DerivedType, typename... Args>
    DerivedType *resetPtr(Args &&...args)
    {
        ptr_member_ = makeShared<DerivedType>(std::forward<Args>(args)...);
        return static_cast<DerivedType *>(ptr_member_.get());
    };

    /** output the observer as derived reference */
    template <class DerivedType, typename... Args>
    DerivedType &resetRef(Args &&...args)
    {
        return *resetPtr<DerivedType>(std::forward<Args>(args)...);
    };

    /** output the observer as pointer */
    BaseType *assignPtr(SharedPtr<BaseType> shared_ptr)
    {
        ptr_member_ = shared_ptr;
        return ptr_member_.get();
    };

    /** output the observer as reference */
    BaseType &assignRef(SharedPtr<BaseType> shared_ptr)
    {
        ptr_member_ = shared_ptr;
        return *ptr_member_.get();
    };

  private:
    SharedPtr<BaseType> ptr_member_;
};

/**
 * @class SharedPtrsKeeper
 * @brief A wrapper to provide an ownership for
 * a vector of base class shared pointers which point to derived objects.
 * This is designed to be a private member of class whose objects are copyable.
 */
template <class BaseType>
class SharedPtrsKeeper
{
  public:
    /** used to create a new derived object in the vector
     * and output its pointer as observer */
    template <class DerivedType, typename... Args>
    DerivedType *createPtr(Args &&...args)
    {
        ptr_keepers_.push_back(SharedPtrKeeper<BaseType>());
        BaseType *observer =
            ptr_keepers_.back().template resetPtr<DerivedType>(std::forward<Args>(args)...);
        return static_cast<DerivedType *>(observer);
    };

    SharedPtrKeeper<BaseType> &operator[](size_t index)
    {
        if (index < ptr_keepers_.size())
        {
            return ptr_keepers_[index];
        }
        std::cout << "\n Error in UniquePtrsKeeper : UniquePtr index is out of bound! \n";
        exit(1);
    }

    size_t size() const
    {
        return ptr_keepers_.size();
    }

  private:
    std::vector<SharedPtrKeeper<BaseType>> ptr_keepers_;
};

} // namespace SPH
#endif // OWNERSHIP_H
