﻿#ifndef DUNE_STUFF_FUNCTION_INTERFACE_HH
#define DUNE_STUFF_FUNCTION_INTERFACE_HH

#include <vector>
#include <memory>
#include <string>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/version.hh>
#include <dune/common/deprecated.hh>

#include <dune/geometry/genericreferenceelements.hh>

#if HAVE_DUNE_GRID
# include <dune/grid/io/file/vtk/vtkwriter.hh>
#endif

#if HAVE_DUNE_FEM
# include <dune/fem/function/common/function.hh>
# include <dune/fem/space/common/functionspace.hh>
#endif

namespace Dune {
namespace Stuff {
#if HAVE_DUNE_GRID
namespace Function {

// forward, include is below
template< class GridViewType, int dimRange >
class VisualizationAdapter;

}
#endif // HAVE_DUNE_GRID


/**
 *  \brief  Interface for a set of matrix valued globalvalued functions, which can be evaluated locally on one Entity.
 *
 *  \note   see specialization for rangeDimCols = 1 for vector and scalar valued localfunction sets.
 */
template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class LocalfunctionSetInterface
{
public:
  typedef EntityImp EntityType;

  typedef DomainFieldImp                                  DomainFieldType;
  static const unsigned int                               dimDomain = domainDim;
  typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;

  typedef RangeFieldImp                                               RangeFieldType;
  static const unsigned int                                           dimRange = rangeDim;
  static const unsigned int                                           dimRangeCols = rangeDimCols;
  typedef Dune::FieldMatrix< RangeFieldType, dimRange, dimRangeCols > RangeType;

  typedef std::string JacobianRangeType; // <- this is yet unclear, but we need a type

  LocalfunctionSetInterface(const EntityType& ent)
    : entity_(ent)
  {}

  virtual ~LocalfunctionSetInterface() {}

  virtual const EntityType& entity() const
  {
    return entity_;
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented.''
   * @{
   **/
  virtual size_t size() const = 0;

  virtual size_t order() const = 0;

  virtual void evaluate(const DomainType& /*xx*/, std::vector< RangeType >& /*ret*/) const = 0;

//  virtual void jacobian(const DomainType& /*xx*/, std::vector< JacobianRangeType >& /*ret*/) const = 0;
  /* @} */

  /**
   * \defgroup provided ´´These methods are provided by the interface.''
   * @{
   **/
  std::vector< RangeType > evaluate(const DomainType& xx) const
  {
    std::vector< RangeType > ret(size(), RangeType(0));
    evaluate(xx, ret);
    return ret;
  }

//  std::vector< JacobianRangeType > jacobian(const DomainType& xx) const
//  {
//    std::vector< JacobianRangeType > ret(size(), JacobianRangeType(0));
//    jacobian(xx, ret);
//    return ret;
//  }
  /* @} */

protected:
  bool is_a_valid_point(const DomainType& xx) const
  {
    const auto& reference_element = GenericReferenceElements< DomainFieldType, dimDomain >::general(entity().type());
    return reference_element.checkInside(xx);
  }

  const EntityType& entity_;
}; // class LocalfunctionSetInterface


/**
 *  \brief  Interface for a set of scalar or vector globalvalued functions, which can be evaluated locally on one Entity.
 */
template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class LocalfunctionSetInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
{
public:
  typedef EntityImp EntityType;

  typedef DomainFieldImp                                  DomainFieldType;
  static const unsigned int                               dimDomain = domainDim;
  typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;

  typedef RangeFieldImp                                 RangeFieldType;
  static const unsigned int                             dimRange = rangeDim;
  static const unsigned int                             dimRangeCols = 1;
  typedef Dune::FieldVector< RangeFieldType, dimRange > RangeType;

  typedef Dune::FieldMatrix< RangeFieldType, dimRange, dimDomain > JacobianRangeType;

  LocalfunctionSetInterface(const EntityType& ent)
    : entity_(ent)
  {}

  virtual ~LocalfunctionSetInterface() {}

  const EntityType& entity() const
  {
    return entity_;
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented.''
   * @{
   **/
  virtual size_t size() const = 0;

  virtual size_t order() const = 0;

  virtual void evaluate(const DomainType& /*xx*/, std::vector< RangeType >& /*ret*/) const = 0;

  virtual void jacobian(const DomainType& /*xx*/, std::vector< JacobianRangeType >& /*ret*/) const = 0;
  /* @} */

  /**
   * \defgroup provided ´´These methods are provided by the interface.''
   * @{
   **/
  std::vector< RangeType > evaluate(const DomainType& xx) const
  {
    std::vector< RangeType > ret(size(), RangeType(0));
    evaluate(xx, ret);
    return ret;
  }

  std::vector< JacobianRangeType > jacobian(const DomainType& xx) const
  {
    std::vector< JacobianRangeType > ret(size(), JacobianRangeType(0));
    jacobian(xx, ret);
    return ret;
  }
  /* @} */

protected:
  bool is_a_valid_point(const DomainType& xx) const
  {
    const auto& reference_element = GenericReferenceElements< DomainFieldType, dimDomain >::general(entity().type());
    return reference_element.checkInside(xx);
  }

  const EntityType& entity_;
}; // class LocalfunctionSetInterface< ...., 1 >


/**
 *  \brief  Interface for matrix valued globalvalued functions, which can be evaluated locally on one Entity.
 *
 *          This is the interface for matrixvalued functions, see the specialization for rangeDimCols = 1 for scalar
 *          and vector valued functions.
 */
template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class LocalfunctionInterface
  : public LocalfunctionSetInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols >
{
  typedef LocalfunctionSetInterface
      < EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols > BaseType;
public:
  typedef EntityImp EntityType;

  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const unsigned int                   dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType       DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int                 dimRange = BaseType::dimRange;
  static const unsigned int                 dimRangeCols = BaseType::dimRangeCols;
  typedef typename BaseType::RangeType      RangeType;

  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  LocalfunctionInterface(const EntityType& ent)
    : BaseType(ent)
  {}

  virtual ~LocalfunctionInterface() {}

  /**
   * \defgroup haveto ´´These methods have to be implemented in addition to the ones required from the BaseType.''
   * @{
   **/
  virtual void evaluate(const DomainType& /*xx*/, RangeType& /*ret*/) const = 0;

  virtual void jacobian(const DomainType& /*xx*/, JacobianRangeType& /*ret*/) const = 0;
  /* @} */

  /**
   * \defgroup providedbase ´´These methods are provided by the interface to please LocalfunctionSetInterface.''
   * @{
   **/
  virtual size_t size() const
  {
    return 1;
  }

  virtual void evaluate(const DomainType& xx, std::vector< RangeType >& ret) const
  {
    assert(ret.size() >= 1);
    evaluate(xx, ret[0]);
  }

  virtual void jacobian(const DomainType& xx, std::vector< JacobianRangeType >& ret) const
  {
    assert(ret.size() >= 1);
    jacobian(xx, ret[0]);
  }
  /* @} */

  /**
   * \defgroup provided ´´These methods are provided by the interface.''
   * @{
   **/
  RangeType evaluate(const DomainType& xx) const
  {
    RangeType ret(0);
    evaluate(xx, ret);
    return ret;
  }

  JacobianRangeType jacobian(const DomainType& xx) const
  {
    JacobianRangeType ret(0);
    jacobian(xx, ret);
    return ret;
  }
  /* @} */
}; // class LocalfunctionInterface


class IsLocalizableFunction {};


/**
 * \brief Interface for functions which provide a LocalfunctionInterface for an entity.
 */
template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class LocalizableFunctionInterface
  : public IsLocalizableFunction
{
  typedef LocalizableFunctionInterface
      < EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols > ThisType;
public:
  typedef EntityImp EntityType;

  typedef DomainFieldImp    DomainFieldType;
  static const unsigned int dimDomain = domainDim;

  typedef RangeFieldImp     RangeFieldType;
  static const unsigned int dimRange = rangeDim;
  static const unsigned int dimRangeCols = rangeDimCols;

  typedef LocalfunctionInterface
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols > LocalfunctionType;

  typedef typename LocalfunctionType::DomainType        DomainType;
  typedef typename LocalfunctionType::RangeType         RangeType;
  typedef typename LocalfunctionType::JacobianRangeType JacobianRangeType;

  virtual ~LocalizableFunctionInterface() {}

  static std::string static_id()
  {
    return "dune.stuff.function";
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented.''
   * @{
   **/
  virtual std::unique_ptr< LocalfunctionType > local_function(const EntityType& /*entity*/) const = 0;

  virtual ThisType* copy() const = 0;
  /* @} */

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string name() const
  {
    return static_id();
  }
  /* @} */

#if HAVE_DUNE_GRID
  template< class GridViewType >
  void visualize(const GridViewType& grid_view, const std::string filename) const
  {
    if (filename.empty()) DUNE_THROW(RangeError, "Empty filename given!");
    auto adapter = std::make_shared< Stuff::Function::VisualizationAdapter< GridViewType, dimRange > >(*this);
    VTKWriter< GridViewType > vtk_writer(grid_view, VTK::nonconforming);
    vtk_writer.addVertexData(adapter);
    vtk_writer.write(filename);
  } // ... visualize(...)
#endif // HAVE_DUNE_GRID
}; // class LocalizableFunctionInterface


/**
 * \brief Interface for scalar and vector valued stationary function.
 */
template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class DUNE_DEPRECATED_MSG("Please derive your functions from LocalizableFunctionInterface in the future!") FunctionInterface
#if HAVE_DUNE_FEM
  : public Dune::Fem::Function< Dune::Fem::FunctionSpace< DomainFieldImp, RangeFieldImp, domainDim, rangeDim >,
                                FunctionInterface< DomainFieldImp, domainDim, RangeFieldImp, rangeDim > >
#endif // HAVE_DUNE_FEM
{
public:
  typedef DomainFieldImp                                  DomainFieldType;
  static const unsigned int                               dimDomain = domainDim;
  typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;

  typedef RangeFieldImp                                 RangeFieldType;
  static const unsigned int                             dimRange = rangeDim;
  typedef Dune::FieldVector< RangeFieldType, dimRange > RangeType;
#if HAVE_DUNE_FEM
  typedef typename Dune::Fem::Function< Dune::Fem::FunctionSpace< DomainFieldImp, RangeFieldImp, domainDim, rangeDim >,
                               FunctionInterface< DomainFieldImp, domainDim, RangeFieldImp, rangeDim > >
        ::JacobianRangeType JacobianRangeType;
#else
  typedef Dune::FieldMatrix< RangeFieldType, dimRange, dimDomain > JacobianRangeType;
#endif

  virtual ~FunctionInterface() {}

  static std::string static_id()
  {
    return "dune.stuff.function";
  }

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string name() const
  {
    return static_id();
  }

  virtual int order() const
  {
    return -1;
  }
  /* @} */

  /** \defgroup must This method has to be implemented.'' */
  /* @{ */
  virtual void evaluate(const DomainType& /*x*/, RangeType& /*ret*/) const = 0;
  /* @} */

  virtual RangeType evaluate(const DomainType& x) const
  {
    RangeType ret;
    evaluate(x, ret);
    return ret;
  }

  virtual void jacobian(const DomainType& /*x*/, JacobianRangeType& /*ret*/) const
  {
    DUNE_THROW(Dune::NotImplemented, "You really have to implement this!");
  }

  virtual JacobianRangeType jacobian(const DomainType& x) const
  {
    JacobianRangeType ret;
    jacobian(x, ret);
    return ret;
  }
}; // class FunctionInterface


/**
 * \brief Interface for scalar and vector valued timedependent functions.
 */
template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class TimedependentFunctionInterface
{
public:
  typedef DomainFieldImp                                  DomainFieldType;
  static const unsigned int                               dimDomain = domainDim;
  typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;

  typedef RangeFieldImp                                 RangeFieldType;
  static const unsigned int                             dimRange = rangeDim;
  typedef Dune::FieldVector< RangeFieldType, dimRange > RangeType;

  virtual ~TimedependentFunctionInterface() {}

  static std::string static_id()
  {
    return "dune.stuff.timedependentfunction";
  }

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string name() const
  {
    return static_id();
  }

  virtual int order() const
  {
    return -1;
  }
  /* @} */

  /** \defgroup must ´´This method has to be implemented.'' */
  /* @{ */
  virtual void evaluate(const DomainType& /*xx*/, const double& /*tt*/, RangeType& /*ret*/) const = 0;
  /* @} */
}; // class TimedependentFunctionInterface


//! use this to throw a stationary function into an algorithm that expects an instationary one
template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows >
struct TimeFunctionAdapter
  : public Dune::Stuff::TimedependentFunctionInterface< DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows >
{
  typedef Dune::Stuff::FunctionInterface< DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows > WrappedType;

  TimeFunctionAdapter(const WrappedType& wr)
    : wrapped_(wr)
  {}

  virtual void evaluate(const typename WrappedType::DomainType& x,
                        typename WrappedType::RangeType& ret) const
  {
    wrapped_(x, ret);
  }

  virtual void evaluate(const typename WrappedType::DomainType& x,
                        const double& /*t*/,
                        typename WrappedType::RangeType& ret) const
  {
    wrapped_(x, ret);
  }

  const WrappedType& wrapped_;
};

template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows >
TimeFunctionAdapter< DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows >
timefunctionAdapted(const Dune::Stuff::FunctionInterface< DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows >& wrapped)
{
  return TimeFunctionAdapter< DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows >(wrapped);
}


} // namespace Stuff
} // namespace Dune

#include "default.hh"

#endif // DUNE_STUFF_FUNCTION_INTERFACE_HH
