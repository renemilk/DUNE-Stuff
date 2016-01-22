// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_FUNCTIONS_SPE10MODEL2_HH
#define DUNE_STUFF_FUNCTIONS_SPE10MODEL2_HH

#include <iostream>
#include <memory>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/fvector.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/functions/global.hh>

namespace Dune {
namespace Stuff {
namespace Functions {
namespace Spe10 {
namespace internal {


static const std::string model2_filename = "spe_perm.dat";
static const size_t model2_x_elements = 60;
static const size_t model2_y_elements = 220;
static const size_t model2_z_elements = 85;
static const double model_2_length_x = 365.76;
static const double model_2_length_y = 670.56;
static const double model_2_length_z = 51.816;


} // namespace internal


// required for FunctionsProvider
template< class E, class D, size_t d, class R, size_t r, size_t rC = 1 >
class Model2
  : public LocalizableFunctionInterface< E, D, d, R, r, rC >
{
  Model2() { static_assert(AlwaysFalse< E >::value, "Not available for these dimensions!"); }
};


/**
 * Grid originally had LL (0,0,0) to UR (365.76, 670.56, 51.816) corners
 *
 */
template <class EntityImp, class DomainFieldImp, class RangeFieldImp>
class Model2<EntityImp, DomainFieldImp, 3, RangeFieldImp, 3, 3>
  : public Stuff::GlobalFunctionInterface<EntityImp, DomainFieldImp, 3, RangeFieldImp, 3, 3> {
  typedef Model2<EntityImp, DomainFieldImp, 3, RangeFieldImp, 3, 3>                         ThisType;
  typedef Stuff::GlobalFunctionInterface<EntityImp, DomainFieldImp, 3, RangeFieldImp, 3, 3> BaseType;
public:
  using typename BaseType::DomainType;
  using BaseType::dimDomain;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".spe10.model2";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["filename"] = internal::model2_filename;
    config["upper_right"] = "[" + Common::toString(internal::model_2_length_x) + " "
                                + Common::toString(internal::model_2_length_y) + " "
                                + Common::toString(internal::model_2_length_z) + "]";
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr< ThisType > create(const Common::Configuration config = default_config(),
                                            const std::string sub_name = static_id())
  {
    const Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    return Common::make_unique< ThisType >(
          cfg.get("filename", default_cfg.get< std::string >("filename")),
          cfg.get("upper_right",  default_cfg.get< DomainType >("upper_right")));
  } // ... create(...)

  Model2(const std::string& data_filename = default_config().get< std::string >("filename"),
         const DomainType& upper_right = default_config().get< DomainType >("upper_right"))
    : num_elements_{{internal::model2_x_elements,
                     internal::model2_y_elements,
                     internal::model2_z_elements}}
    , deltas_{{upper_right[0]/internal::model2_x_elements,
               upper_right[1]/internal::model2_y_elements,
               upper_right[2]/internal::model2_z_elements}}
    , permeability_(nullptr)
    , permMatrix_(0.0)
    , filename_(data_filename)
  {
    readPermeability();
  }

  virtual ~Model2() {
    delete permeability_;
    permeability_ = nullptr;
  }

  virtual std::string name() const override final
  {
    return static_id();
  }

  //! currently used in gdt assembler
  virtual void evaluate(const DomainType& x, typename BaseType::RangeType& diffusion) const final override
  {
    if (!permeability_) {
      DSC_LOG_ERROR_0 << "The SPE10-permeability data file could not be opened. This file does\n"
                      << "not come with the dune-multiscale repository due to file size. To download it\n"
                      << "execute\n"
                      << "wget http://www.spe.org/web/csp/datasets/por_perm_case2a.zip\n"
                      << "unzip the file and move the file 'spe_perm.dat' to\n"
                      << "dune-multiscale/dune/multiscale/problems/spe10_permeability.dat!\n";
      DUNE_THROW(IOError, "Data file for Groundwaterflow permeability could not be opened!");
    }

    // 3 is the maximum space dimension
    for (size_t dim = 0; dim < dimDomain; ++dim)
      permIntervalls_[dim] = std::min(size_t(std::floor(x[dim] / deltas_[dim])), num_elements_[dim] -1);

    const int offset = permIntervalls_[0] + permIntervalls_[1] * internal::model2_x_elements
                       + permIntervalls_[2] * internal::model2_y_elements * internal::model2_x_elements;
    for (size_t dim = 0; dim < dimDomain; ++dim) {
      const auto idx = offset + dim * (internal::model2_x_elements*internal::model2_z_elements*internal::model2_z_elements);
      diffusion[dim][dim] = permeability_[idx];
    }
  }

  virtual size_t order() const
  {
    return 0u;
  }

private:
  void readPermeability()
  {
    std::ifstream file(filename_);
    double val;
    if (!file) { // file couldn't be opened
      return;
    }
    file >> val;
    int counter = 0;
    permeability_ = new double[dimDomain*internal::model2_x_elements*internal::model2_y_elements*internal::model2_z_elements];
    while (!file.eof()) {
      // keep reading until end-of-file
      permeability_[counter++] = val;
      file >> val; // sets EOF flag if no value found
    }
    file.close();
  }

  std::array<size_t, dimDomain> num_elements_;
  std::array<double, dimDomain> deltas_;
  double* permeability_; //! TODO automatic memory
  mutable typename BaseType::DomainType permIntervalls_;
  mutable Dune::FieldMatrix<double, BaseType::DomainType::dimension, BaseType::DomainType::dimension> permMatrix_;
  const std::string filename_;
};


} // namespace Spe10
} // namespace Functions
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTIONS_SPE10MODEL2_HH

