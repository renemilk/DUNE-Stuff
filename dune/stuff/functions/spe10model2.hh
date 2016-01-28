// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_FUNCTIONS_SPE10_MODEL2_HH
#define DUNE_STUFF_FUNCTIONS_SPE10_MODEL2_HH

#include "checkerboard.hh"


namespace Dune {
namespace Stuff {
namespace Exceptions {


class spe10_model2_data_file_missing : public Dune::IOError {};


} // namespace Exceptions
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
static const double model_2_min_value = 6.65e-8; // isotropic: 0.000665
static const double model_2_max_value = 20000;


} // namespace internal


// default, to allow for specialization
template< class E, class D, size_t d, class R, size_t r, size_t rC = 1 >
class Model2
  : public LocalizableFunctionInterface< E, D, d, R, r, rC >
{
  Model2() { static_assert(AlwaysFalse< E >::value, "Not available for these dimensions!"); }
};


template< class EntityImp, class DomainFieldImp, class RangeFieldImp >
class Model2< EntityImp, DomainFieldImp, 3, RangeFieldImp, 3, 3 >
  : public Checkerboard< EntityImp, DomainFieldImp, 3, RangeFieldImp, 3, 3 >
{
  typedef Checkerboard< EntityImp, DomainFieldImp, 3, RangeFieldImp, 3, 3 > BaseType;
  typedef Model2< EntityImp, DomainFieldImp, 3, RangeFieldImp, 3, 3 >       ThisType;
public:
  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;
  using typename BaseType::DomainType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;

  static const bool available = true;

  static std::string static_id()
  {
    return LocalizableFunctionInterface
        < EntityImp, DomainFieldImp, 3, RangeFieldImp, 1, 1 >::static_id() + ".spe10.model2";
  } // ... static_id(...)
public:
  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["filename"] = internal::model2_filename;
    config["name"] = static_id();
    config["lower_left"] = "[0 0 0]";
    config["upper_right"] = "[" + Common::toString(internal::model_2_length_x)
                            + " " + Common::toString(internal::model_2_length_y)
                            + " " + Common::toString(internal::model_2_length_z)+ "]";
    config["anisotropic"] = "true";
    config["min"] = Common::toString(internal::model_2_min_value);
    config["max"] = Common::toString(internal::model_2_max_value);
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr< ThisType > create(const Common::Configuration config = BaseType::default_config(),
                                            const std::string sub_name = BaseType::static_id())
  {
    // get correct config
    const Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    // create
    return Common::make_unique< ThisType >(
          cfg.get("filename",     default_cfg.get< std::string >("filename")),
          cfg.get("name",         default_cfg.get< std::string >("name")),
          cfg.get("lower_left",   default_cfg.get< DomainType >("lower_left")),
          cfg.get("upper_right",  default_cfg.get< DomainType >("upper_right")),
          cfg.get("anisotropic",  default_cfg.get< bool >("anisotropic")),
          cfg.get("min",          default_cfg.get< RangeFieldType >("min")),
          cfg.get("max",          default_cfg.get< RangeFieldType >("max")));
  } // ... create(...)

  Model2(const std::string& filename = default_config().get< std::string >("filename"),
         const std::string nm = default_config().get< std::string >("name"),
         const Common::FieldVector< DomainFieldType, dimDomain >& lower_left = default_config().get< DomainType >("lower_left"),
         const Common::FieldVector< DomainFieldType, dimDomain >& upper_right = default_config().get< DomainType >("upper_right"),
         const bool anisotropic = true,
         const RangeFieldType min = default_config().get< RangeFieldType >("min"),
         const RangeFieldType max = default_config().get< RangeFieldType >("max"))
    : BaseType(lower_left,
               upper_right,
               {internal::model2_x_elements, internal::model2_y_elements, internal::model2_z_elements},
               read_values_from_file(filename, anisotropic, min, max),
               nm)
  {}

  virtual std::string type() const override final
  {
    return LocalizableFunctionInterface
        < EntityImp, DomainFieldImp, 3, RangeFieldImp, 1, 1 >::static_id() + ".spe10.model2";
  }

private:
  static std::vector< RangeType > read_values_from_file(const std::string& filename,
                                                        const bool anisotropic,
                                                        const RangeFieldType& min,
                                                        const RangeFieldType& max)

  {
    std::ifstream datafile(filename);
    if (!datafile.is_open())
      DUNE_THROW(Exceptions::spe10_model2_data_file_missing, "could not open '" << filename << "'!");
    if (!(max > min))
      DUNE_THROW(Dune::RangeError,
                 "max (is " << max << ") has to be larger than min (is " << min << ")!");
    const RangeFieldType scale = (max - min) / (internal::model_2_max_value - internal::model_2_min_value);
    const RangeFieldType shift = min - scale*internal::model_2_min_value;
    // read all the data from the file
    const size_t entries_per_dim = internal::model2_x_elements*internal::model2_y_elements*internal::model2_z_elements;
    std::vector< double > data(3*entries_per_dim);
    double tmp = 0;
    size_t counter = 0;
    while (datafile >> tmp && counter < data.size())
      data[counter++] = (tmp*scale) + shift;
    datafile.close();
    if (counter != data.size())
      DUNE_THROW(Dune::IOError,
                 "wrong number of entries in '" << filename << "' (are " << counter << ", should be "
                 << data.size() << ")!");
    // migrate data to appropriate format
    std::vector<RangeType> ret(entries_per_dim, RangeType(0.));
    for (size_t ii = 0; ii < ret.size(); ++ii) {
      ret[ii][0][0] = data[ii];
      ret[ii][1][1] = data[anisotropic ?   entries_per_dim + ii : ii];
      ret[ii][2][2] = data[anisotropic ? 2*entries_per_dim + ii : ii];
    }
    return ret;
  } // ... read_values_from_file(...)
}; // class Model2


} // namespace Spe10
} // namespace Functions
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTIONS_SPE10_MODEL2_HH
