// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_GRID_PART_HH
#define DUNE_STUFF_GRID_PART_HH

#include <dune/common/static_assert.hh>

#if HAVE_DUNE_FEM
# include <dune/fem/gridpart/levelgridpart.hh>
# include <dune/fem/gridpart/common/capabilities.hh>
#endif // HAVE_DUNE_FEM

namespace Dune {
namespace Stuff {
namespace Grid {

#if HAVE_DUNE_FEM


template< class GridImp >
class LevelPartWithGridsOriginalGridView
  : public Dune::Fem::LevelGridPart< GridImp >
{
public:
  typedef Dune::Fem::LevelGridPart< GridImp > OriginalGridPartType;

  typedef typename OriginalGridPartType::LevelGridView GridViewType;

  explicit LevelPartWithGridsOriginalGridView(GridImp& g)
    : OriginalGridPartType(g)
  {}

  LevelPartWithGridsOriginalGridView(GridImp& g, const int l)
    : OriginalGridPartType(g, l)
  {}

  GridViewType gridView() const
  {
    return GridViewType(OriginalGridPartType::grid(), OriginalGridPartType::level());
  }
}; // class LevelPartWithGridsOriginalGridView


#else // HAVE_DUNE_FEM


template< class GridImp >
class LevelPartWithGridsOriginalGridView
{
  static_assert(AlwaysFalse< GridImp >::value, "You are missing dune-fem!");
}; // class LevelPartWithGridsOriginalGridView


#endif // HAVE_DUNE_FEM

} // namespace Grid
} // namespace Stuff

#if HAVE_DUNE_FEM

namespace Fem {
namespace GridPartCapabilities {


template< class GridType >
struct hasGrid< Stuff::Grid::LevelPartWithGridsOriginalGridView< GridType > >
{
  static const bool v = hasGrid< typename Stuff::Grid::LevelPartWithGridsOriginalGridView< GridType >::OriginalGridPartType >::v;
};

template< class GridType >
struct hasSingleGeometryType< Stuff::Grid::LevelPartWithGridsOriginalGridView< GridType > >
{
  static const bool v = hasSingleGeometryType< typename Stuff::Grid::LevelPartWithGridsOriginalGridView< GridType >::OriginalGridPartType >::v;
  static const unsigned int topologyId = hasSingleGeometryType< typename Stuff::Grid::LevelPartWithGridsOriginalGridView< GridType >::OriginalGridPartType >::topologyId;
};

template< class GridType >
struct isCartesian< Stuff::Grid::LevelPartWithGridsOriginalGridView< GridType > >
{
  static const bool v = isCartesian< typename Stuff::Grid::LevelPartWithGridsOriginalGridView< GridType >::OriginalGridPartType >::v;
};

template< class GridType, int codim  >
struct hasEntity< Stuff::Grid::LevelPartWithGridsOriginalGridView< GridType >, codim >
{
  static const bool v = hasEntity< typename Stuff::Grid::LevelPartWithGridsOriginalGridView< GridType >::OriginalGridPartType, codim >::v;
};

template< class GridType >
struct isParallel< Stuff::Grid::LevelPartWithGridsOriginalGridView< GridType > >
{
  static const bool v = isParallel< typename Stuff::Grid::LevelPartWithGridsOriginalGridView< GridType >::OriginalGridPartType >::v;
};

template< class GridType, int codim >
struct canCommunicate< Stuff::Grid::LevelPartWithGridsOriginalGridView< GridType >, codim >
{
  static const bool v = canCommunicate< typename Stuff::Grid::LevelPartWithGridsOriginalGridView< GridType >::OriginalGridPartType, codim >::v;
};

template< class GridType >
struct isConforming< Stuff::Grid::LevelPartWithGridsOriginalGridView< GridType > >
{
  static const bool v = isConforming< typename Stuff::Grid::LevelPartWithGridsOriginalGridView< GridType >::OriginalGridPartType >::v;
};


} // namespace GridPartCapabilities
} // namespace Fem

#endif // HAVE_DUNE_FEM

} // namespace Dune

#endif // DUNE_STUFF_GRID_PART_HH
