/**
 * @file Abundance.h
 * @author  C.Marsh
 * @date 15/07/2018
 * @section LICENSE
 *
 * @section DESCRIPTION
 *
 * This derived quantity will calculate the amount of Abundance
 * in the partition with a selectivity
 */
#ifndef DERIVEDQUANTITIES_ABUNDANCE_H_
#define DERIVEDQUANTITIES_ABUNDANCE_H_

// headers
#include "DerivedQuantities/DerivedQuantity.h"

#include "Layers/Children/IntLayer.h"

// namespaces
namespace niwa {
class IntLayer;
class Selectivity;
namespace derivedquantities {

// classes
class Abundance : public niwa::DerivedQuantity {
public:
  // methods
  explicit Abundance(Model* model);
  virtual                     ~Abundance() = default;
  void                        PreExecute() override final;  // TODO play with this concept, but might be a bit too computationally demanding
  void                        Execute() override final;
  void                        DoValidate() override final;
  void                        DoBuild() override final;

protected:
  string                      abundance_layer_label_;
  niwa::layers::IntLayer*     abundance_layer_ = nullptr;
  string                      selectivity_label_;
  Selectivity*                selectivity_ = nullptr;
  vector<vector<vector<float>>> cell_offset_for_selectivity_;
  bool                        length_based_selectivity_ = false;
};

} /* namespace derivedquantities */
} /* namespace niwa */
#endif /* DERIVEDQUANTITIES_ABUNDANCE_H_ */
