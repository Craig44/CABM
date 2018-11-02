/**
 * @file Biomass.h
 * @author  C.Marsh
 * @date 15/07/2018
 * @section LICENSE
 *
 * @section DESCRIPTION
 *
 * This derived quantity will calculate the amount of biomass
 * in the partition with a selectivity
 */
#ifndef DERIVEDQUANTITIES_BIOMASS_H_
#define DERIVEDQUANTITIES_BIOMASS_H_

// headers
#include "DerivedQuantities/DerivedQuantity.h"

#include "Layers/Children/IntLayer.h"

// namespaces
namespace niwa {
class IntLayer;
class Selectivity;
namespace derivedquantities {

// classes
class Biomass : public niwa::DerivedQuantity {
public:
  // methods
  explicit Biomass(Model* model);
  virtual                     ~Biomass() = default;
  void                        PreExecute() override final;  // TODO play with this concept, but might be a bit too computationally demanding
  void                        Execute() override final;
  void                        DoValidate() override final;
  void                        DoBuild() override final;

protected:
  string                      biomass_layer_label_;
  niwa::layers::IntLayer*     biomass_layer_ = nullptr;
  vector<string>                      selectivity_label_;
  vector<Selectivity*>                selectivity_;
  bool                        length_based_selectivity_ = false;
};

} /* namespace derivedquantities */
} /* namespace niwa */
#endif /* DERIVEDQUANTITIES_MATURE_BIOMASS_H_ */
