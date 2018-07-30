/**
 * @file BiomassByCell.h
 * @author  C.Marsh
 * @date 30/07/2018
 * @section LICENSE
 * @description
 *
 */
#ifndef DERIVEDQUANTITIES_BIOMASS_BY_CELL_H_
#define DERIVEDQUANTITIES_BIOMASS_BY_CELL_H_

// headers
#include "DerivedQuantities/DerivedQuantity.h"


namespace niwa {
class Selectivity;
namespace derivedquantities {

//**********************************************************************
//
//
//**********************************************************************
class BiomassByCell : public DerivedQuantity {
public:
  // Functions
  BiomassByCell(Model* model);
  virtual                     ~BiomassByCell() {};
  void                        PreExecute() override final;  // TODO play with this concept, but might be a bit too computationally demanding
  void                        Execute() override final;
  void                        DoValidate() override final;
  void                        DoBuild() override final;

protected:
  // Variables
  string                      selectivity_label_;
  Selectivity*                selectivity_ = nullptr;
  vector<vector<vector<float>>> cell_offset_for_selectivity_;
  bool                        length_based_selectivity_ = false;
  vector<vector<float>>       cache_in_space_;
  vector<vector<float>>       value_in_space_;


};
} /* namespace derivedquantities */
} /* namespace niwa */
#endif /* DERIVEDQUANTITIES_BIOMASS_BY_CELL_H_ */
