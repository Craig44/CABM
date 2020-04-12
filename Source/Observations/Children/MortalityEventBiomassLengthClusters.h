/**
 * @file MortalityEventBiomassLengthClusters.h
 * @author  C.Marsh github.com/Craig44
 * @date 4/11/2018
 * @section LICENSE
 *
 *
 * @section DESCRIPTION
 *
 * This class takes census information from a mortality process, generates clusters (tow/trip level tags) then randomly selects a sub sample of ages and lengths
 *  to then generate an age frequency either via age length key or direct ageing.
 */
#ifndef OBSERVATIONS_MORTALITY_EVENT_LENGTH_BIOMASS_CLUSTERS_H_
#define OBSERVATIONS_MORTALITY_EVENT_LENGTH_BIOMASS_CLUSTERS_H_

// headers
#include "Observations/Observation.h"
#include "Layers/Children/CategoricalLayer.h"

#include "AgeingErrors/AgeingError.h"
//#include "Processes/Children/Mortality.h"
#include "Processes/Children/Mortality/MortalityEventBiomass.h"
#include <boost/math/distributions/normal.hpp>


// namespaces
namespace niwa {
namespace observations {

using processes::MortalityEventBiomass;

/**
 * class definition
 */
class MortalityEventBiomassLengthClusters : public niwa::Observation {
public:
  // methods
  MortalityEventBiomassLengthClusters(Model* model);
  virtual                     ~MortalityEventBiomassLengthClusters();
  void                        DoValidate() override final;
  virtual void                DoBuild() override;
  void                        DoReset() override final { };
  void                        PreExecute() override final;
  void                        Execute() override final;
  void                        Simulate() override final;
  bool                        HasYear(unsigned year) const override final { return std::find(years_.begin(), years_.end(), year) != years_.end(); }
  virtual void                FillReportCache(ostringstream& cache);

protected:
  // Methods
  void                            ResetPreSimulation();

  // members
  vector<unsigned>                years_;
  vector<string>                  cells_;
  layers::CategoricalLayer*       layer_ = nullptr;
  string                          layer_label_;
  string                          allocation_ = PARAM_RANDOM;
  AllocationType                  allocation_type_ = AllocationType::kRandom;

  bool                            sexed_flag_;
  vector<unsigned>                length_bins_;

  // cluster info
  float                           average_cluster_weight_;
  float                           cluster_cv_;
  float                           cluster_sigma_;
  bool                            normal_clusters_ = true;
  unsigned                        length_samples_per_clusters_;
  vector<float>                   cluster_length_freq_;
  vector<unsigned>                agent_ndx_for_length_subsample_within_cluster_;
  vector<unsigned>                agent_ndx_for_age_subsample_within_cluster_;

  vector<unsigned>                agent_ndx_cluster_;
  vector<vector<unsigned>>        census_ndx_cluster_;
  vector<float>                   stratum_lf_;
  vector<unsigned>                expected_aged_length_freq_;
  vector<unsigned>                sampled_aged_length_freq_;
  //unsigned                        number_of_bootstraps_ = 0;

  //
  vector<float>                   target_length_distribution_;
  unsigned                        n_length_bins_;
  unsigned                        model_len_bins_;
  string                          obs_type_;
  bool                            are_obs_props_ = true;
  string                          fishery_label_;
  MortalityEventBiomass*          mortality_process_ = nullptr;
  vector<vector<processes::census_data>>* fishery_census_data_ = nullptr; // TODO: consider, the object that we want here, is dynamic in memory at DoExecute, so can't assign address until DoSimulate()

  string                          process_label_;
  string                          stratum_weight_method_ = PARAM_NONE;
  // TODO change from string -> unsigned int for a little speed up
  map<string,vector<unsigned>>    stratum_rows_;
  map<string,vector<unsigned>>    stratum_cols_;
  map<string,float>               stratum_area_;
  map<string,float>               stratum_biomass_;
  WorldView*                      world_ = nullptr;
  vector<unsigned>                fishery_years_;
  map<string,vector<float>>       stratum_length_frequency_;
  vector<vector<float>>           age_length_key_;

  parameters::Table*              cluster_sample_table_ = nullptr;

  map<unsigned, map<string, unsigned>>  cluster_by_year_and_stratum_;

  map<unsigned,map<string,vector<vector<float>>>>     age_length_key_by_year_stratum_;
  map<unsigned,map<string,vector<float>>>             lf_by_year_stratum_;
  vector<float>                                       model_length_mid_points_;

  // Storing cluster based information
  vector<vector<vector<float>>>               cluster_weight_; // n_years x n_stratum x n_clusters (weight of cluster)
  vector<vector<vector<float>>>               cluster_length_sample_weight_;
  vector<vector<vector<vector<unsigned>>>>    cluster_length_samples_; // n_years x n_stratum x n_clusters x n_samples
  vector<vector<vector<vector<unsigned>>>>    cluster_length_samples_sex_; // n_years x n_stratum x n_clusters x n_samples

  vector<vector<vector<float>>>               length_target_; // n_years x n_stratum x age_bins


};

} /* namespace observations */
} /* namespace niwa */

#endif /* OBSERVATIONS_MORTALITY_EVENT_LENGTH_BIOMASS_CLUSTERS_H_ */
