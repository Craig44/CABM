/**
 * @file Factory.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 13/12/2012
 * @section LICENSE
 *
 * Copyright NIWA Science �2012 - www.niwa.co.nz
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */

// Headers
#include "Factory.h"

#include "Model/Model.h"
#include "Model/Managers.h"
#include "Processes/Manager.h"

#include "Children/Nop.h"
#include "Children/Ageing.h"
#include "Children/Recruitment/RecruitmentBevertonHolt.h"
#include "Children/Recruitment/RecruitmentConstant.h"
#include "Children/Growth/GrowthVonBertalanffyWithBasic.h"
#include "Children/Growth/GrowthSchnuteWithBasic.h"
#include "Children/Mortality/MortalityConstantRate.h"
#include "Children/Mortality/MortalityCull.h"
#include "Children/Mortality/MortalityExploitation.h"
#include "Children/Mortality/MortalityEventBiomass.h"
#include "Children/Mortality/MortalityEventHybrid.h"
#include "Children/Mortality/MortalityEffortBased.h"
#include "Children/Mortality/MortalityEffortBasedWithCovar.h"

#include "Children/Movement/MovementBoxTransfer.h"
#include "Children/Movement/MovementBoxTransferDensity.h"
#include "Children/Movement/MovementPreference.h"
#include "Children/Tagging.h"
#include "Children/TagShedding.h"
#include "Children/Maturity.h"
#include "Children/Mortality/MortalityBaranovSimple.h"
#include "Children/Mortality/MortalityBaranov.h"

// Namespaces
namespace niwa {
namespace processes {

/**
 * Create the instance of our object as defined by the two parameters
 * object_type and sub_type.
 *
 * @param object_type The type of object to create (e.g age_size, process)
 * @param sub_type The child type of the object to create (e.g ageing, schnute)
 * @return shared_ptr to the object we've created
 */
Process* Factory::Create(Model* model, const string& object_type, const string& sub_type) {
  Process* result = nullptr;

  string object = object_type;
  string sub    = sub_type;

  /**
   * If object_type is not "process" or "processes" then we're using a special
   * declaration (e.g @maturation) and we want to modify this back to the standard
   * method so we only need 1 set of conditional statements.
   */
  if (object != PARAM_PROCESS && object != PARAM_PROCESSES) {
    LOG_FINE() << "Changing object_type (" << object << ") and sub_type (" << ") to the standard declaration format";
    if (sub != "")
      sub = object_type + "_" + sub_type;
    else
      sub = object_type;

    object = PARAM_PROCESS;

    LOG_FINE() << "Finished modification of object_type (" << object << ") and sub_type (" << sub << ")";
  }

  if (object == PARAM_PROCESS || object == PARAM_PROCESSES) {
    if (sub == PARAM_NOP)
          result = new Nop(model);
    else if (sub == PARAM_RECRUITMENT_BEVERTON_HOLT)
      result = new RecruitmentBevertonHolt(model);
    else if (sub == PARAM_AGEING)
      result = new Ageing(model);
    else if (sub == PARAM_RECRUITMENT_CONSTANT)
      result = new RecruitmentConstant(model);
    else if (sub == PARAM_GROWTH_VON_BERTALANFFY_WITH_BASIC)
      result = new GrowthVonBertalanffyWithBasic(model);
    else if (sub == PARAM_GROWTH_SCHNUTE_WITH_BASIC)
      result = new GrowthSchnuteWithBasic(model);
    else if (sub == PARAM_MORTALITY_BARANOV)
      result = new MortalityBaranov(model);
    else if (sub == PARAM_MORTALITY_EXPLOITATION)
      result = new MortalityExploitation(model);
    else if (sub == PARAM_MORTALITY_BARANOV_SIMPLE)
      result = new MortalityBaranovSimple(model);
    else if (sub == PARAM_MORTALITY_CONSTANT_RATE)
      result = new MortalityConstantRate(model);
    else if (sub == PARAM_MORTALITY_EVENT_BIOMASS)
      result = new MortalityEventBiomass(model);
    else if (sub == PARAM_MORTALITY_EVENT_HYBRID)
      result = new MortalityEventHybrid(model);
    else if (sub == PARAM_MORTALITY_CULL)
      result = new MortalityCull(model);
    else if (sub == PARAM_MORTALITY_EFFORT_BASED)
      result = new MortalityEffortBased(model);
    else if (sub == PARAM_MORTALITY_EFFORT_WITH_COVARIATES)
      result = new MortalityEffortBasedWithCovar(model);
    else if (sub == PARAM_MOVEMENT_BOX_TRANSFER)
      result = new MovementBoxTransfer(model);
    else if (sub == PARAM_MOVEMENT_BOX_TRANSFER_DENSITY_TRIGGER)
      result = new MovementBoxTransferDensity(model);
    else if (sub == PARAM_PREFERENCE_MOVEMENT)
      result = new MovementPreference(model);
    else if (sub == PARAM_MATURATION)
      result = new Maturity(model);
    else if (sub == PARAM_TAGGING)
      result = new Tagging(model);
    else if (sub == PARAM_TAG_SHEDDING)
      result = new TagShedding(model);
    if (result)
      model->managers().process()->AddObject(result);
  }

  return result;
}

} /* namespace processes */
} /* namespace niwa */
