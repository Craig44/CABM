/*
 * Manager.cpp
 *
 *  Created on: 13/12/2012
 *      Author: Admin
 */

#include "Manager.h"

#include "Model/Model.h"
#include "Logging/Logging.h"


namespace niwa {
namespace processes {

Manager::Manager() {
}

Manager::~Manager() noexcept(true) {
}

/**
 * Validate any loaded processes we have.
 */
void Manager::Validate() {
  LOG_TRACE();
  LOG_CODE_ERROR() << "This method is not supported";
}

void Manager::Validate(Model* model) {
  LOG_TRACE();
  base::Manager<niwa::processes::Manager, niwa::Process>::Validate();

  if (objects_.size() == 0)
    LOG_FATAL() << "The configuration file requires you specify at least one type of process. E.g @recruitment, @mortality, @ageing";

  for (auto process : objects_) {
    LOG_FINEST() << "processes managed" << process->label();
  }
}


/**
 * override Base classes and Build non mortality and growth processes
 */
void Manager::BuildRemainingProcesses() {
  LOG_FINEST() << "Starting Build... with " << objects_.size() << " objects";
  for(auto stored_object : objects_) {
    if ((stored_object->process_type() != ProcessType::kGrowth) && (stored_object->process_type() != ProcessType::kMortality)) {
      LOG_FINEST() << "Building process = " << stored_object->label();
      stored_object->Build();
    }
  }
  LOG_FINEST() << "Build Finished";
}

/**
 * override Base classes and Build mortality and growth processes
 */
void Manager::BuildGrowthAndMortalityProcesses() {
  LOG_FINEST() << "Starting Build... with " << objects_.size() << " objects";
  for(auto stored_object : objects_) {
    if ((stored_object->process_type() == ProcessType::kGrowth) || (stored_object->process_type() == ProcessType::kMortality)) {
      stored_object->Build();
      LOG_FINEST() << "Building process = " << stored_object->label();
    }
  }
  LOG_FINEST() << "Build Finished";
}

/**
 * Return the process with the name passed in as a parameter.
 * If no process is found then an empty pointer will
 * be returned.
 *
 * @param label The name of the process to find
 * @return A pointer to the process or empty pointer
 */
Process* Manager::GetProcess(const string& label) {
  for (auto process : objects_) {
    if (process->label() == label)
      return process;
  }

  return nullptr;
}

/**
 * Return the growth process with the name passed in as a parameter.
 * If no process is found then an empty pointer will
 * be returned.
 *
 * @param label The name of the process to find
 * @return A pointer to the process or empty pointer
 */
Growth* Manager::GetGrowthProcess(const string& label) {
  LOG_TRACE();
  Growth* growth_ptr = nullptr;
  for (auto process : objects_) {
    if ((process->label() == label) && (process->process_type() == ProcessType::kGrowth)) {
      growth_ptr = dynamic_cast<Growth*>(process);
      return growth_ptr;
    }
  }

  return growth_ptr;
}

/**
 * Return the growth process with the name passed in as a parameter.
 * If no process is found then an empty pointer will
 * be returned.
 *
 * @param label The name of the process to find
 * @return A pointer to the process or empty pointer
 */
Mortality* Manager::GetMortalityProcess(const string& label) {
  LOG_TRACE();
  Mortality* mortality_ptr = nullptr;
  for (auto process : objects_) {
    if ((process->label() == label) && (process->process_type() == ProcessType::kMortality)) {
      mortality_ptr = dynamic_cast<Mortality*>(process);
      return mortality_ptr;
    }
  }

  return mortality_ptr;
}
} /* namespace processes */
} /* namespace niwa */
