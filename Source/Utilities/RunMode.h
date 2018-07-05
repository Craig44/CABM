/*
 * RunMode.h
 *
 *  Created on: 25/10/2013
 *      Author: Admin
 */

#ifndef UTILITIES_RUNMODE_H_
#define UTILITIES_RUNMODE_H_

// Enumerated Types
namespace RunMode {
enum Type {
  kInvalid      = 1,
  kLicense      = 2,
  kVersion      = 3,
  kHelp         = 4,
  kBasic        = 8,
};
}

#endif /* RUNMODE_H_ */
