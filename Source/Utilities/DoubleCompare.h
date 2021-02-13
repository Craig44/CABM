/*
 * DoubleCompare.h
 *
 *  Created on: 20/12/2012
 *      Author: Admin
 */

#ifndef DOUBLECOMPARE_H_
#define DOUBLECOMPARE_H_

#include "Utilities/Types.h"

// Namespaces
namespace niwa {
namespace utilities {
namespace doublecompare {

using niwa::utilities::Double;

// Defines
#define TEMP_ONE  1.0
#define TRUE_ZERO 0.0
#define TEMP_ZERO 1e-15
#define DELTA 1e-11

inline bool IsZero(const float &value) { return (value < TEMP_ZERO && value > -TEMP_ZERO); }
//inline bool IsInfinite(const Double &value) { return (isinf(value));}
inline bool IsTrueZero(const float &value) { return (value < TRUE_ZERO && value > -TRUE_ZERO); }
inline bool IsOne(const float &value) { return ( ((value-TEMP_ONE) < TEMP_ZERO) && ((value-TEMP_ONE) > -TEMP_ZERO) ); }
inline bool IsEqual(float A, Double B) { return ( ((A-B) < TEMP_ZERO) && ((A-B) > -TEMP_ZERO) ); }

inline niwa::utilities::Double ZeroFun(float x) {
  if (x >= TEMP_ZERO)
    return x;

  return TEMP_ZERO / (2.0 - (x / TEMP_ZERO));
}

inline niwa::utilities::Double ZeroFun(float x, float delta) {
  if (x >= delta)
    return x;

  return delta / (2.0 - (x / delta));
}


} /* namespace doublecompare */
} /* namespace utilities */
} /* namespace niwa */



//  static bool     isZero(double A) { return ( (A < ZERO) && (A > -ZERO) ); }
//  static bool     isTrueZero(double A) { return ( (A < TRUE_ZERO) && (A > -TRUE_ZERO) ); }
//  static bool     isNonNegative(double A) { return ( 0.0 <= A ); }
//  static bool     isPositive(double A) { return ( 0.0 < A ); }
//  static bool     isEqual(double A, double B) { return ( ((A-B) < ZERO) && ((A-B) > -ZERO) ); }
//  static bool     isBetween(double A, double B, double C) {
//    return ( ((A-B) > -ZERO) && ((A-C) < ZERO) );
//  }


#endif /* DOUBLECOMPARE_H_ */
