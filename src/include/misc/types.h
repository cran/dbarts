/* src/include/misc/types.h.  Generated from types.h.in by configure.  */
#ifndef MISC_TYPES_H
#define MISC_TYPES_H

#include <stdint.h> // type should be replaced with a uintXX_t

#ifndef XINT_TYPE
/* Define to the integer type to be used to map the floating point x to
   integer types. Essentially limits the number of cut points allowable. */
#define XINT_TYPE uint16_t
#endif

typedef XINT_TYPE misc_xint_t;

#endif // MISC_TYPES_H

