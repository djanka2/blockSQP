/***************************************************************************
   VPLAN Source Code
   Version: $Id$
   Author: Cetin Sert
 ***************************************************************************/

/**
 * \file fornorm.h
 * \author Stefan Koerkel
 * \brief Macros for fortran function name normalizations
 */

#ifndef FORNORM_H
#define FORNORM_H

#define ___CONCAT(x,y) x##y

#if !defined NO_UNDERSCORING
#define __APPEND_UNDERSCORE(x) ___CONCAT(x,_)
#else
#define __APPEND_UNDERSCORE(x) x
#endif

#ifdef UPPERCASE
#define __FORNORM2(l,U) __APPEND_UNDERSCORE(l)
#else
#define __FORNORM2(l,U) __APPEND_UNDERSCORE(U)
#endif

#define __FORNORM(l) __APPEND_UNDERSCORE(l)

#define __FORNORML(l,U) __APPEND_UNDERSCORE(l)
#define __FORNORMB(l,U) __APPEND_UNDERSCORE(l)

//#define __FORNORML(l,U) U
//#define __FORNORMB(l,U) U

#endif
