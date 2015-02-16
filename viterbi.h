/*
 * vcftools.h
 */

#ifndef VCFTOOLS_H_
#define VCFTOOLS_H_

#include <algorithm>
#include <functional>   // std::plus
#include <set>
#include <sstream>
#include <thread>

#include "output_log.h"
#include "parameters.h"
#include "bcf_file.h"
#include "vcf_file.h"
#include "variant_file.h"
#include "header.h"

#include <ctime>
#include <sys/time.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif


using namespace std;

#endif /* VCFTOOLS_H_ */
