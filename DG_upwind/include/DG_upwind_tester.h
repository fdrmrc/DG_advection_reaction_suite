#ifndef DG_upwind_tester_h
#define DG_upwind_tester_h

#include <gtest/gtest.h>

#include <fstream>

#include "DG_upwind.h"

using namespace dealii;

// Test Fixture for Poisson problem, using integralconstant
template <class Integral>
class DG_upwind_Tester : public ::testing::Test, public AdvectionProblem<Integral::value>
{
public:

	DG_upwind_Tester() = default;
};

#endif
