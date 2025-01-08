//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2022, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================

/*
  Unit tests for the CKKS scheme
 */

#include "UnitTestUtils.h"
#include "UnitTestCCParams.h"
#include "UnitTestCryptoContext.h"
#include "scheme/ckksrns/ckksrns-utils.h"
#include "cryptocontext-ser.h"
#include "scheme/ckksrns/ckksrns-ser.h"

#include <iostream>
#include <vector>
#include "gtest/gtest.h"
#include <iterator>
#include <string>

using namespace lbcrypto;

//===========================================================================================================
enum TEST_CASE_TYPE {
    BOOTSTRAP_FULL = 0,
    BOOTSTRAP_EDGE,
    BOOTSTRAP_SPARSE,
    BOOTSTRAP_KEY_SWITCH,
    BOOTSTRAP_ITERATIVE,
    BOOTSTRAP_NUM_TOWERS,
    BOOTSTRAP_SERIALIZE,
};

static std::ostream& operator<<(std::ostream& os, const TEST_CASE_TYPE& type) {
    std::string typeName;
    switch (type) {
        case BOOTSTRAP_FULL:
            typeName = "BOOTSTRAP_FULL";
            break;
        case BOOTSTRAP_EDGE:
            typeName = "BOOTSTRAP_EDGE";
            break;
        case BOOTSTRAP_SPARSE:
            typeName = "BOOTSTRAP_SPARSE";
            break;
        case BOOTSTRAP_KEY_SWITCH:
            typeName = "BOOTSTRAP_KEY_SWITCH";
            break;
        case BOOTSTRAP_ITERATIVE:
            typeName = "BOOTSTRAP_ITERATIVE";
            break;
        case BOOTSTRAP_NUM_TOWERS:
            typeName = "BOOTSTRAP_NUM_TOWERS";
            break;
        case BOOTSTRAP_SERIALIZE:
            typeName = "BOOTSTRAP_SERIALIZE";
            break;
        default:
            typeName = "UNKNOWN";
            break;
    }
    return os << typeName;
}
//===========================================================================================================
struct TEST_CASE_UTCKKSRNSCS_BOOT {
    TEST_CASE_TYPE testCaseType;
    // test case description - MUST BE UNIQUE
    std::string description;

    UnitTestCCParams params;

    // additional test case data
    // ........
    std::vector<uint32_t> levelBudget;
    std::vector<uint32_t> dim1;
    uint32_t slots;

    std::string buildTestName() const {
        std::stringstream ss;
        ss << testCaseType << "_" << description;
        return ss.str();
    }
    std::string toString() const {
        std::stringstream ss;
        ss << "testCaseType [" << testCaseType << "], " << params.toString();
        return ss.str();
    }
};

// this lambda provides a name to be printed for every test run by INSTANTIATE_TEST_SUITE_P.
// the name MUST be constructed from digits, letters and '_' only
static auto testName = [](const testing::TestParamInfo<TEST_CASE_UTCKKSRNSCS_BOOT>& test) {
    return test.param.buildTestName();
};

static std::ostream& operator<<(std::ostream& os, const TEST_CASE_UTCKKSRNSCS_BOOT& test) {
    return os << test.toString();
}
//===========================================================================================================
constexpr uint32_t MULT_DEPTH   = 25;
constexpr uint32_t RDIM         = 64;
constexpr uint32_t NUM_LRG_DIGS = 3;

// to test for composite scaling degree d = 3
constexpr uint32_t SMODSIZED3 = 78;
constexpr uint32_t FMODSIZED3 = 89;

// to test for composite scaling degree d = 2
constexpr uint32_t SMODSIZED2 = 59;
constexpr uint32_t FMODSIZED2 = 60;

// clang-format off
static std::vector<TEST_CASE_UTCKKSRNSCS_BOOT> testCases = {
    // TestType,     Descr, Scheme,          RDim, MultDepth,  SModSize,     DSize, BatchSz, SecKeyDist,      MaxRelinSkDeg, FModSize,  SecLvl,       KSTech, ScalTech,        LDigits,      PtMod, StdDev, EvalAddCt, KSCt, MultTech, EncTech, PREMode, multipartyMode, decryptionNoiseMode, executionMode, noiseEstimate, registerWordSize, LvlBudget, Dim1,       Slots
#if NATIVEINT == 64
    { BOOTSTRAP_FULL, "01", {CKKSRNS_SCHEME, RDIM, MULT_DEPTH, SMODSIZED2,     DFLT,  DFLT,    UNIFORM_TERNARY, DFLT,          FMODSIZED2,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 1, 1 },  { 32, 32 }, RDIM/2 },
    { BOOTSTRAP_FULL, "02", {CKKSRNS_SCHEME, RDIM, MULT_DEPTH, SMODSIZED3,     DFLT,  DFLT,    UNIFORM_TERNARY, DFLT,          FMODSIZED3,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 1, 1 },  { 32, 32 }, RDIM/2 },
#endif
    // ==========================================
    // TestType,     Descr, Scheme,          RDim, MultDepth,  SModSize,     DSize, BatchSz, SecKeyDist,      MaxRelinSkDeg, FModSize,  SecLvl,       KSTech, ScalTech,        LDigits,      PtMod, StdDev, EvalAddCt, KSCt, MultTech, EncTech, PREMode, multipartyMode, decryptionNoiseMode, executionMode, noiseEstimate, registerWordSize, LvlBudget, Dim1,     Slots
#if NATIVEINT == 64
    { BOOTSTRAP_FULL, "11", {CKKSRNS_SCHEME, RDIM, MULT_DEPTH, SMODSIZED2,     DFLT,  DFLT,    UNIFORM_TERNARY, DFLT,          FMODSIZED2,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 3, 3 },  { 0, 0 }, RDIM/2 },
    { BOOTSTRAP_FULL, "12", {CKKSRNS_SCHEME, RDIM, MULT_DEPTH, SMODSIZED3,     DFLT,  DFLT,    UNIFORM_TERNARY,  DFLT,         FMODSIZED3,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 3, 3 },  { 0, 0 }, RDIM/2 },
#endif
    // ==========================================
    // TestType,      Descr, Scheme,          RDim, MultDepth,  SModSize,     DSize, BatchSz, SecKeyDist,      MaxRelinSkDeg, FModSize,  SecLvl,       KSTech, ScalTech,        LDigits,      PtMod, StdDev, EvalAddCt, KSCt, MultTech, EncTech, PREMode, multipartyMode, decryptionNoiseMode, executionMode, noiseEstimate, registerWordSize, LvlBudget, Dim1,     Slots
#if NATIVEINT == 64
    { BOOTSTRAP_EDGE, "01", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED2,     DFLT,  DFLT,    UNIFORM_TERNARY, DFLT,          FMODSIZED2,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 1, 1 },  { 0, 0 }, RDIM/4 },
    { BOOTSTRAP_EDGE, "02", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED3,     DFLT,  DFLT,    UNIFORM_TERNARY, DFLT,          FMODSIZED3,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 1, 1 },  { 0, 0 }, RDIM/4 },
#endif
    // ==========================================
    // TestType,      Descr, Scheme,          RDim, MultDepth,  SModSize,     DSize, BatchSz, SecKeyDist,      MaxRelinSkDeg, FModSize,  SecLvl,       KSTech, ScalTech,        LDigits,      PtMod, StdDev, EvalAddCt, KSCt, MultTech, EncTech, PREMode, multipartyMode, decryptionNoiseMode, executionMode, noiseEstimate, registerWordSize, LvlBudget, Dim1,     Slots
#if NATIVEINT == 64
    { BOOTSTRAP_EDGE, "11", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED2,     DFLT,  DFLT,    UNIFORM_TERNARY, DFLT,          FMODSIZED2,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 2, 2 },  { 0, 0 }, RDIM/4 },
    { BOOTSTRAP_EDGE, "12", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED3,     DFLT,  DFLT,    UNIFORM_TERNARY,  DFLT,         FMODSIZED3,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 2, 2 },  { 0, 0 }, RDIM/4 },
#endif
    // ==========================================
    // TestType,        Descr, Scheme,          RDim, MultDepth,  SModSize,     DSize, BatchSz, SecKeyDist,      MaxRelinSkDeg, FModSize,  SecLvl,       KSTech, ScalTech,        LDigits,      PtMod, StdDev, EvalAddCt, KSCt, MultTech, EncTech, PREMode, multipartyMode, decryptionNoiseMode, executionMode, noiseEstimate, registerWordSize, LvLBudget, Dim1,     Slots
#if NATIVEINT == 64
    { BOOTSTRAP_SPARSE, "01", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED2,     DFLT,  DFLT,    UNIFORM_TERNARY, DFLT,          FMODSIZED2,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 1, 1 },  { 8, 8 }, 8 },
    { BOOTSTRAP_SPARSE, "02", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED3,     DFLT,  DFLT,    UNIFORM_TERNARY,  DFLT,         FMODSIZED3,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 1, 1 },  { 8, 8 }, 8 },
#endif
    // ==========================================
    // TestType,        Descr, Scheme,          RDim, MultDepth,  SModSize,     DSize, BatchSz, SecKeyDist,      MaxRelinSkDeg, FModSize,  SecLvl,       KSTech, ScalTech,        LDigits,      PtMod, StdDev, EvalAddCt, KSCt, MultTech, EncTech, PREMode, multipartyMode, decryptionNoiseMode, executionMode, noiseEstimate, registerWordSize, LvLBudget, Dim1,     Slots
#if NATIVEINT == 64
    { BOOTSTRAP_SPARSE, "11", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED2,     DFLT,  DFLT,    UNIFORM_TERNARY, DFLT,          FMODSIZED2,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 2, 2 },  { 0, 0 }, 8 },
    { BOOTSTRAP_SPARSE, "12", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED3,     DFLT,  DFLT,    UNIFORM_TERNARY,  DFLT,          FMODSIZED3,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 2, 2 },  { 0, 0 }, 8 },
 #endif
    // ==========================================
    // TestType,        Descr, Scheme,          RDim, MultDepth,  SModSize,     DSize, BatchSz, SecKeyDist,      MaxRelinSkDeg, FModSize,  SecLvl,       KSTech, ScalTech,        LDigits,      PtMod, StdDev, EvalAddCt, KSCt, MultTech, EncTech, PREMode, multipartyMode, decryptionNoiseMode, executionMode, noiseEstimate, registerWordSize, LvLBudget, Dim1,     Slots
#if NATIVEINT == 64
    { BOOTSTRAP_SPARSE, "21", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED2,     DFLT,  DFLT,    UNIFORM_TERNARY, DFLT,          FMODSIZED2,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 1, 1 },  { 0, 0 }, 8 },
    { BOOTSTRAP_SPARSE, "22", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED3,     DFLT,  DFLT,    UNIFORM_TERNARY,  DFLT,          FMODSIZED3,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 1, 1 },  { 0, 0 }, 8 },
#endif
    // ==========================================
    // TestType,        Descr, Scheme,          RDim, MultDepth,  SModSize,     DSize, BatchSz, SecKeyDist,      MaxRelinSkDeg, FModSize,  SecLvl,       KSTech, ScalTech,        LDigits,      PtMod, StdDev, EvalAddCt, KSCt, MultTech, EncTech, PREMode, multipartyMode, decryptionNoiseMode, executionMode, noiseEstimate, registerWordSize, LvLBudget, Dim1,     Slots
#if NATIVEINT == 64
    { BOOTSTRAP_SPARSE, "31", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED2,     DFLT,  DFLT,    UNIFORM_TERNARY, DFLT,          FMODSIZED2,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 3, 2 },  { 0, 0 }, 8 },
    { BOOTSTRAP_SPARSE, "32", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED3,     DFLT,  DFLT,    UNIFORM_TERNARY,  DFLT,         FMODSIZED3,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 3, 2 },  { 0, 0 }, 8 },
#endif

#if NATIVEINT == 64
    { BOOTSTRAP_SPARSE, "39", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED2,     DFLT,  DFLT,    UNIFORM_TERNARY, DFLT,          FMODSIZED2,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 1, 1 },  { 0, 0 }, 1 },
    { BOOTSTRAP_SPARSE, "40", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED3,     DFLT,  DFLT,    UNIFORM_TERNARY,  DFLT,          FMODSIZED3,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 1, 1 },  { 0, 0 }, 1 },
#endif

    // ==========================================
    // TestType,            Descr, Scheme,          RDim, MultDepth,  SModSize,     DSize, BatchSz, SecKeyDist,      MaxRelinSkDeg, FModSize,  SecLvl,       KSTech, ScalTech,        LDigits,      PtMod, StdDev, EvalAddCt, KSCt, MultTech, EncTech, PREMode, multipartyMode, decryptionNoiseMode, executionMode, noiseEstimate, registerWordSize, LvLBudget, Dim1
#if NATIVEINT == 64
    { BOOTSTRAP_KEY_SWITCH, "01", {CKKSRNS_SCHEME,  2048, MULT_DEPTH, SMODSIZED2,     DFLT,  8,       UNIFORM_TERNARY,  DFLT,          FMODSIZED2,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 3, 2 },  { 0, 0 } },
    { BOOTSTRAP_KEY_SWITCH, "02", {CKKSRNS_SCHEME,  2048, MULT_DEPTH, SMODSIZED3,     DFLT,  8,       UNIFORM_TERNARY, DFLT,          FMODSIZED3,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 3, 2 },  { 0, 0 } },
#endif
    // ==========================================
    // TestType,           Descr, Scheme,          RDim, MultDepth,  SModSize,     DSize, BatchSz, SecKeyDist,      MaxRelinSkDeg, FModSize,  SecLvl,       KSTech, ScalTech,        LDigits,      PtMod, StdDev, EvalAddCt, KSCt, MultTech, EncTech, PREMode, multipartyMode, decryptionNoiseMode, executionMode, noiseEstimate, registerWordSize, LvLBudget, Dim1,     Slots
#if NATIVEINT == 64
    { BOOTSTRAP_ITERATIVE, "01", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED2,     DFLT,  8,       UNIFORM_TERNARY,  DFLT,          FMODSIZED2,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 3, 3 },  { 0, 0 }, 8},
    { BOOTSTRAP_ITERATIVE, "02", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED3,     DFLT,  8,       UNIFORM_TERNARY, DFLT,          FMODSIZED3,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 3, 3 },  { 0, 0 }, 8},
#endif
    // TestType,           Descr, Scheme,          RDim, MultDepth,  SModSize,     DSize, BatchSz, SecKeyDist,      MaxRelinSkDeg, FModSize,  SecLvl,       KSTech, ScalTech,        LDigits,      PtMod, StdDev, EvalAddCt, KSCt, MultTech, EncTech, PREMode, multipartyMode, decryptionNoiseMode, executionMode, noiseEstimate, registerWordSize, LvLBudget, Dim1,     Slots
#if NATIVEINT == 64
    { BOOTSTRAP_ITERATIVE, "09", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED2,     DFLT,  DFLT,    UNIFORM_TERNARY,  DFLT,         FMODSIZED2,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 3, 3 },  { 0, 0 }, RDIM/2},
    { BOOTSTRAP_ITERATIVE, "10", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED3,     DFLT,  DFLT,    UNIFORM_TERNARY, DFLT,          FMODSIZED3,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 3, 3 },  { 0, 0 }, RDIM/2},
#endif
    // ==========================================
    // TestType,           Descr, Scheme,          RDim, MultDepth,  SModSize,     DSize, BatchSz, SecKeyDist,      MaxRelinSkDeg, FModSize,  SecLvl,       KSTech, ScalTech,        LDigits,      PtMod, StdDev, EvalAddCt, KSCt, MultTech, EncTech, PREMode, multipartyMode, decryptionNoiseMode, executionMode, noiseEstimate, registerWordSize, LvLBudget, Dim1,     Slots
#if NATIVEINT == 64
    { BOOTSTRAP_NUM_TOWERS, "01", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED2,     DFLT,  8,       UNIFORM_TERNARY, DFLT,         FMODSIZED2,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 3, 2 },  { 0, 0 }, 8},
    { BOOTSTRAP_NUM_TOWERS, "02", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED3,     DFLT,  8,       UNIFORM_TERNARY, DFLT,          FMODSIZED3,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 3, 2 },  { 0, 0 }, 8},
#endif
    // TestType,            Descr, Scheme,          RDim, MultDepth,  SModSize,     DSize, BatchSz, SecKeyDist,      MaxRelinSkDeg, FModSize,  SecLvl,       KSTech, ScalTech,        LDigits,      PtMod, StdDev, EvalAddCt, KSCt, MultTech, EncTech, PREMode, multipartyMode, decryptionNoiseMode, executionMode, noiseEstimate, registerWordSize, LvLBudget, Dim1,     Slots
#if NATIVEINT == 64
    { BOOTSTRAP_NUM_TOWERS, "09", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED2,     DFLT,  DFLT,    UNIFORM_TERNARY, DFLT,         FMODSIZED2,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 3, 2 },  { 0, 0 }, RDIM/2},
    { BOOTSTRAP_NUM_TOWERS, "10", {CKKSRNS_SCHEME,  RDIM, MULT_DEPTH, SMODSIZED3,     DFLT,  DFLT,    UNIFORM_TERNARY, DFLT,          FMODSIZED3,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 3, 2 },  { 0, 0 }, RDIM/2},
#endif
    // ==========================================
    // TestType,           Descr, Scheme,         RDim, MultDepth,  SModSize,     DSize, BatchSz, SecKeyDist,      MaxRelinSkDeg, FModSize,  SecLvl,       KSTech, ScalTech,        LDigits,      PtMod, StdDev, EvalAddCt, KSCt, MultTech, EncTech, PREMode, multipartyMode, decryptionNoiseMode, executionMode, noiseEstimate, registerWordSize, LvLBudget, Dim1,       Slots
#if NATIVEINT == 64
    { BOOTSTRAP_SERIALIZE, "01", {CKKSRNS_SCHEME, RDIM, MULT_DEPTH, SMODSIZED2,     DFLT,  DFLT,    UNIFORM_TERNARY, DFLT,          FMODSIZED2,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 1, 1 },  { 32, 32 }, RDIM/2 },
    { BOOTSTRAP_SERIALIZE, "02", {CKKSRNS_SCHEME, RDIM, MULT_DEPTH, SMODSIZED3,     DFLT,  DFLT,    UNIFORM_TERNARY, DFLT,          FMODSIZED3,  HEStd_NotSet, HYBRID, COMPOSITESCALINGAUTO,       NUM_LRG_DIGS, DFLT,  DFLT,   DFLT,      DFLT, DFLT,     DFLT,    DFLT, DFLT, DFLT, DFLT, DFLT, 32},   { 1, 1 },  { 32, 32 }, RDIM/2 },
#endif
    // ==========================================
};
// clang-format on
//===========================================================================================================
class UTCKKSRNSCS_BOOT : public ::testing::TestWithParam<TEST_CASE_UTCKKSRNSCS_BOOT> {
    using Element = DCRTPoly;

    // The precision after which we consider two values equal.
    // This is necessary because CKKS works for approximate numbers.
    const double eps = 0.0001;

    // CalculateApproximationError() calculates the precision number (or approximation error).
    // The higher the precision, the less the error.
    double CalculateApproximationError(const std::vector<std::complex<double>>& result,
                                       const std::vector<std::complex<double>>& expectedResult) {
        if (result.size() != expectedResult.size())
            OPENFHE_THROW("Cannot compare vectors with different numbers of elements");

        // using the infinity norm
        double maxError = 0;
        for (size_t i = 0; i < result.size(); ++i) {
            double error = std::abs(result[i].real() - expectedResult[i].real());
            if (maxError < error)
                maxError = error;
        }

        return std::abs(std::log2(maxError));
    }

protected:
    void SetUp() {}

    void TearDown() {
        CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
    }

    void UnitTest_Bootstrap(const TEST_CASE_UTCKKSRNSCS_BOOT& testData, const std::string& failmsg = std::string()) {
        try {
            CryptoContext<Element> cc(UnitTestGenerateContext(testData.params));

            cc->EvalBootstrapSetup(testData.levelBudget, testData.dim1, testData.slots);

            auto keyPair = cc->KeyGen();
            cc->EvalBootstrapKeyGen(keyPair.secretKey, testData.slots);
            cc->EvalAtIndexKeyGen(keyPair.secretKey, {6});
            cc->EvalMultKeyGen(keyPair.secretKey);

            std::vector<std::complex<double>> input;
            if (testData.slots < 8) {
                input = Fill({0.1415926}, testData.slots);
            }
            else {
                input = Fill({0.111111, 0.222222, 0.333333, 0.444444, 0.555555, 0.666666, 0.777777, 0.888888},
                             testData.slots);
            }

            size_t encodedLength = input.size();

            Plaintext plaintext1 = cc->MakeCKKSPackedPlaintext(input, 1, MULT_DEPTH - 1, nullptr, testData.slots);
            auto ciphertext1     = cc->Encrypt(keyPair.publicKey, plaintext1);
            auto ciphertextAfter = cc->EvalBootstrap(ciphertext1);

            Plaintext result;
            cc->Decrypt(keyPair.secretKey, ciphertextAfter, &result);
            result->SetLength(encodedLength);
            plaintext1->SetLength(encodedLength);
            checkEquality(result->GetCKKSPackedValue(), plaintext1->GetCKKSPackedValue(), eps,
                          failmsg + " Bootstrapping for fully packed ciphertexts fails");

            int32_t rotIndex = (testData.slots < 8) ? 0 : 6;
            auto temp6       = input;
            std::rotate(temp6.begin(), temp6.begin() + rotIndex, temp6.end());

            auto ciphertext6 = cc->EvalAtIndex(ciphertextAfter, rotIndex);
            Plaintext result6;
            cc->Decrypt(keyPair.secretKey, ciphertext6, &result6);
            result6->SetLength(encodedLength);
            checkEquality(result6->GetCKKSPackedValue(), temp6, eps,
                          failmsg + " EvalAtIndex after Bootstrapping for fully packed ciphertexts fails");
        }
        catch (std::exception& e) {
            std::cerr << "Exception thrown from " << __func__ << "(): " << e.what() << std::endl;
            // make it fail
            EXPECT_TRUE(0 == 1) << failmsg;
        }
        catch (...) {
            UNIT_TEST_HANDLE_ALL_EXCEPTIONS;
        }
    }

    void UnitTest_Bootstrap_KeySwitching(const TEST_CASE_UTCKKSRNSCS_BOOT& testData,
                                         const std::string& failmsg = std::string()) {
        try {
            CryptoContext<Element> cc(UnitTestGenerateContext(testData.params));

            // We set the cout precision to 8 decimal digits for a nicer output.
            // If you want to see the error/noise introduced by CKKS, bump it up
            // to 15 and it should become visible.
            // std::cout.precision(8);
            auto keyPair = cc->KeyGen();
            cc->EvalAtIndexKeyGen(keyPair.secretKey, {1});

            double eps                          = 0.00000001;
            std::vector<std::complex<double>> a = {0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 4.0, 5.0};
            std::vector<std::complex<double>> b = {0.5, 0.75, 1.0, 2.0, 3.0, 4.0, 5.0};
            Plaintext plaintext_a               = cc->MakeCKKSPackedPlaintext(a);
            auto comp_a                         = plaintext_a->GetCKKSPackedValue();
            Plaintext plaintext_b               = cc->MakeCKKSPackedPlaintext(b);
            auto comp_b                         = plaintext_b->GetCKKSPackedValue();

            // Test for KeySwitchExt + KeySwitchDown
            auto ciphertext = cc->Encrypt(keyPair.publicKey, plaintext_a);
            ciphertext      = cc->KeySwitchExt(ciphertext, true);
            ciphertext      = cc->KeySwitchDown(ciphertext);

            Plaintext result;
            cc->Decrypt(keyPair.secretKey, ciphertext, &result);
            result->SetLength(a.size());
            checkEquality(result->GetCKKSPackedValue(), comp_a, eps,
                          failmsg + " Bootstrapping for KeySwitchExt + KeySwitchDown failed");

            // Test for EvalFastRotationExt
            ciphertext  = cc->Encrypt(keyPair.publicKey, plaintext_a);
            auto digits = cc->EvalFastRotationPrecompute(ciphertext);
            ciphertext  = cc->EvalFastRotationExt(ciphertext, 1, digits, true);
            ciphertext  = cc->KeySwitchDown(ciphertext);

            cc->Decrypt(keyPair.secretKey, ciphertext, &result);
            result->SetLength(b.size());
            checkEquality(result->GetCKKSPackedValue(), comp_b, eps,
                          failmsg + " Bootstrapping for EvalFastRotationExt failed");

            // Test for KeySwitchExt + KeySwitchDown w/o first element
            ciphertext        = cc->Encrypt(keyPair.publicKey, plaintext_a);
            auto firstCurrent = ciphertext->GetElements()[0];
            ciphertext        = cc->KeySwitchExt(ciphertext, false);
            ciphertext        = cc->KeySwitchDown(ciphertext);
            auto elements     = ciphertext->GetElements();
            elements[0] += firstCurrent;
            ciphertext->SetElements(elements);

            cc->Decrypt(keyPair.secretKey, ciphertext, &result);
            result->SetLength(a.size());
            checkEquality(result->GetCKKSPackedValue(), comp_a, eps,
                          failmsg + " Bootstrapping for KeySwitchExt + KeySwitchDown w/o first element failed");

            // Test for EvalFastRotationExt w/o first element
            ciphertext   = cc->Encrypt(keyPair.publicKey, plaintext_a);
            firstCurrent = ciphertext->GetElements()[0];
            // Find the automorphism index that corresponds to rotation index index.
            usint autoIndex = FindAutomorphismIndex2nComplex(1, 4096);
            std::vector<usint> map(4096 / 2);
            PrecomputeAutoMap(4096 / 2, autoIndex, &map);
            firstCurrent = firstCurrent.AutomorphismTransform(autoIndex, map);
            digits       = cc->EvalFastRotationPrecompute(ciphertext);
            ciphertext   = cc->EvalFastRotationExt(ciphertext, 1, digits, false);
            ciphertext   = cc->KeySwitchDown(ciphertext);
            elements     = ciphertext->GetElements();
            elements[0] += firstCurrent;
            ciphertext->SetElements(elements);

            cc->Decrypt(keyPair.secretKey, ciphertext, &result);
            result->SetLength(b.size());
            checkEquality(result->GetCKKSPackedValue(), comp_b, eps,
                          failmsg + " Bootstrapping for EvalFastRotationExt w/o first element failed");
        }
        catch (std::exception& e) {
            std::cerr << "Exception thrown from " << __func__ << "(): " << e.what() << std::endl;
            // make it fail
            EXPECT_TRUE(0 == 1) << failmsg;
        }
        catch (...) {
            UNIT_TEST_HANDLE_ALL_EXCEPTIONS;
        }
    }

    void UnitTest_Bootstrap_Iterative(const TEST_CASE_UTCKKSRNSCS_BOOT& testData,
                                      const std::string& failmsg = std::string()) {
        try {
            CryptoContext<Element> cc(UnitTestGenerateContext(testData.params));

            cc->EvalBootstrapSetup(testData.levelBudget, testData.dim1, testData.slots);

            auto keyPair = cc->KeyGen();
            cc->EvalBootstrapKeyGen(keyPair.secretKey, testData.slots);
            cc->EvalAtIndexKeyGen(keyPair.secretKey, {6});
            cc->EvalMultKeyGen(keyPair.secretKey);

            auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters());

            std::vector<std::complex<double>> input(
                Fill({0.111111, 0.222222, 0.333333, 0.444444, 0.555555, 0.666666, 0.777777, 0.888888}, testData.slots));
            size_t encodedLength = input.size();

            Plaintext plaintext = cc->MakeCKKSPackedPlaintext(
                input, 1, cryptoParams->GetCompositeDegree() * (MULT_DEPTH - 1), nullptr, testData.slots);
            auto ciphertext      = cc->Encrypt(keyPair.publicKey, plaintext);
            auto ciphertextAfter = cc->EvalBootstrap(ciphertext);

            Plaintext result;
            cc->Decrypt(keyPair.secretKey, ciphertextAfter, &result);
            result->SetLength(encodedLength);
            uint32_t precision =
                std::floor(CalculateApproximationError(result->GetCKKSPackedValue(), plaintext->GetCKKSPackedValue()));

            // Give buffer for precision to be lower than one measured result.
            const double precisionBuffer = 5;
            precision -= precisionBuffer;

            // Add numIterations as a parameter.
            uint32_t numIterations       = 2;
            auto ciphertextTwoIterations = cc->EvalBootstrap(ciphertext, numIterations, precision);

            Plaintext resultTwoIterations;
            cc->Decrypt(keyPair.secretKey, ciphertextTwoIterations, &resultTwoIterations);
            result->SetLength(encodedLength);
            auto actualResult = resultTwoIterations->GetCKKSPackedValue();
            checkEquality(actualResult, plaintext->GetCKKSPackedValue(), eps,
                          failmsg + " Bootstrapping with " + std::to_string(numIterations) + " iterations failed");
            double precisionMultipleIterations =
                CalculateApproximationError(actualResult, plaintext->GetCKKSPackedValue());

            EXPECT_GE(precisionMultipleIterations + precisionBuffer, numIterations * precision);

            auto temp6 = input;
            std::rotate(temp6.begin(), temp6.begin() + 6, temp6.end());

            auto ciphertext6 = cc->EvalAtIndex(ciphertextAfter, 6);
            Plaintext result6;
            cc->Decrypt(keyPair.secretKey, ciphertext6, &result6);
            result6->SetLength(encodedLength);
            checkEquality(result6->GetCKKSPackedValue(), temp6, eps,
                          failmsg + " EvalAtIndex after Bootstrapping for ciphertexts fails");
        }
        catch (std::exception& e) {
            std::cerr << "Exception thrown from " << __func__ << "(): " << e.what() << std::endl;
            // make it fail
            EXPECT_TRUE(0 == 1) << failmsg;
        }
        catch (...) {
            UNIT_TEST_HANDLE_ALL_EXCEPTIONS;
        }
    }

    void UnitTest_Bootstrap_NumTowers(const TEST_CASE_UTCKKSRNSCS_BOOT& testData,
                                      const std::string& failmsg = std::string()) {
        // This test checks to make sure that we return the original ciphertext if we
        // start with more towers than the number of towers we would end up with by
        // bootstrapping.
        try {
            CryptoContext<Element> cc(UnitTestGenerateContext(testData.params));

            cc->EvalBootstrapSetup(testData.levelBudget, testData.dim1, testData.slots);

            auto keyPair = cc->KeyGen();
            cc->EvalBootstrapKeyGen(keyPair.secretKey, testData.slots);
            cc->EvalAtIndexKeyGen(keyPair.secretKey, {6});
            cc->EvalMultKeyGen(keyPair.secretKey);

            std::vector<std::complex<double>> input(
                Fill({0.111111, 0.222222, 0.333333, 0.444444, 0.555555, 0.666666, 0.777777, 0.888888}, testData.slots));
            size_t encodedLength = input.size();

            // We start with a ciphertext with 0 levels consumed.
            Plaintext plaintext  = cc->MakeCKKSPackedPlaintext(input);
            auto ciphertext      = cc->Encrypt(keyPair.publicKey, plaintext);
            auto ciphertextAfter = cc->EvalBootstrap(ciphertext);

            auto initNumTowers          = ciphertext->GetElements()[0].GetNumOfElements();
            auto bootstrappingNumTowers = ciphertextAfter->GetElements()[0].GetNumOfElements();
            // Check to make sure we don't lose any towers.
            EXPECT_EQ(initNumTowers, bootstrappingNumTowers);

            Plaintext result;
            cc->Decrypt(keyPair.secretKey, ciphertextAfter, &result);
            result->SetLength(encodedLength);
            auto actualResult = result->GetCKKSPackedValue();
            checkEquality(actualResult, plaintext->GetCKKSPackedValue(), eps, failmsg + " Bootstrapping failed");

            auto ciphertextTwoIterations             = cc->EvalBootstrap(ciphertext);
            auto bootstrappingNumTowersTwoIterations = ciphertextTwoIterations->GetElements()[0].GetNumOfElements();
            // Check to make sure we don't lose any towers with double-iteration bootstrapping.
            EXPECT_EQ(initNumTowers, bootstrappingNumTowersTwoIterations);

            Plaintext result2;
            cc->Decrypt(keyPair.secretKey, ciphertextTwoIterations, &result2);
            result->SetLength(encodedLength);
            auto actualResult2 = result2->GetCKKSPackedValue();
            checkEquality(actualResult2, plaintext->GetCKKSPackedValue(), eps,
                          failmsg + " Bootstrapping with two iterations failed");
        }
        catch (std::exception& e) {
            std::cerr << "Exception thrown from " << __func__ << "(): " << e.what() << std::endl;
            // make it fail
            EXPECT_TRUE(0 == 1) << failmsg;
        }
        catch (...) {
            UNIT_TEST_HANDLE_ALL_EXCEPTIONS;
        }
    }

    void UnitTest_Bootstrap_Serialize(const TEST_CASE_UTCKKSRNSCS_BOOT& testData,
                                      const std::string& failmsg = std::string()) {
        try {
            CryptoContextImpl<DCRTPoly>::ClearEvalMultKeys();
            CryptoContextImpl<DCRTPoly>::ClearEvalSumKeys();
            CryptoContextImpl<DCRTPoly>::ClearEvalAutomorphismKeys();
            CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();

            CryptoContext<Element> ccInit(UnitTestGenerateContext(testData.params));
            ccInit->EvalBootstrapSetup(testData.levelBudget, testData.dim1, testData.slots, 0, false);
            ccInit->EvalBootstrapSetup(testData.levelBudget, testData.dim1, testData.slots / 2, 0, false);

            auto keyPairInit = ccInit->KeyGen();
            ccInit->EvalMultKeyGen(keyPairInit.secretKey);
            ccInit->EvalBootstrapKeyGen(keyPairInit.secretKey, testData.slots);
            ccInit->EvalBootstrapKeyGen(keyPairInit.secretKey, testData.slots / 2);
            //==============================================================
            // Serialize all necessary objects
            std::stringstream cc_stream;
            Serial::Serialize(ccInit, cc_stream, SerType::BINARY);

            std::stringstream secretKey_stream;
            Serial::Serialize(keyPairInit.secretKey, secretKey_stream, SerType::BINARY);

            std::stringstream publicKey_stream;
            Serial::Serialize(keyPairInit.publicKey, publicKey_stream, SerType::BINARY);

            std::stringstream automorphismKey_stream;
            CryptoContextImpl<DCRTPoly>::SerializeEvalAutomorphismKey(automorphismKey_stream, SerType::BINARY);

            std::stringstream evalMultKey_stream;
            CryptoContextImpl<DCRTPoly>::SerializeEvalMultKey(evalMultKey_stream, SerType::BINARY);
            //====================================================================================================
            // Removed the serialized objects from the memory
            CryptoContextImpl<DCRTPoly>::ClearEvalMultKeys();
            CryptoContextImpl<DCRTPoly>::ClearEvalSumKeys();
            CryptoContextImpl<DCRTPoly>::ClearEvalAutomorphismKeys();
            CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
            //====================================================================================================
            // Deserialize all necessary objects
            CryptoContext<Element> cc;
            Serial::Deserialize(cc, cc_stream, SerType::BINARY);

            KeyPair<Element> keyPair;
            Serial::Deserialize(keyPair.secretKey, secretKey_stream, SerType::BINARY);
            Serial::Deserialize(keyPair.publicKey, publicKey_stream, SerType::BINARY);
            CryptoContextImpl<DCRTPoly>::DeserializeEvalAutomorphismKey(automorphismKey_stream, SerType::BINARY);
            CryptoContextImpl<DCRTPoly>::DeserializeEvalMultKey(evalMultKey_stream, SerType::BINARY);

            cc->EvalBootstrapPrecompute(testData.slots);
            cc->EvalBootstrapPrecompute(testData.slots / 2);
            //====================================================================================================
            std::vector<std::complex<double>> input(
                Fill({0.111111, 0.222222, 0.333333, 0.444444, 0.555555, 0.666666, 0.777777, 0.888888}, testData.slots));
            size_t encodedLength = input.size();

            Plaintext plaintext1  = cc->MakeCKKSPackedPlaintext(input, 1, MULT_DEPTH - 1, nullptr, testData.slots);
            auto ciphertext1      = cc->Encrypt(keyPair.publicKey, plaintext1);
            auto ciphertext1After = cc->EvalBootstrap(ciphertext1);

            Plaintext result;
            cc->Decrypt(keyPair.secretKey, ciphertext1After, &result);
            result->SetLength(encodedLength);
            plaintext1->SetLength(encodedLength);
            checkEquality(result->GetCKKSPackedValue(), plaintext1->GetCKKSPackedValue(), eps,
                          failmsg + " Bootstrapping for fully packed ciphertexts fails");

            //====================================================================================================
            std::vector<std::complex<double>> input2(
                Fill({0.111111, 0.222222, 0.333333, 0.444444}, testData.slots / 2));
            size_t encodedLength2 = input2.size();

            Plaintext plaintext2  = cc->MakeCKKSPackedPlaintext(input2, 1, MULT_DEPTH - 1, nullptr, testData.slots / 2);
            auto ciphertext2      = cc->Encrypt(keyPair.publicKey, plaintext2);
            auto ciphertext2After = cc->EvalBootstrap(ciphertext2);

            cc->Decrypt(keyPair.secretKey, ciphertext2After, &result);
            result->SetLength(encodedLength2);
            plaintext2->SetLength(encodedLength2);
            checkEquality(result->GetCKKSPackedValue(), plaintext2->GetCKKSPackedValue(), eps,
                          failmsg + " Bootstrapping for fully packed ciphertexts fails");
            //====================================================================================================
            EXPECT_TRUE(1 == 1) << failmsg;
        }
        catch (std::exception& e) {
            std::cerr << "Exception thrown from " << __func__ << "(): " << e.what() << std::endl;
            // make it fail
            EXPECT_TRUE(0 == 1) << failmsg;
        }
        catch (...) {
            UNIT_TEST_HANDLE_ALL_EXCEPTIONS;
        }
    }
};

//===========================================================================================================
TEST_P(UTCKKSRNSCS_BOOT, CKKSRNS) {
    setupSignals();
    auto test = GetParam();

    switch (test.testCaseType) {
        case BOOTSTRAP_FULL:
        case BOOTSTRAP_EDGE:
        case BOOTSTRAP_SPARSE:
            UnitTest_Bootstrap(test, test.buildTestName());
            break;
        case BOOTSTRAP_KEY_SWITCH:
            UnitTest_Bootstrap_KeySwitching(test, test.buildTestName());
            break;
        case BOOTSTRAP_ITERATIVE:
            UnitTest_Bootstrap_Iterative(test, test.buildTestName());
            break;
        case BOOTSTRAP_NUM_TOWERS:
            UnitTest_Bootstrap_NumTowers(test, test.buildTestName());
            break;
        case BOOTSTRAP_SERIALIZE:
            UnitTest_Bootstrap_Serialize(test, test.buildTestName());
            break;
        default:
            break;
    }
}

INSTANTIATE_TEST_SUITE_P(UnitTests, UTCKKSRNSCS_BOOT, ::testing::ValuesIn(testCases), testName);
