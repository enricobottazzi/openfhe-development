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
  Description:
  This code benchmarks multiplication performance between BigInteger and DCRTPoly approaches.
 */

#define _USE_MATH_DEFINES

#include "benchmark/benchmark.h"
#include "lattice/lat-hal.h"
#include "math/discreteuniformgenerator.h"

#include <iostream>
#include <map>
#include <memory>
#include <vector>
#include <chrono>

using namespace lbcrypto;

void BM_ScalarMultiplication_BigIntVsDCRTPoly(benchmark::State& state) {
    usint n      = 8192;  // cyclotomic order
    size_t kRes  = 51;
    size_t size  = 12;

    auto params = std::make_shared<ILDCRTParams<BigInteger>>(2 * n, size, kRes);

    // Create a uniform random polynomial
    DCRTPoly::DugType dug;
    DCRTPoly u(dug, params, Format::EVALUATION);

    // TEST 1: Multiplication with BigInteger 1
    BigInteger one(1);
    DCRTPoly result1;

    for (auto _ : state) {
        auto start = std::chrono::high_resolution_clock::now();
        result1 = u * one;
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        state.counters["BigInteger_Mult_Time"] = elapsed.count();
    }

    // TEST 2: Multiplication with DCRTPoly 1
    std::string value = {"1"};
    BigInteger bigInt(value);
    BigVector bigVec(params->GetRingDimension(), params->GetModulus());
    bigVec[0] = bigInt;
    
    // Create a Poly from the BigInteger
    lbcrypto::PolyImpl<lbcrypto::BigVector> polyLarge(params, Format::COEFFICIENT);
    polyLarge.SetValues(bigVec, Format::COEFFICIENT);

    // Convert polyLarge to a DCRTPoly and switch format
    lbcrypto::DCRTPoly onePoly(polyLarge, params);
    onePoly.SwitchFormat();
    
    DCRTPoly result2;
    
    for (auto _ : state) {
        auto start2 = std::chrono::high_resolution_clock::now();
        result2 = u * onePoly;
        auto end2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed2 = end2 - start2;
        
        state.counters["DCRTPoly_Mult_Time"] = elapsed2.count();
    }
    
    // Check if the results match
    bool resultsMatch = (result1 == result2);
    state.counters["Results_Match"] = resultsMatch ? 1 : 0;
    
    if (!resultsMatch) {
        std::cout << "WARNING: Results from BigInteger and DCRTPoly multiplication do not match!" << std::endl;
    } else {
        std::cout << "SUCCESS: Results from BigInteger and DCRTPoly multiplication match." << std::endl;
    }
}

BENCHMARK(BM_ScalarMultiplication_BigIntVsDCRTPoly)
    ->Unit(benchmark::kMicrosecond)
    ->MinTime(0.5);

// execute the benchmarks
BENCHMARK_MAIN();