#include "openfhecore.h"
#include <iostream>
#include <iomanip>
#include <mach/mach.h>
#include <chrono>  // For timing measurements

// For trapdoor operations
#include "lattice/trapdoor.h"
#include "lattice/dgsampling.h"

using namespace lbcrypto;

// Function to get current memory usage in bytes
size_t getMemoryUsageBytes() {
    task_basic_info info;
    mach_msg_type_number_t infoCount = TASK_BASIC_INFO_COUNT;

    if (task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t)&info, &infoCount) != KERN_SUCCESS)
        return 0;

    return info.resident_size;
}

int main() {
    std::cout << "Memory Consumption Measurement for TrapdoorGenSquareMat and GaussSampSquareMat in OpenFHE" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    
    // Parameters for trapdoor generation
    usint n         = 8192;
    size_t dcrtBits = 51;
    size_t size     = 7;
    double sigma    = 4.578;
    size_t d        = 4;

    // Create parameters
    auto params = std::make_shared<ILDCRTParams<BigInteger>>(2 * n, size, dcrtBits);

    double val    = params->GetModulus().ConvertToDouble();
    double logTwo = std::ceil(log2(val));
    usint k       = (usint)(logTwo);

    auto zero_alloc    = DCRTPoly::Allocator(params, Format::EVALUATION);
    auto uniform_alloc = DCRTPoly::MakeDiscreteUniformAllocator(params, Format::EVALUATION);
    
    // Record memory and time before TrapdoorGenSquareMat
    size_t memBefore = getMemoryUsageBytes();
    auto startTime = std::chrono::high_resolution_clock::now();
    
    std::cout << "Memory before TrapdoorGenSquareMat: " 
              << memBefore << " bytes" << std::endl;
    
    // Perform TrapdoorGenSquareMat operation
    std::pair<Matrix<DCRTPoly>, RLWETrapdoorPair<DCRTPoly>> trapPair =
        RLWETrapdoorUtility<DCRTPoly>::TrapdoorGenSquareMat(params, sigma, d);
    
    // Record memory and time after TrapdoorGenSquareMat
    auto endTime = std::chrono::high_resolution_clock::now();
    size_t memAfter = getMemoryUsageBytes();
    
    // Calculate memory usage and time
    size_t memUsed = memAfter - memBefore;
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
    
    std::cout << "Memory after TrapdoorGenSquareMat: " 
              << memAfter << " bytes" << std::endl;
    std::cout << "Memory difference for TrapdoorGenSquareMat: " 
              << memUsed << " bytes" << std::endl;
    std::cout << "Time for TrapdoorGenSquareMat: " 
              << duration.count() << " seconds" << std::endl;
    
    // Now measure memory for GaussSampSquareMat
    std::cout << std::string(80, '-') << std::endl;
    std::cout << "Measuring memory for GaussSampSquareMat:" << std::endl;
    
    // Set up parameters for GaussSampSquareMat
    uint32_t base = 2;
    double c      = (base + 1) * sigma;
    double s      = SPECTRAL_BOUND_D(n, k, base, d);
    DCRTPoly::DggType dgg(sigma);
    DCRTPoly::DggType dggLargeSigma(sqrt(s * s - c * c));
    
    Matrix<DCRTPoly> U(zero_alloc, d, d, uniform_alloc);
    
    // Record memory and time before GaussSampSquareMat
    size_t memBeforeGauss = getMemoryUsageBytes();
    auto startTimeGauss = std::chrono::high_resolution_clock::now();
    
    std::cout << "Memory before GaussSampSquareMat: " 
              << memBeforeGauss << " bytes" << std::endl;
    
    // Perform GaussSampSquareMat operation
    Matrix<DCRTPoly> z = RLWETrapdoorUtility<DCRTPoly>::GaussSampSquareMat(n, k, trapPair.first,
                                                                           trapPair.second, U, dgg, dggLargeSigma);
    
    // Record memory and time after GaussSampSquareMat
    auto endTimeGauss = std::chrono::high_resolution_clock::now();
    size_t memAfterGauss = getMemoryUsageBytes();
    
    // Calculate memory usage and time
    size_t memUsedGauss = memAfterGauss - memBeforeGauss;
    auto durationGauss = std::chrono::duration_cast<std::chrono::duration<double>>(endTimeGauss - startTimeGauss);
    
    std::cout << "Memory after GaussSampSquareMat: " 
              << memAfterGauss << " bytes" << std::endl;
    std::cout << "Memory difference for GaussSampSquareMat: " 
              << memUsedGauss << " bytes" << std::endl;
    std::cout << "Time for GaussSampSquareMat: " 
              << durationGauss.count() << " seconds" << std::endl;
    
    std::cout << std::string(80, '-') << std::endl;
    std::cout << "Note: Memory is measured using the Mach API (resident_size) on macOS in bytes." << std::endl;
    std::cout << "      Negative values may indicate memory was freed during the operation." << std::endl;
    std::cout << "      Time is measured in seconds using std::chrono::high_resolution_clock." << std::endl;
    
    return 0;
}
