#ifndef SOBOL_SEQUENCE_H
#define SOBOL_SEQUENCE_H

#include <vector>
#include <cstdint>
#include <stdexcept>
#include <iostream>
#include <iomanip>
/**
 * @class SobolSequence
 * @brief Generates multi-dimensional Sobol sequences using Gray code optimization
 * 
 * This class implements the Sobol sequence generation algorithm as described by
 * Sobol (1967) with Gray code optimization for efficient computation.
 * 
 * The implementation supports up to 1111 dimensions using pre-computed primitive
 * polynomials and initial direction numbers.
 */
class SobolSequence {
private:
    // Maximum supported dimensions
    static constexpr uint32_t MAX_DIMENSIONS = 10; // Simplified for this example
    static constexpr uint32_t MAX_BITS = 32;
    
    /**
     * @struct PrimitivePolynomial
     * @brief Stores primitive polynomial coefficients and degree
     */
    struct PrimitivePolynomial {
        uint32_t degree;        // Degree of the polynomial
        uint32_t coefficients;  // Binary coefficients (bit-packed)
        std::vector<uint32_t> initial_m; // Initial direction numbers
    };
    
    /**
     * Pre-computed primitive polynomials for first 10 dimensions
     * Format: {degree, coefficients, {initial_m_values...}}
     * 
     * Coefficients are stored as binary numbers where bit i represents
     * the coefficient of x^i term (excluding the leading x^degree term)
     */
    static const std::vector<PrimitivePolynomial> primitive_polynomials;
    
    uint32_t dimensions_;           // Number of dimensions
    uint64_t current_index_;       // Current sequence index
    std::vector<std::vector<uint32_t>> direction_numbers_; // Direction numbers for each dimension
    std::vector<uint32_t> current_point_;  // Current point (as integers)
    
    /**
     * @brief Initializes direction numbers for all dimensions
     * 
     * Direction numbers are computed using the recurrence relation based on
     * primitive polynomials. The first few direction numbers are given as
     * initial values, and the rest are computed using the recurrence.
     */
    void initializeDirectionNumbers() {
        direction_numbers_.resize(dimensions_);
        
        for (uint32_t dim = 0; dim < dimensions_; ++dim) {
            direction_numbers_[dim].resize(MAX_BITS);
            
            if (dim == 0) {
                // First dimension is special - uses Van der Corput sequence in base 2
                // Direction numbers are powers of 2: 2^31, 2^30, 2^29, ...
                for (uint32_t i = 0; i < MAX_BITS; ++i) {
                    direction_numbers_[dim][i] = 1U << (31 - i);
                }
            } else {
                // For other dimensions, use primitive polynomials
                const auto& poly = primitive_polynomials[dim - 1];
                
                // Set initial direction numbers (given values, left-shifted appropriately)
                // The initial direction numbers are scaled by 2^(32-i)
                // effectively, the recurrence relation that was:
                // m_i = 2*a_1*m_{i-1} XOR 2^2*a_2*m_{i-2} XOR ... XOR 2^{k-1}*a_{k-1}*m_{i-(k-1)} XOR 2^k*m_{i-k} XOR m_{i-k}
                // is now converted to:
                // m_i = a_1*m_{i-1} XOR a_2*m_{i-2} XOR ... XOR m_{i-k} XOR m_{i-k}/2^k
                // k is the degree of the polynomial
                // Therefore, m_i = u_i * 2^(32-i) where u_i is the original direction number
                // we skip the step v_i = m_i / 2^i so that we can efficiently store
                // the direction numbers as integers in [0, 2^32)
                // Once the real points are computed, we scale them back to [0, 1) {see line 167}
                for (uint32_t i = 0; i < poly.initial_m.size() && i < MAX_BITS; ++i) {
                    direction_numbers_[dim][i] = poly.initial_m[i] << (31 - i);
                }
                
                // Compute remaining direction numbers using recurrence relation
                // m_i = a_1*m_{i-1} XOR a_2*m_{i-2} XOR ... XOR m_{i-k} XOR m_{i-k}/2^k
                // k is the degree of the polynomial
                for (uint32_t i = poly.degree; i < MAX_BITS; ++i) {
                    // m_{i-k} XOR m_{i-k}/2^k
                    uint32_t m = direction_numbers_[dim][i - poly.degree] ^ (direction_numbers_[dim][i - poly.degree] >> poly.degree);
                    
                    for (uint32_t j = 1; j < poly.degree; ++j) {
                        if (poly.coefficients & (1U << (j - 1))) {
                            // XOR m_{i-1} ... XOR m_{i-(k-1)}
                            m ^= direction_numbers_[dim][i - j];
                        }
                    }
                    
                    direction_numbers_[dim][i] = m;
                }
            }
        }
    }
public:
    /**
     * @brief Constructor
     * @param dimensions Number of dimensions (must be <= MAX_DIMENSIONS)
     * @throws std::invalid_argument if dimensions is 0 or exceeds maximum
     */
    explicit SobolSequence(uint32_t dimensions) 
        : dimensions_(dimensions), current_index_(0) {
        if (dimensions == 0) {
            throw std::invalid_argument("Dimensions must be at least 1");
        }
        if (dimensions > MAX_DIMENSIONS) {
            throw std::invalid_argument("Too many dimensions requested");
        }
        
        initializeDirectionNumbers();
        current_point_.resize(dimensions_, 0);
    }

    /** 
     * @brief Function to provide Rightmost zero bit
     * @param N Number of points to generate
     * @return Array of uint32_t containing the rightmost zero bit for each index
    */
    uint32_t* rightmostZeroBit(uint32_t N) const {
        // Count trailing ones to find the rightmost zero
        uint32_t *C = new uint32_t[N];
        C[0] = 1;
        for (uint32_t i = 1; i <= N - 1; i++)
        {
            C[i] = 1;
            uint32_t value = i;
            while (value & 1)
            {
                value >>= 1;
                C[i]++;
            }
        }
        return C;
    }
    
    /**
     * @brief Generates the next point in the Sobol sequence
     * @return Vector of doubles in [0, 1) representing the next point
     * 
     * Uses Gray code optimization: instead of computing the entire point from scratch,
     * we XOR the previous point with a single direction number corresponding to
     * the bit that changed in the Gray code representation.
     */
    std::vector<double> nextPoint(uint32_t changed_bit) {
        if (current_index_ == 0) {
            // First point is always the origin
            current_index_++;
            return std::vector<double>(dimensions_, 0.0);
        }
        
        // Update each dimension by XORing with the appropriate direction number
        std::vector<double> point(dimensions_);
        for (uint32_t dim = 0; dim < dimensions_; ++dim) {
            current_point_[dim] ^= direction_numbers_[dim][changed_bit];
            // Convert to double in [0, 1) by dividing by 2^32
            point[dim] = static_cast<double>(current_point_[dim]) / (1ULL << 32);
        }
        
        current_index_++;
        return point;
    }
    
    /**
     * @brief Resets the sequence to the beginning
     */
    void reset() {
        current_index_ = 0;
        std::fill(current_point_.begin(), current_point_.end(), 0);
    }
    
    /**
     * @brief Gets the current index in the sequence
     * @return Current index (0-based)
     */
    uint64_t getCurrentIndex() const {
        return current_index_;
    }
    
    /**
     * @brief Gets the number of dimensions
     * @return Number of dimensions
     */
    uint32_t getDimensions() const {
        return dimensions_;
    }
    
    /**
     * @brief Utility function to print direction numbers (for debugging)
     * @param dim Dimension to print (0-indexed)
     * @param count Number of direction numbers to print
     */
    void printDirectionNumbers(uint32_t dim, uint32_t count = 10) const {
        if (dim >= dimensions_) {
            std::cout << "Invalid dimension" << std::endl;
            return;
        }
        
        std::cout << "Direction numbers for dimension " << dim << ":" << std::endl;
        for (uint32_t i = 0; i < std::min(count, MAX_BITS); ++i) {
            std::cout << "v[" << i << "] = " << std::hex << std::setw(8) << std::setfill('0') 
                      << direction_numbers_[dim][i] << std::dec << std::endl;
        }
    }
};
// Static member definition - primitive polynomials for dimensions 2-11
// Each entry contains: {degree, coefficients, {initial direction numbers}}
const std::vector<SobolSequence::PrimitivePolynomial> SobolSequence::primitive_polynomials = {
    // Dimension 2: x^1 + 1 (coefficients = 0, since no middle terms)
    {1, 0, {1}},
    
    // Dimension 3: x^2 + x + 1 (coefficients = 1, representing the x term)
    {2, 1, {1, 3}},
    
    // Dimension 4: x^3 + x + 1
    {3, 1, {1, 3, 1}},
    
    // Dimension 5: x^3 + x^2 + 1
    {3, 2, {1, 1, 1}},
    
    // Dimension 6: x^4 + x + 1
    {4, 1, {1, 1, 3, 3}},
    
    // Dimension 7: x^4 + x^3 + 1
    {4, 4, {1, 3, 5, 13}},
    
    // Dimension 8: x^5 + x^2 + 1
    {5, 2, {1, 1, 5, 5, 17}},
    
    // Dimension 9: x^5 + x^4 + x^2 + x + 1
    {5, 11, {1, 1, 5, 5, 5}},
    
    // Dimension 10: x^5 + x^4 + x^3 + x + 1
    {5, 13, {1, 3, 15, 17, 63}}
};
#endif // SOBOL_SEQUENCE_H
/**
 * Example usage and testing
 */
#include <iostream>
#include <vector>

// Example: Integrate f(x,y) = x*y over [0,1]²
double monteCarloIntegrate() {
    SobolSequence sobol(2);
    double sum = 0.0;
    int n = 10000;

    uint32_t* C = sobol.rightmostZeroBit(n);
    
    for (int i = 0; i < n; ++i) {
        auto point = sobol.nextPoint(C[i]-1);
        sum += point[0] * point[1]; // f(x,y) = x*y
    }
    
    return sum / n; // Should converge to 0.25
}

// Example: Integrate f(x,y) = x over [0,1]
double monteCarloIntegrate_1D() {
    SobolSequence sobol(1);
    double sum = 0.0;
    int n = 10000;
    
    uint32_t* C = sobol.rightmostZeroBit(n);

    for (int i = 0; i < n; ++i) {
        auto point = sobol.nextPoint(C[i]-1);
        sum += point[0]; // f(x,y) = x
    }
    
    return sum / n; // Should converge to 0.5
}

int main() {
    try {
        // Create a 2D Sobol sequence generator
        SobolSequence sobol(2);
        
        std::cout << "First 16 points of 2D Sobol sequence:" << std::endl;
        std::cout << "Index\tX\t\tY" << std::endl;
        std::cout << "-----\t--------\t--------" << std::endl;

        int n;
        
        n = 16; // Number of points to generate
        uint32_t* C = sobol.rightmostZeroBit(n);
        for (int i = 0; i < n; ++i) {
            auto point = sobol.nextPoint(C[i]);
            std::cout << i << "\t" << std::fixed << std::setprecision(6) 
                      << point[0] << "\t" << point[1] << std::endl;
        }
        
        // Reset and demonstrate 3D sequence
        std::cout << "\nFirst 8 points of 3D Sobol sequence:" << std::endl;
        std::cout << "Index\tX\t\tY\t\tZ" << std::endl;
        std::cout << "-----\t--------\t--------\t--------" << std::endl;
        
        SobolSequence sobol3d(3);
        n = 8; // Number of points to generate
        uint32_t* C3d = sobol3d.rightmostZeroBit(n);
        for (int i = 0; i < n; ++i) {
            auto point = sobol3d.nextPoint(C3d[i]);
            std::cout << i << "\t" << std::fixed << std::setprecision(6) 
                      << point[0] << "\t" << point[1] << "\t" << point[2] << std::endl;
        }
        
        // Print direction numbers for first dimension (for educational purposes)
        std::cout << "\nDirection numbers for dimension 0:" << std::endl;
        sobol.printDirectionNumbers(0, 8);

        double result = monteCarloIntegrate();
        std::cout << "\nMonte Carlo integration result for f(x,y) = x*y over [0,1]²: " 
                  << std::fixed << std::setprecision(6) << result << std::endl;

        double result_1D = monteCarloIntegrate_1D();
        std::cout << "Monte Carlo integration result for f(x) = x over [0,1]: " 
                  << std::fixed << std::setprecision(6) << result_1D << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
