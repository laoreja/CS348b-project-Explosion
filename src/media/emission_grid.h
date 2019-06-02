
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MEDIA_EMISSIONGRID_H
#define PBRT_MEDIA_EMISSIONGRID_H

// media/emission_grid.h*
#include "medium.h"
#include "transform.h"
#include "stats.h"
//#include "spectrum.h"

namespace pbrt {

STAT_MEMORY_COUNTER("Memory/Volume density grid", densityEmissionBytes);
STAT_MEMORY_COUNTER("Memory/Volume temperature grid", temperatureBytes);

//static const Matrix3x3 Ma (0.4002400,  0.7076000, -0.0808100,
//                    -0.2263000, 1.1653200,  0.0457000,
//                    0.0000000,  0.0000000,  0.9182200);
//
//static const Matrix3x3 invMa  (1.8599364, -1.1293816, 0.2198974,
//                        0.3611914, 0.6388125, -0.0000064,
//                        0.0000000, 0.0000000,  1.0890636);

//static Vector3f tempToXYZ(Float single_T){
//    // xyz initialized to (0, 0, 0) by vector
//    Vector3f xyz(0., 0., 0.);
//    Float Le[nCIESamples];
//    BlackbodyNormalized(CIE_lambda, nCIESamples, single_T, Le);
////    for(int i = 0; i < nCIESamples; i++){
////        std::cout << Le[i] << " ";
////    }
////    std::cout << std::endl;
//
//    for (int i = 0; i < nCIESamples; ++i) {
////        Float val = InterpolateSpectrumSamples(CIE_lambda, Le, nCIESamples, CIE_lambda[i]);
//        Float val = Le[i];
//        xyz[0] += val * CIE_X[i];
//        xyz[1] += val * CIE_Y[i];
//        xyz[2] += val * CIE_Z[i];
//    }
//    Float scale = Float(CIE_lambda[nCIESamples - 1] - CIE_lambda[0]) /
//                  Float(CIE_Y_integral * nCIESamples);
//    xyz[0] *= scale;
//    xyz[1] *= scale;
//    xyz[2] *= scale;
//    return xyz;
//}

// EmissionGridMedium Declarations
class EmissionGridMedium : public Medium {
  public:
    // EmissionGridMedium Public Methods
    EmissionGridMedium(const Spectrum &sigma_a, const Spectrum &sigma_s, Float g,
                      int nx, int ny, int nz, const Transform &mediumToWorld,
                      const Float *d, const Float *t, Float tempScale)
        : sigma_a(sigma_a),
          sigma_s(sigma_s),
          emission_ratio(sigma_a/(sigma_a+sigma_s)),
          g(g),
          nx(nx),
          ny(ny),
          nz(nz),
          WorldToMedium(Inverse(mediumToWorld)),
          density(new Float[nx * ny * nz]),
          temperature(new Float[nx * ny * nz]),
          tempScale(tempScale) {
        densityEmissionBytes += nx * ny * nz * sizeof(Float);
        temperatureBytes += nx * ny * nz * sizeof(Float);
        memcpy((Float *)density.get(), d, sizeof(Float) * nx * ny * nz);
        memcpy((Float *)temperature.get(), t, sizeof(Float) * nx * ny * nz);
        // TODO: if temperature is 0, then no emission
        // TODO: after things work, try sigmas as a function of density & temperature

        // Precompute values for Monte Carlo sampling of _EmissionGridMedium_
        sigma_t = (sigma_a + sigma_s)[0];
        if (Spectrum(sigma_t) != sigma_a + sigma_s)
            Error(
                "EmissionGridMedium requires a spectrally uniform attenuation "
                "coefficient!");
        Float maxDensity = 0;
        for (int i = 0; i < nx * ny * nz; ++i)
            maxDensity = std::max(maxDensity, density[i]);
        std::cout << "maxDensity: " << maxDensity << std::endl;
        invMaxDensity = 1 / maxDensity;
//        Float maxTemperature = 0;
//        for (int i = 0; i < nx * ny * nz; ++i)
//            maxTemperature = std::max(maxTemperature, temperature[i]);
//        std::cout << "max temperature: " << maxTemperature << " tempScale: " << tempScale << std::endl;
//        Vector3f maxTempXYZ = tempToXYZ(maxTemperature * tempScale);
//        std::cout << "maxTempXYZ: " << maxTempXYZ << std::endl;
//        Vector3f whitePointLMS = Matrix3x3::Mat3x3DotVec(Ma, maxTempXYZ);
//        std::cout << "whitePointLMS: " << whitePointLMS << std::endl;
//        CHECK_NE(whitePointLMS.x, 0);
//        CHECK_NE(whitePointLMS.y, 0);
//        CHECK_NE(whitePointLMS.z, 0);
//        Matrix3x3 invWPLMS(1./whitePointLMS.x, 0., 0.,
//                           0., 1./whitePointLMS.y, 0.,
//                           0., 0., 1./whitePointLMS.z);
//        colorAdaptMat = Matrix3x3::Mul(invWPLMS, Ma);
//        colorAdaptMat = Matrix3x3::Mul(invMa, colorAdaptMat);
//        std::cout << "colorAdaptMat: " << colorAdaptMat << std::endl;

        LeGrid = std::unique_ptr<Spectrum[]>(new Spectrum[nx * ny * nz]);
        for (int i = 0; i < nx * ny * nz; ++i) {
            if (density[i] == 0 || temperature[i] * tempScale < 1.) {
                LeGrid[i] = Spectrum(0);
            } else {
                Medium::tempToRGBSpectrum(temperature[i] * tempScale, LeGrid[i]);
            }
        }
    }

    Float Density(const Point3f &p) const;
    Float Temperature(const Point3f &p) const;
    Spectrum Le(const Point3f &p) const;

    Float D(const Point3i &p) const {
        Bounds3i sampleBounds(Point3i(0, 0, 0), Point3i(nx, ny, nz));
        if (!InsideExclusive(p, sampleBounds)) return 0;
        return density[(p.z * ny + p.y) * nx + p.x];
    }
    Float T(const Point3i &p) const {
        Bounds3i sampleBounds(Point3i(0, 0, 0), Point3i(nx, ny, nz));
        if (!InsideExclusive(p, sampleBounds)) return 0;
        return temperature[(p.z * ny + p.y) * nx + p.x];
    }
    Spectrum LeFromGrid(const Point3i &p) const {
        Bounds3i sampleBounds(Point3i(0, 0, 0), Point3i(nx, ny, nz));
        if (!InsideExclusive(p, sampleBounds)) return 0;
        return LeGrid[(p.z * ny + p.y) * nx + p.x];
    }

//    RGBSpectrum tempToRGBSpectrum(Float single_T);

    Spectrum Sample(const Ray &ray, Sampler &sampler, MemoryArena &arena,
                    MediumInteraction *mi) const;
    Spectrum Tr(const Ray &ray, Sampler &sampler) const;


  private:
    // EmissionGridMedium Private Data
    const Spectrum sigma_a, sigma_s, emission_ratio;
    const Float g;
    const int nx, ny, nz;
    const Transform WorldToMedium;
    std::unique_ptr<Float[]> density;
    std::unique_ptr<Float[]> temperature;
    std::unique_ptr<Spectrum[]> LeGrid;
    Float sigma_t;
    Float invMaxDensity;

    Float tempScale;
//    Matrix3x3 colorAdaptMat;
};

}  // namespace pbrt

#endif  // PBRT_MEDIA_EMISSIONGRID_H
