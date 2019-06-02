
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

#ifndef PBRT_MEDIA_EMISSIONHOMOGENEOUS_H
#define PBRT_MEDIA_EMISSIONHOMOGENEOUS_H

// media/emission_homogeneous.h*
#include "medium.h"

namespace pbrt {

// EmissionHomogeneousMedium Declarations
class EmissionHomogeneousMedium : public Medium {
  public:
    // EmissionHomogeneousMedium Public Methods
    EmissionHomogeneousMedium(const Spectrum &sigma_a, const Spectrum &sigma_s, Float g, Float T)
        : sigma_a(sigma_a),
          sigma_s(sigma_s),
          sigma_t(sigma_s + sigma_a),
          emission_ratio(sigma_a / (sigma_s + sigma_a)),
          scattering_ratio(sigma_s / (sigma_s + sigma_a)),
          g(g)
//          T(T)
          {
            // initialize constant emission radiance
            Medium::tempToRGBSpectrum(T, emission_radiance);
//            Float Le[nCIESamples];
//            BlackbodyNormalized(CIE_lambda, nCIESamples, T, Le);
//            Float xyz[3] = {0, 0, 0};
//            for (int i = 0; i < nCIESamples; ++i) {
//                xyz[0] += Le[i] * CIE_X[i];
//                xyz[1] += Le[i] * CIE_Y[i];
//                xyz[2] += Le[i] * CIE_Z[i];
//            }
////            std::cout << "xyz before scale: " << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;
//            Float scale = Float(CIE_lambda[nCIESamples - 1] - CIE_lambda[0]) /
//                          Float(nCIESamples);
//            xyz[0] *= scale;
//            xyz[1] *= scale;
//            xyz[2] *= scale;
////            std::cout << "xyz after scale: " << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;
//
//            emission_radiance = RGBSpectrum::FromXYZ(xyz);
////            std::cout << "raw rgb: " << emission_radiance << std::endl;
//            Float biggest = emission_radiance.maxComponent();
//            if (biggest > 1.){
//                emission_radiance /= biggest;
//            }
////            std::cout << "rgb after normalize: " << emission_radiance << std::endl;
//            emission_radiance.clampNegative();
            std::cout << "T: " << T << " emission_radiance: " << emission_radiance << std::endl;
    }
    Spectrum Tr(const Ray &ray, Sampler &sampler) const;
    Spectrum Sample(const Ray &ray, Sampler &sampler, MemoryArena &arena,
                    MediumInteraction *mi) const;

  private:
    // EmissionHomogeneousMedium Private Data
    const Spectrum sigma_a, sigma_s, sigma_t, emission_ratio, scattering_ratio;
    const Float g;
//    const Float T;
    Spectrum emission_radiance;
};

}  // namespace pbrt

#endif  // PBRT_MEDIA_EMISSIONHOMOGENEOUS_H
