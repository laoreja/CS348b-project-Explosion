
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


// media/emission_vdb.cpp*
#include "media/emission_vdb.h"
#include "paramset.h"
#include "sampler.h"
#include "stats.h"
#include "interaction.h"

//#include "spectrum.h"

namespace pbrt {

STAT_RATIO("Media/EmissionOpenVDB steps per Tr() call", nEmissionOpenVDBTrSteps, nEmissionOpenVDBTrCalls);

// EmissionVDBMedium Method Definitions

// The following implementation is slower but definitely more memory efficient
Float EmissionVDBMedium::Density(const Point3f &p) const {
    // Compute voxel coordinates and offsets for _p_
    return densitySampler->isSample(openvdb::Vec3R(p.x * nx - .5f + sx, p.y * ny - .5f + sy, p.z * nz - .5f + sz));
}

Spectrum EmissionVDBMedium::Le(const Point3f &p) const {
    // Compute voxel coordinates and offsets for _p_
    openvdb::Vec3f Le_vec3f = LeSampler->isSample(openvdb::Vec3R(p.x * nx - .5f + sx, p.y * ny - .5f + sy, p.z * nz - .5f + sz));
    return Spectrum(Le_vec3f.x(), Le_vec3f.y(), Le_vec3f.z());
}

Spectrum EmissionVDBMedium::Sample(const Ray &rWorld, Sampler &sampler,
                                   MemoryArena &arena,
                                   MediumInteraction *mi) const {
    ProfilePhase _(Prof::MediumSample);
    Ray ray = WorldToMedium(
        Ray(rWorld.o, Normalize(rWorld.d), rWorld.tMax * rWorld.d.Length()));
    // Compute $[\tmin, \tmax]$ interval of _ray_'s overlap with medium bounds
    const Bounds3f b(Point3f(0, 0, 0), Point3f(1, 1, 1));
    Float tMin, tMax;
    if (!b.IntersectP(ray, &tMin, &tMax)) return Spectrum(1.f);

    // Run delta-tracking iterations to sample a medium interaction
    Float t = tMin;
    while (true) {
        t -= std::log(1 - sampler.Get1D()) * invMaxDensity / sigma_t;
        if (t >= tMax) break;
        if (Density(ray(t)) * invMaxDensity > sampler.Get1D()) {
            // TODO: fix here
            bool sampledEmission = sampler.Get1D() < emission_ratio[0];
            if (sampledEmission){
                // Populate _mi_ with medium interaction information and return
                *mi = MediumInteraction(rWorld(t), -rWorld.d, rWorld.time, this,
                                        nullptr, Le(ray(t)));
                return Spectrum(1.f);
            } else {
                // Populate _mi_ with medium interaction information and return
                PhaseFunction *phase = ARENA_ALLOC(arena, HenyeyGreenstein)(g);
                *mi = MediumInteraction(rWorld(t), -rWorld.d, rWorld.time, this,
                                        phase);
                return Spectrum(1.f);
            }
        }
    }
    return Spectrum(1.f);
}

Spectrum EmissionVDBMedium::Tr(const Ray &rWorld, Sampler &sampler) const {
    ProfilePhase _(Prof::MediumTr);
    ++nEmissionOpenVDBTrCalls;

    Ray ray = WorldToMedium(
        Ray(rWorld.o, Normalize(rWorld.d), rWorld.tMax * rWorld.d.Length()));
    // Compute $[\tmin, \tmax]$ interval of _ray_'s overlap with medium bounds
    const Bounds3f b(Point3f(0, 0, 0), Point3f(1, 1, 1));
    Float tMin, tMax;
    if (!b.IntersectP(ray, &tMin, &tMax)) return Spectrum(1.f);

    // Perform ratio tracking to estimate the transmittance value
    Float Tr = 1, t = tMin;
    while (true) {
        ++nEmissionOpenVDBTrSteps;
        t -= std::log(1 - sampler.Get1D()) * invMaxDensity / sigma_t;
        if (t >= tMax) break;
        Float density = Density(ray(t));
        Tr *= 1 - std::max((Float)0, density * invMaxDensity);
        // Added after book publication: when transmittance gets low,
        // start applying Russian roulette to terminate sampling.
        const Float rrThreshold = .1;
        if (Tr < rrThreshold) {
            Float q = std::max((Float).05, 1 - Tr);
            if (sampler.Get1D() < q) return 0;
            Tr /= 1 - q;
        }
    }
    return Spectrum(Tr);
}

}  // namespace pbrt
