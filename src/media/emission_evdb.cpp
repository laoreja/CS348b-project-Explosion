
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


// media/emission_evdb.cpp*
#include "media/emission_evdb.h"
#include "paramset.h"
#include "sampler.h"
#include "stats.h"
#include "interaction.h"
#include "light.h"

//#include "spectrum.h"

namespace pbrt {

STAT_RATIO("Media/EmissionOpenVDB steps per Tr() call", nEmissionOpenVDBTrSteps, nEmissionOpenVDBTrCalls);

// EmissionEVDBMedium Method Definitions

// The following implementation is slower but definitely more memory efficient
Float EmissionEVDBMedium::Density(const Point3f &p) const {
    // Compute voxel coordinates and offsets for _p_
    Point3f pSamples(p.x * nx - .5f, p.y * ny - .5f, p.z * nz - .5f);
    Point3i pi = (Point3i)Floor(pSamples);
    Vector3f d = pSamples - (Point3f)pi;

    // Trilinearly interpolate density values to compute local density
    Float d00 = Lerp(d.x, D(pi), D(pi + Vector3i(1, 0, 0)));
    Float d10 = Lerp(d.x, D(pi + Vector3i(0, 1, 0)), D(pi + Vector3i(1, 1, 0)));
    Float d01 = Lerp(d.x, D(pi + Vector3i(0, 0, 1)), D(pi + Vector3i(1, 0, 1)));
    Float d11 = Lerp(d.x, D(pi + Vector3i(0, 1, 1)), D(pi + Vector3i(1, 1, 1)));
    Float d0 = Lerp(d.y, d00, d10);
    Float d1 = Lerp(d.y, d01, d11);
    return Lerp(d.z, d0, d1);
}

Spectrum EmissionEVDBMedium::EmissionRadiance(const Point3f &p) const {
    // Compute voxel coordinates and offsets for _p_
    Point3f pSamples(p.x * nx - .5f, p.y * ny - .5f, p.z * nz - .5f);
    Point3i pi = (Point3i)Floor(pSamples);
    Vector3f d = pSamples - (Point3f)pi;

    // Trilinearly interpolate density values to compute local density
    Spectrum d00 = Lerp(d.x, LeFromGrid(pi), LeFromGrid(pi + Vector3i(1, 0, 0)));
    Spectrum d10 = Lerp(d.x, LeFromGrid(pi + Vector3i(0, 1, 0)), LeFromGrid(pi + Vector3i(1, 1, 0)));
    Spectrum d01 = Lerp(d.x, LeFromGrid(pi + Vector3i(0, 0, 1)), LeFromGrid(pi + Vector3i(1, 0, 1)));
    Spectrum d11 = Lerp(d.x, LeFromGrid(pi + Vector3i(0, 1, 1)), LeFromGrid(pi + Vector3i(1, 1, 1)));
    Spectrum d0 = Lerp(d.y, d00, d10);
    Spectrum d1 = Lerp(d.y, d01, d11);
    return Lerp(d.z, d0, d1);
}

Spectrum EmissionEVDBMedium::Sample_Li(const Interaction &ref, Sampler &sampler, Vector3f *wi, Float *xPdf, Float *tPdfGivenWi,
                                       VisibilityTester *vis) const {
    ProfilePhase _(Prof::LightSample);

    Point3f p_unit_cube;
    Float pdf_unit_cube;
    p_unit_cube = distribution->SampleContinuous(Point3f(sampler.Get1D(), sampler.Get1D(), sampler.Get1D()),
                                                 &pdf_unit_cube);

    Bounds3f worldBound = (MediumToWorld)(Bounds3f(Point3f(0., 0., 0.), Point3f(1., 1., 1.)));
    *xPdf = pdf_unit_cube / worldBound.Volume();

    Interaction lightIt;
    lightIt.p = MediumToWorld(p_unit_cube, Vector3f(0, 0, 0), &lightIt.pError);
    lightIt.n = Normalize(Normal3f(ref.p - lightIt.p));
    Float lightToRefDistSquared = (lightIt.p - ref.p).LengthSquared();
    if (*xPdf == 0 || lightToRefDistSquared == 0) {
        *xPdf = 0;
        *tPdfGivenWi = 0;
        return Spectrum(0.f);
    }
    Float lightToRefDist = (lightIt.p - ref.p).Length();
    *wi = Normalize(lightIt.p - ref.p);
    *vis = VisibilityTester(ref, lightIt);  // p0: ref, p1: lightIt
    *tPdfGivenWi = Tr(Ray(ref.p, *wi, lightToRefDist), sampler)[0];

    return EmissionRadiance(p_unit_cube) / lightToRefDistSquared;


// Diffuse::Sample_Li
//    ProfilePhase _(Prof::LightSample);
//    Interaction pShape = shape->Sample(ref, u, pdf);
//    pShape.mediumInterface = mediumInterface;
//    if (*pdf == 0 || (pShape.p - ref.p).LengthSquared() == 0) {
//        *pdf = 0;
//        return 0.f;
//    }
//    *wi = Normalize(pShape.p - ref.p);
//    *vis = VisibilityTester(ref, pShape);
//    return L(pShape, -*wi);

// disk::sample
//    Interaction Disk::Sample(const Point2f &u, Float *pdf) const {
//        Point2f pd = ConcentricSampleDisk(u);
//        Point3f pObj(pd.x * radius, pd.y * radius, height);
//        Interaction it;
//        it.n = Normalize((*ObjectToWorld)(Normal3f(0, 0, 1)));
//        if (reverseOrientation) it.n *= -1;
//        it.p = (*ObjectToWorld)(pObj, Vector3f(0, 0, 0), &it.pError);
//        *pdf = 1 / Area();
//        return it;
//    }
}

//Float EmissionEVDBMedium::Pdf_t_given_Wi(const Interaction &ref, Sampler &sampler, const Vector3f &wi) const {
//    ProfilePhase _(Prof::LightPdf);
//    Ray ray = WorldToMedium(Ray(ref.p, wi));
//    // Compute $[\tmin, \tmax]$ interval of _ray_'s overlap with medium bounds
//    const Bounds3f b(Point3f(0, 0, 0), Point3f(1, 1, 1));
//    Float tMin, tMax;
//    if (!b.IntersectP(ray, &tMin, &tMax)) {
//        return 0.f;
//    }
//
//    Float t = tMin;
//    Float pdf = 1.;
//    while (true) {
//        t -= std::log(1 - sampler.Get1D()) * invMaxDensity / sigma_t;
//        if (t >= tMax) {
//            pdf = 0.;
//            break;
//        }
//        Float density = Density(ray(t));
//        pdf *= 1 - std::max((Float)0, density * invMaxDensity);
//        if (density * invMaxDensity > sampler.Get1D()) {
//            break;
//        }
//    }
//    return pdf;
//}

Spectrum EmissionEVDBMedium::Sample_t_given_Wi(const Interaction &ref, Sampler &sampler, const Vector3f &wi, Float *xPdf, Float *tPdfGivenWi, Float *tRay) const {
    ProfilePhase _(Prof::LightPdf);
    Ray ray = WorldToMedium(Ray(ref.p, wi));
    // Compute $[\tmin, \tmax]$ interval of _ray_'s overlap with medium bounds
    const Bounds3f b(Point3f(0, 0, 0), Point3f(1, 1, 1));
    Bounds3f worldBound = (MediumToWorld)(b);

    Float tMin, tMax;
    if (!b.IntersectP(ray, &tMin, &tMax)) {
        *xPdf = 0.;
        *tPdfGivenWi = 0.;
        return Spectrum(0.f);
    }

    Float t = tMin;
    *xPdf = 0;
    *tPdfGivenWi = 1.;
    Spectrum res(0.f);
    while (true) {
        t -= std::log(1 - sampler.Get1D()) * invMaxDensity / sigma_t;
        if (t >= tMax) {
            *tPdfGivenWi = 0.;
            break;
        }
        Float density = Density(ray(t));
        *tPdfGivenWi *= 1 - std::max((Float)0, density * invMaxDensity);
        if (density * invMaxDensity > sampler.Get1D()) {
            Point3f x = ray(t);
            res = EmissionRadiance(x);
            *xPdf = distribution->Pdf(x) / worldBound.Volume();
            *tRay = t;
            break;
        }
    }
    return res;
}

Spectrum EmissionEVDBMedium::Sample(const Ray &rWorld, Sampler &sampler,
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
                                        nullptr, EmissionRadiance(ray(t)));
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

Spectrum EmissionEVDBMedium::Tr(const Ray &rWorld, Sampler &sampler) const {
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
