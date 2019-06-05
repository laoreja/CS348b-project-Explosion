
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

#ifndef PBRT_MEDIA_EMISSIONEVDB_H
#define PBRT_MEDIA_EMISSIONEVDB_H

// media/emission_vdb.h*
#include <sys/stat.h>
#include "medium.h"
#include "transform.h"
#include "stats.h"
#include "sampling.h"
#include <openvdb/openvdb.h>

namespace pbrt {

// EmissionEVDBMedium Declarations
class EmissionEVDBMedium : public Medium {
  public:
    // EmissionEVDBMedium Public Methods
    EmissionEVDBMedium(const Spectrum &sigma_a, const Spectrum &sigma_s, Float g,
                      const Transform &mediumToWorld,
                      Float tempScale, Float densityScale,
                      const std::string & vdbpath, Float sampleProb)
        : Medium(sampleProb),
          sigma_a(sigma_a),
          sigma_s(sigma_s),
          emission_ratio(sigma_a/(sigma_a+sigma_s)),
          g(g),
          MediumToWorld(mediumToWorld),
          WorldToMedium(Inverse(mediumToWorld)),
          tempScale(tempScale),
          densityScale(densityScale)
          {
        // Precompute values for Monte Carlo sampling of _EmissionEVDBMedium_
        sigma_t = (sigma_a + sigma_s)[0];
        if (Spectrum(sigma_t) != sigma_a + sigma_s)
            Error(
                "EmissionEVDBMedium requires a spectrally uniform attenuation "
                "coefficient!");

// This checks whether the file exists.
// OpenVDB will report IoErrors on other issues (like not open correctly) even if the file exists
        struct stat buffer;
        CHECK(stat (vdbpath.c_str(), &buffer) == 0) << "IoError: no file exists at vdbpath: " << vdbpath << std::endl;
        std::cout << "vdbpath: " << vdbpath << " sampleProb: " << sampleProb << std::endl;

        openvdb::initialize();
        openvdb::io::File vdbfile(vdbpath);
        vdbfile.open();
        openvdb::GridBase::Ptr baseGrid;

        bool findDensity = false;
        for (openvdb::io::File::NameIterator nameIter = vdbfile.beginName();
             nameIter != vdbfile.endName(); ++nameIter) {
            std::cout << nameIter.gridName() << std::endl;
            if (nameIter.gridName() == "density") {
                baseGrid = vdbfile.readGrid(nameIter.gridName());
                openvdb::CoordBBox bbox = baseGrid->evalActiveVoxelBoundingBox();

                sx = bbox.min()[0];
                sy = bbox.min()[1];
                sz = bbox.min()[2];
                nx = bbox.extents()[0];
                ny = bbox.extents()[1];
                nz = bbox.extents()[2];

                densityGrid = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);
                findDensity = true;
            } else if (nameIter.gridName() == "temperature" || nameIter.gridName() == "heat") {
                std::cout << "find temperature grid with name: " << nameIter.gridName() << std::endl;
                if (nameIter.gridName() == "temperature") useHeat = false;
                else useHeat = true;
                baseGrid = vdbfile.readGrid(nameIter.gridName());
                temperatureGrid = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);
            }
        }
        CHECK(findDensity) << "No density grid in the vdb file at " << vdbpath << std::endl;
        vdbfile.close();
        if (temperatureGrid == nullptr) {
            temperatureGrid = openvdb::FloatGrid::create(0.);
        }

        density = std::unique_ptr<Float[]>(new Float[nx * ny * nz]());
        LeGrid = std::unique_ptr<Spectrum[]>(new Spectrum[nx * ny * nz]());
        std::unique_ptr<Float[]> forDistribution(new Float[nx * ny * nz]());

        openvdb::FloatGrid ::Accessor tempAccessor = temperatureGrid->getAccessor();
        openvdb::Coord xyz;

        Float curDensity;
        Float maxDensity = 0;
        Float maxRawTemp = 0;
        Spectrum temp;
        int i;

        for (openvdb::FloatGrid::ValueOnIter iter = densityGrid->beginValueOn(); iter; ++iter) {
            curDensity = iter.getValue();
            maxDensity = std::max(maxDensity, curDensity);

            xyz = iter.getCoord();
            i = ( (xyz.z() - sz) * ny + xyz.y() - sy) * nx + xyz.x() - sx;
            density[i] = curDensity * densityScale;

            Float curTemp = tempAccessor.getValue(xyz);
            if (curTemp > maxRawTemp) {
                maxRawTemp = curTemp;
            }
            if (curTemp * tempScale < 1.) {
                continue;
            }
//            if(useHeat)
//                Medium::tempToRGBSpectrum(std::max((curTemp * 0.5f + 1.25f), 6500.f) * tempScale, LeGrid[i]);
//            else
            Float realTemp = curTemp * tempScale;
            if (realTemp > 6000.f) {
                realTemp = 6000.f;
            }
            Medium::tempToRGBSpectrum(realTemp, LeGrid[i]);
            forDistribution[i] = LeGrid[i].y();
        }
        std::cout << "maxDensity: " << maxDensity << std::endl;
        std::cout << "maxRawTemp: " << maxRawTemp << std::endl;
        invMaxDensity = 1 / maxDensity;

        distribution.reset(new Distribution3D(forDistribution.get(), nx, ny, nz));
    }

    Float Density(const Point3f &p) const;
    Spectrum EmissionRadiance(const Point3f &p) const;

    Float D(const Point3i &p) const {
        Bounds3i sampleBounds(Point3i(0, 0, 0), Point3i(nx, ny, nz));
        if (!InsideExclusive(p, sampleBounds)) return 0;
        return density[(p.z * ny + p.y) * nx + p.x];
    }
    Spectrum LeFromGrid(const Point3i &p) const {
        Bounds3i sampleBounds(Point3i(0, 0, 0), Point3i(nx, ny, nz));
        if (!InsideExclusive(p, sampleBounds)) return 0;
        return LeGrid[(p.z * ny + p.y) * nx + p.x];
    }

    Spectrum Sample(const Ray &ray, Sampler &sampler, MemoryArena &arena,
                    MediumInteraction *mi) const;
    Spectrum Tr(const Ray &ray, Sampler &sampler) const;

    Spectrum Sample_Li(const Interaction &ref, Sampler &sampler, Vector3f *wi, Float *xPdf, Float *tPdfGivenWi,
                                           VisibilityTester *vis) const;
    Spectrum Sample_t_given_Wi(const Interaction &ref, Sampler &sampler, const Vector3f &wi, Float *xPdf, Float *tPdfGivenWi, Float *tRay) const;
//    Float Pdf_t_given_Wi(const Interaction &ref, Sampler &sampler, const Vector3f &wi) const;

  private:
    // EmissionEVDBMedium Private Data
    openvdb::FloatGrid::Ptr densityGrid = nullptr;
    openvdb::FloatGrid::Ptr temperatureGrid = nullptr;

    std::unique_ptr<Float[]> density;
    std::unique_ptr<Spectrum[]> LeGrid;

    const Spectrum sigma_a, sigma_s, emission_ratio;
    const Float g;
    int sx, sy, sz;
    int nx, ny, nz;
    const Transform WorldToMedium, MediumToWorld;
    Float sigma_t;
    Float invMaxDensity;

    const Float tempScale, densityScale;

    // For light sampling in MIS
    std::unique_ptr<Distribution3D> distribution;
    bool useHeat;
};

}  // namespace pbrt

#endif  // PBRT_MEDIA_EMISSIONEVDB_H
