
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

#ifndef PBRT_MEDIA_EMISSIONVDB_H
#define PBRT_MEDIA_EMISSIONVDB_H

// media/emission_vdb.h*
//#include <sys/stat.h>
#include "medium.h"
#include "transform.h"
#include "stats.h"
#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>

//#include "spectrum.h"

namespace pbrt {

// EmissionVDBMedium Declarations
class EmissionVDBMedium : public Medium {
  public:
    // EmissionVDBMedium Public Methods
    EmissionVDBMedium(const Spectrum &sigma_a, const Spectrum &sigma_s, Float g,
                      const Transform &mediumToWorld,
                      Float tempScale, const std::string & vdbpath)
        : sigma_a(sigma_a),
          sigma_s(sigma_s),
          emission_ratio(sigma_a/(sigma_a+sigma_s)),
          g(g),
          WorldToMedium(Inverse(mediumToWorld)),
          tempScale(tempScale) {
        // Precompute values for Monte Carlo sampling of _EmissionVDBMedium_
        sigma_t = (sigma_a + sigma_s)[0];
        if (Spectrum(sigma_t) != sigma_a + sigma_s)
            Error(
                "EmissionVDBMedium requires a spectrally uniform attenuation "
                "coefficient!");

// Such check has already been done in the library
//        struct stat buffer;
//        CHECK(stat (vdbpath.c_str(), &buffer) == 0) << "IoError: no file exists at vdbpath: " << vdbpath << std::endl;
        std::cout << "vdbpath: " << vdbpath << std::endl;
        openvdb::initialize();
        openvdb::io::File vdbfile(vdbpath);
        vdbfile.open();
        openvdb::GridBase::Ptr baseGrid;

        bool findDensity = false;
        for (openvdb::io::File::NameIterator nameIter = vdbfile.beginName();
             nameIter != vdbfile.endName(); ++nameIter) {
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
                baseGrid = vdbfile.readGrid(nameIter.gridName());
                temperatureGrid = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);
            }
        }
        CHECK(findDensity) << "No density grid in the vdb file at " << vdbpath << std::endl;
        vdbfile.close();
        if (temperatureGrid == nullptr) {
            temperatureGrid = openvdb::FloatGrid::create(0.);
        }
        LeGrid = openvdb::Vec3fGrid::create(openvdb::Vec3f(0., 0., 0.));
        openvdb::Vec3fGrid::Accessor LeAccessor = LeGrid->getAccessor();
//        openvdb::FloatGrid ::Accessor DensityAccessor = densityGrid->getAccessor();
        openvdb::FloatGrid ::Accessor tempAccessor = temperatureGrid->getAccessor();
        openvdb::Coord xyz;

        Float curDensity;
        Float maxDensity = 0;
        Spectrum temp;

        for (openvdb::FloatGrid::ValueOnIter iter = densityGrid->beginValueOn(); iter; ++iter) {
            curDensity = iter.getValue();
            maxDensity = std::max(maxDensity, curDensity);

            xyz = iter.getCoord();
            Float curTemp = tempAccessor.getValue(xyz);
            if (curTemp * tempScale < 1.) {
                continue;
            }
            Medium::tempToRGBSpectrum(curTemp * tempScale, temp);
            LeAccessor.setValue(xyz, openvdb::Vec3f(temp[0], temp[1], temp[2]));
        }
        std::cout << "maxDensity: " << maxDensity << std::endl;
        invMaxDensity = 1 / maxDensity;

        densitySampler = std::make_shared<openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler>>(openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler>(*densityGrid));
        LeSampler = std::make_shared<openvdb::tools::GridSampler<openvdb::Vec3fGrid , openvdb::tools::BoxSampler>>(openvdb::tools::GridSampler<openvdb::Vec3fGrid , openvdb::tools::BoxSampler>(*LeGrid));


    }

    Float Density(const Point3f &p) const;
//    Float Temperature(const Point3f &p) const;
    Spectrum Le(const Point3f &p) const;

    Spectrum Sample(const Ray &ray, Sampler &sampler, MemoryArena &arena,
                    MediumInteraction *mi) const;
    Spectrum Tr(const Ray &ray, Sampler &sampler) const;

    Spectrum emission_ratio;
  private:
    // EmissionVDBMedium Private Data
//    openvdb::io::File vdbfile;
    openvdb::FloatGrid::Ptr densityGrid = nullptr;
    openvdb::FloatGrid::Ptr temperatureGrid = nullptr;
    openvdb::Vec3fGrid::Ptr LeGrid;

    openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler>::Ptr densitySampler;
    openvdb::tools::GridSampler<openvdb::Vec3fGrid, openvdb::tools::BoxSampler>::Ptr LeSampler;

    const Spectrum sigma_a, sigma_s;
    const Float g;
    int sx, sy, sz;
    int nx, ny, nz;
    const Transform WorldToMedium;
    // Seems no spectrum grid, so have to use the old method
//    std::unique_ptr<Spectrum[]> LeGrid;
    Float sigma_t;
    Float invMaxDensity;

    Float tempScale;
};

}  // namespace pbrt

#endif  // PBRT_MEDIA_EMISSIONVDB_H
