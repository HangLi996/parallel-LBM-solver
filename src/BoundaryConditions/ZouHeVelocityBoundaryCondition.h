#ifndef INCLUDED_BOUNDARY_CONDITIONS_ZOU_HE_VELOCITY_BOUDNARY
#define INCLUDED_BOUNDARY_CONDITIONS_ZOU_HE_VELOCITY_BOUDNARY

#include "PostProcessor.h"
#include "../LBM/node.h"

// Forward declaration
struct Domain;

namespace BoundaryConditions {
    class ZouHeVelocityBoundary : public PostProcessor
    {
        protected:
            std::vector<size_t> d_acts_on_indices;
            std::vector<double> d_velocity;
            Domain *d_domain;

        public:
            ZouHeVelocityBoundary(std::vector<double> velocity, std::vector<size_t> acts_on_indices, Domain *domain);
    };

    class ZouHeVelocityNorthBoundary : public ZouHeVelocityBoundary
    {
        public:
            ZouHeVelocityNorthBoundary(std::vector<double> velocity, std::vector<size_t> acts_on_indices, Domain *domain);
            void process();
    };

    // class ZouHeVelocityEastBoundary : public ZouHeVelocityBoundary
    // {
    //     public:
    //         void process() override;
    // };

    // class ZouHeVelocitySouthBoundary : public ZouHeVelocityBoundary
    // {
    //     public:
    //         void process() override;
    // };

    // class ZouHeVelocityWestBoundary : public ZouHeVelocityBoundary
    // {
    //     public:
    //         void process() override;
    // };

    // Corner
    // class ZouHeVelocityNorthEastCornerBoundary : public ZouHeVelocityBoundary {};
    // class ZouHeVelocityNorthWestCornerBoundary : public ZouHeVelocityBoundary {};

    // class ZouHeVelocitySouthEastCornerBoundary : public ZouHeVelocityBoundary {};
    // class ZouHeVelocitySouthWestCornerBoundary : public ZouHeVelocityBoundary {};
}

#endif