#ifndef INCLUDED_LBM_SIMULATION
#define INCLUDED_LBM_SIMULATION

#include <memory>

#include "../Domains/DomainInitializer.h"
#include "../VelocitySets/VelocitySet.h"
#include "node.h"
#include "../LBM/Messenger.h"

#include "../Reporting/reporting.h"
#include "../Reporting/MatlabReporter.h"

namespace LBM {
    class Simulation
    {
        using Initializer_Ptr = std::unique_ptr<Domains::DomainInitializer>;

        std::unique_ptr<Domain> d_domain;
        
        // Performance timing accumulators
        double d_stream_time;
        double d_communicate_time;
        double d_collision_time;
        double d_poststream_time;
        size_t d_step_count;

        public:
            Simulation(Initializer_Ptr initializer);
            ~Simulation();
            void step();
            void report(::Reporting::MatlabReporter reporter);
            void report();
            // Get simulation configuration for performance reporting
            double getOmega() const { return d_domain->omega; }
            std::vector<size_t> getDomainSize() const { return d_domain->domain_size; }
            // Get performance timing statistics
            double getStreamTime() const { return d_stream_time; }
            double getCommunicateTime() const { return d_communicate_time; }
            double getCollisionTime() const { return d_collision_time; }
            double getPostStreamTime() const { return d_poststream_time; }
            size_t getStepCount() const { return d_step_count; }
            void resetTiming() {
                d_stream_time = 0.0;
                d_communicate_time = 0.0;
                d_collision_time = 0.0;
                d_poststream_time = 0.0;
                d_step_count = 0;
            }

        private:
            void stream();
            void postStreamProcess();
            void collission();
            void communicate();
            void cacheHaloRecvCounts();

    };
}
#endif