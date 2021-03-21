/**
 * @file
 * @brief Implementation of Geant4 logging destination
 * @copyright Copyright (c) 2021 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "G4LoggingDestination.hpp"

using namespace allpix;

G4LoggingDestination* G4LoggingDestination::instance = nullptr;
allpix::LogLevel G4LoggingDestination::reporting_level_g4cout = allpix::LogLevel::TRACE;
allpix::LogLevel G4LoggingDestination::reporting_level_g4cerr = allpix::LogLevel::DEBUG;

void G4LoggingDestination::setG4coutReportingLevel(LogLevel level) {
    G4LoggingDestination::reporting_level_g4cout = level;
}

void G4LoggingDestination::setG4cerrReportingLevel(LogLevel level) {
    G4LoggingDestination::reporting_level_g4cerr = level;
}

G4int G4LoggingDestination::ReceiveG4cout(const G4String& msg) {
    if(!msg.empty() && G4LoggingDestination::reporting_level_g4cout <= allpix::Log::getReportingLevel() &&
       !allpix::Log::getStreams().empty()) {
        // Remove line-break always added to G4String
        const_cast<G4String&>(msg).pop_back();
        allpix::Log().getStream(G4LoggingDestination::reporting_level_g4cout,
                                __FILE_NAME__,
                                std::string(static_cast<const char*>(__func__)),
                                __LINE__)
            << msg;
    }
    return 0;
}

G4int G4LoggingDestination::ReceiveG4cerr(const G4String& msg) {
    if(!msg.empty() && G4LoggingDestination::reporting_level_g4cerr <= allpix::Log::getReportingLevel() &&
       !allpix::Log::getStreams().empty()) {
        // Remove line-break always added to G4String
        const_cast<G4String&>(msg).pop_back();
        allpix::Log().getStream(G4LoggingDestination::reporting_level_g4cerr,
                                __FILE_NAME__,
                                std::string(static_cast<const char*>(__func__)),
                                __LINE__)
            << msg;
    }
    return 0;
}
