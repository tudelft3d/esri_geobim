#ifndef SETTINGS_H
#define SETTINGS_H

#include <ifcgeom/settings.h>
#include <ifcparse/IfcFile.h>

#include <boost/optional.hpp>
#include <boost/program_options.hpp>

#include <string>
#include <vector>

// Global settings that are derived from the command line invocation
// parsed by Boost.ProgramOptions
struct geobim_settings {
	std::string input_filename, output_filename;
	std::vector<double> radii;
	bool apply_openings, apply_openings_posthoc, debug, exact_segmentation, minkowski_triangles;
	ifcopenshell::geometry::settings settings;
	boost::optional<std::set<std::string>> entity_names;
	bool entity_names_included;
	IfcParse::IfcFile* file;
};

// Parse the command line settings and do basic initialization
int parse_command_line(geobim_settings& settings, int argc, char** argv);

#endif
