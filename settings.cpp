#include "settings.h"

#include <boost/program_options.hpp>

// Parse the command line settings and do basic initialization
int parse_command_line(geobim_settings& settings, int argc, char ** argv) {
	namespace po = boost::program_options;

	std::string entities;

	typedef po::command_line_parser command_line_parser;
	po::options_description options("Command line options");
	options.add_options()
		("help,h", "display usage information")
		("version,v", "display version information")
		("debug,d", "more verbose log messages")
		("openings,o", "whether to process opening subtractions")
		("timings,t", "print timings after execution")
		//("exact-segmentation,e", "use exact kernel in proximity tree for semantic association (not recommended)")
		("openings-posthoc,O", "whether to process opening subtractions posthoc")
		("minkowski-triangles,T", "force minkowski sum on individual triangles, slow..")
		("entities,e", new po::typed_value<std::string, char>(&entities), "semicolon separated list of IFC entities to include or exclude")
		("exclude,x", "entities are to be excluded")
		("input-file", new po::typed_value<std::string, char>(&settings.input_filename), "input IFC file")
		("output-file", new po::typed_value<std::string, char>(&settings.output_filename), "output OBJ file")
		("radii", boost::program_options::value<std::vector<std::string>>()->multitoken());

	po::positional_options_description positional_options;
	positional_options.add("input-file", 1);
	positional_options.add("output-file", 1);
	positional_options.add("radii", -1);

	po::variables_map vmap;
	try {
		po::store(command_line_parser(argc, argv).
			options(options).positional(positional_options).run(), vmap);
	} catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
	}

	po::notify(vmap);

	if (vmap.count("version")) {
		return 0;
	} else if (vmap.count("help")) {
		return 0;
	} else if (!vmap.count("input-file")) {
		return 1;
	}

	if (!vmap.count("output-file")) {
		settings.output_filename = settings.input_filename;
	}

	// using Kernel_::FT creates weird segfaults, probably due to how ifopsh constructs the box, coordinates components cannot be shared?
	if (vmap.count("radii")) {
		const auto& radii_str = vmap["radii"].as<std::vector<std::string>>();
		std::transform(radii_str.begin(), radii_str.end(), std::back_inserter(settings.radii), [](const std::string& s) {
			// radii.push_back(CGAL::Gmpq(argv[i]));
			return boost::lexical_cast<double>(s);
		});
		std::sort(settings.radii.begin(), settings.radii.end());
	}

	if (settings.radii.empty()) {
		// Binary search will be used.
	}

	settings.apply_openings = vmap.count("openings");
	settings.apply_openings_posthoc = vmap.count("openings-posthoc");
	settings.exact_segmentation = vmap.count("exact-segmentation");
	settings.debug = vmap.count("debug");
	settings.minkowski_triangles = vmap.count("minkowski-triangles");

	settings.file = new IfcParse::IfcFile(settings.input_filename);

	if (!settings.file->good()) {
		std::cerr << "[Error] Unable to parse input file '" << settings.input_filename << "'";
		return 1;
	}

	settings.settings.set(ifcopenshell::geometry::settings::USE_WORLD_COORDS, false);
	settings.settings.set(ifcopenshell::geometry::settings::WELD_VERTICES, false);
	settings.settings.set(ifcopenshell::geometry::settings::SEW_SHELLS, true);
	settings.settings.set(ifcopenshell::geometry::settings::CONVERT_BACK_UNITS, true);
	settings.settings.set(ifcopenshell::geometry::settings::DISABLE_TRIANGULATION, true);
	settings.settings.set(ifcopenshell::geometry::settings::DISABLE_OPENING_SUBTRACTIONS, !settings.apply_openings);

	if (vmap.count("entities")) {
		std::vector<std::string> tokens;
		boost::split(tokens, vmap["entities"].as<std::string>(), boost::is_any_of(";"));
		if (!tokens.empty()) {
			settings.entity_names.emplace();
			settings.entity_names->clear();
			settings.entity_names->insert(tokens.begin(), tokens.end());
		}
		settings.entity_names_included = !vmap.count("exclude");
	}

	return 0;
}
