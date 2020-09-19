#ifndef OPENING_COLLECTOR_H
#define OPENING_COLLECTOR_H

#include "processing.h"

// Collects openings in the file and maps them to their parent elements
// for processing them after applying the Minkowski sum, as that performs
// much less well on non-convex inputs.
struct opening_collector : public execution_context {
	std::list<shape_callback_item> list;
	std::map<IfcUtil::IfcBaseEntity*, IfcUtil::IfcBaseEntity*> opening_to_elem;
	std::multimap<IfcUtil::IfcBaseEntity*, shape_callback_item*> map;

	void init(IfcParse::IfcFile* f);
	opening_collector(IfcParse::IfcFile* f);
	opening_collector(const std::vector<IfcParse::IfcFile*>& f);
	void operator()(shape_callback_item& item);
};

#endif
