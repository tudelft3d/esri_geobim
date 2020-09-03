#include "opening_collector.h"

opening_collector::opening_collector(IfcParse::IfcFile * f) {
	auto rels = f->instances_by_type("IfcRelVoidsElement");
	if (!rels) {
		return;
	}
	for (auto& rel : *rels) {
		auto be = (IfcUtil::IfcBaseEntity*) ((IfcUtil::IfcBaseEntity*)rel)->get_value<IfcUtil::IfcBaseClass*>("RelatingBuildingElement");
		auto op = (IfcUtil::IfcBaseEntity*) ((IfcUtil::IfcBaseEntity*)rel)->get_value<IfcUtil::IfcBaseClass*>("RelatedOpeningElement");
		opening_to_elem.insert({ op, be });
	}
}

void opening_collector::operator()(shape_callback_item& item) {
	auto opit = opening_to_elem.find(item.src);
	if (opit != opening_to_elem.end()) {
		list.push_back(item);
		map.insert({ opit->second, &list.back() });
	}
}
