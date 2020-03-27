#include "DatabaseSpecies.h"

DatabaseSpecies::DatabaseSpecies()
{
	ns_ = 0;
	is_active_ = true;
}

void DatabaseSpecies::SetFromXML(const boost::filesystem::path file_name)
{
	// Read database
	boost::property_tree::ptree ptree;
	boost::property_tree::read_xml(file_name.string(), ptree);

	BOOST_FOREACH(boost::property_tree::ptree::value_type const& node, ptree.get_child("database"))
	{
		boost::property_tree::ptree subtree = node.second;

		if (node.first == "species")
		{
			names_.push_back(subtree.get<std::string>("<xmlattr>.name"));
			chem_names_.push_back(subtree.get<std::string>("<xmlattr>.chemName"));
			CAS_.push_back(subtree.get<std::string>("<xmlattr>.CAS"));
		}
	}

	// Number of species
	ns_ = names_.size();
	is_active_ = true;

	// Print on the screen
	Summary();
}

void DatabaseSpecies::Summary()
{
	for (unsigned int i = 0; i < ns_; i++)
		std::cout << names_[i] << " " << chem_names_[i] << " " << CAS_[i] << std::endl;
}
