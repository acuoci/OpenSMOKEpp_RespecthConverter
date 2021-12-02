#include "Report.h"

Report::Report()
{
	std::vector<std::string> fileName_ = {};
	std::vector<std::string> warningMessage_ = {};

}

void Report::UpdateErrors(std::string type, std::string error) {

	if (type == "WarningMessage")
		warningMessage_.push_back(error);
	else if(type == "GenericWarning")
		genericMessage_.push_back(error);
	else if (type == "ErrorMessage")
		errorMessage_.push_back(error);

}

void Report::UpdateFileName(std::string name) {
	
	fileName_.push_back(name);

}

void Report::WriterReport(const boost::filesystem::path file_name, const boost::filesystem::path output_folder)
{
	std::cout << " * Writing report for the generation of OpenSMOKE++ file(s)..." << std::endl;
	std::ofstream fOut((output_folder / file_name).string(), std::ios::out);
	fOut.setf(std::ios::scientific);
	
	WriteReportOnASCII(fOut);
	
	fOut.close();
}

void Report::WriteReportOnASCII(std::ofstream& fOut)
{
	//if (warningMessage_.size() != 0) {
	std::string ciao = "CUCCIOLO ZINGARELLO";
	//UpdateErrors("WarningMessage", ciao);
	for (int i = 0; i < errorMessage_.size(); i++)
		fOut << errorMessage_[i] << std::endl;
	//}
	//else
		//fOut << "Zinagrello pasticcione" << std::endl;
}
