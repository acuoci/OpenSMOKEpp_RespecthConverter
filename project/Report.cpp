#include "Report.h"


Report::Report() {

	fileName_ = {};
	experimentType_ = {};
	warningMessage_ = {};
	errorMessage_ = {};

}

void Report::WriteOnASCII(std::ofstream& fOut)
{
}
