#include "CReadData.hpp"

CReadData::CReadData(std::string fileName,
                    char delimiter) :
                        fileName(fileName),
                        delimiter(delimiter)
{
    file.open(fileName.c_str());
}

CReadData::~CReadData()
{
}

void CReadData::readAllLines()
{
    lineNum = 0;
    
//    printf("reading file: %s ...\n", fileName.c_str());
    while(std::getline(file, line))
    {
        std::stringstream lineStream(line);
        //data.push_back(std::vector<int>());
        
        readLine(lineStream);
 
        lineNum++;
    }
}

void CReadData::readLine(std::stringstream &lineStream)
{
    while (std::getline(lineStream, valueStr, delimiter))
    {
        std::stringstream ss(valueStr);
        ss >> dataValue;
        //data[lineNum].push_back(dataValue);
        data.push_back(dataValue);
    }
}

void CReadData::printLines(int numLines)
{
    printf("Data: (size: %lu)\n", data.size());
    for (int i = 0; i < data.size(); i++)
    {
        printf(" %d", (int)data[i]);
    }
    printf("\n");
}

//std::vector< std::vector<int> >* CReadData::getData()
const std::vector<double>& CReadData::getData() const
{
    return data;
    //return &data;
}

