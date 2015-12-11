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
    
    printf("reading file: %s ...\n", fileName.c_str());
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
    for (int i = 0; i < 11; i++)
    //for (int i = 0; i < data.size(); i++)
    {
        //printf(" %d ", data[i]);
        printf("line %d: ", i);
        //for (int j = 0; j < data[i].size(); j++)
        for (int j = 0; j < 17; j++)
        {
            printf(" %d ", data[i * 17 + j]);
        }
        printf("\n");
    }
    printf("\n");
}

//std::vector< std::vector<int> >* CReadData::getData()
const std::vector<int>& CReadData::getData() const
{
    return data;
    //return &data;
}

