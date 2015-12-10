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
        data.push_back(std::vector<int>());
        
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
        data[lineNum].push_back(dataValue);
    }
}

void CReadData::printLines(int numLines)
{
    printf("Data: \n");
    for (int i = 0; i < numLines; i++)
    {
        printf("line %d: ", i);
        for (int j = 0; j < data[i].size(); j++)
        {
            printf(" %d ", data[i][j]);
        }
        printf("\n");
    }
}

std::vector< std::vector<int> >* CReadData::getData()
{
    return &data;
}

