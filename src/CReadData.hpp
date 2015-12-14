#ifndef CREADDATA_HPP
#define CREADDATA_HPP

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

class CReadData
{
private:
    std::string fileName;
    std::fstream file;
    std::string line;
    std::string valueStr;
    char delimiter;
    double dataValue;
    int lineNum;
//    std::vector< std::vector<int> > data;
    std::vector<double> data;

    void readLine(std::stringstream &lineStream);
public:
    CReadData(std::string fileName, char delimiter);
    ~CReadData();

    void readAllLines();
    void printLines(int numLines);
//    std::vector< std::vector<int> >* getData();
    const std::vector<double>& getData() const;
};

#endif
