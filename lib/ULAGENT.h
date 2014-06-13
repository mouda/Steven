/*
  File: ULAGENT.h
  Brief:1.each object,ULAGENT is one machine/node.
        2.ULAGENT only called/managed by ULSA3b

  Author: Steven
  Date: 2012/08/07

  Latest Progress:

  Abbreviation:

*/
#ifndef ULAGENT_H
#define ULAGENT_H
#include<vector>
#include<cstddef>//For keyworkd like NULL
class ULAGENT
{
  //locations are set in readLocation();
  //Their groups are set in setIniStru();
  public:
    ULAGENT();
    ULAGENT(const ULAGENT &);
    void aryConstructor(int inputIndex, float X, float Y); // Triggered and Set Location
    bool existed;
    int nodeIndex; //nodeindex := nodeName
    float locX;
    float locY;
    double dibSQ;
    int* ptrHead;  //In which group, the number is its cluster head's name.
    float power;
    double T2ndtier;

};
#endif
